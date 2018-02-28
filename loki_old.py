#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.
#       Author: Francesco Grigoli


import os, sys
import numpy as num
import LatLongUTMconversion
import location
import tt_processing
from obspy.core   import read
from obspy.signal import filter
from scipy.signal import hilbert


class Loki:

    def __init__(self, data_path, output_path, db_path, hdr_filename):
        self.data_path=data_path
        self.output_path=output_path
        self.db_path=db_path
        self.hdr_filename=hdr_filename
        self.data_tree, self.events=self.data_struct(self.data_path, self.output_path)

    def data_struct(self, data_path, output_path):
        events=[]
        data_tree=[]
        for root,dirs,files in os.walk(data_path):
           if not dirs:
              data_tree.append(root)
              events.append(root.split('/')[-1])
        for event in events:
           if not os.path.isdir(output_path):
              os.mkdir(output_path)
        return data_tree, events



#load waveforms
#load Traveltimes
#prepare data to loki format
#process data




    def catalogue_creation(self, event, lat0, lon0, ntrial, refell=23):
        (zorig,eorig,norig)=LatLongUTMconversion.LLtoUTM(refell, lat0, lon0) #da adeguare in python 3
        ev_file=self.output_path+'/'+event+'/'+event+'.loc'
        data=num.loadtxt(ev_file)
        if (ntrial>1):
           w=num.sum(data[:,4])
           xb= ((num.dot(data[:,1],data[:,4])/w)*1000)+eorig
           yb= ((num.dot(data[:,2],data[:,4])/w)*1000)+norig
           late,lone=LatLongUTMconversion.UTMtoLL(refell, yb, xb, zorig)
           zb= num.dot(data[:,3],data[:,4])/w
           cb= num.mean(data[:,4])
           cmax=num.max(data[:,4])
           merr=num.vstack((data[:,1],data[:,2],data[:,3]))
           err=num.cov(merr)
           errmax= num.sqrt(num.max(num.linalg.eigvals(err)))
        else:
           late,lone=LatLongUTMconversion.UTMtoLL(refell, (data[2]*1000)+norig, (data[1]*1000)+eorig, zorig)
           zb= data[3]; errmax='NA'; cb=data[4]; cmax=data[4];
        f=open(self.output_path+'/'+'catalogue','a')
        f.write(event+'    '+str(late)+'   '+str(lone)+'   '+str(zb)+'   '+str(errmax)+'   '+str(cb)+'   '+str(cmax)+'\n')
        f.close()


    def location(self, extension='*', comp=['E','N','Z'], precision='single', *input):
        nshortmin=input[0]; nshortmax=input[1]; slrat=input[2]
        npr=input[3]
        ntrial=input[4]
        traveldb=traveltimes(self.db_path, self.hdr_filename)
        tp=traveldb.load_traveltimes('P', precision)
        ts=traveldb.load_traveltimes('S', precision)
        for event_path in self.data_tree:
            loc=waveforms(event_path, extension, comp)
            event=loc.evid
            print('accessing to the event folder: ', event_path, event)
            if os.path.isdir(self.output_path+'/'+event):
               continue
            else:
               os.mkdir(self.output_path+'/'+event)
            loc.process_data(traveldb.db_stations, epsilon=0.001)
            tp_mod, ts_mod=self.time_extractor(tp, ts, loc.data_stations, traveldb.db_stations)
            tp_mod, ts_mod=tt_processing.tt_f2i(loc.deltat,tp_mod,ts_mod, npr)
            for i in range(ntrial):
                gamma=num.random.random_sample()
                nshort=num.round(nshortmin(1.-gamma*nshortmax)
                nlong=nshort*slrat
                obs_dataP, obs_dataS=loc.recstalta(nshort, nlong) #da cambiare
                corrmatrix=location.stacking(tp_mod, ts_mod, obs_dataP, obs_dataS, npr)
                cmax=num.max(corrmatrix)
                corrmatrix=num.reshape(corrmatrix,(traveldb.nx,traveldb.ny,traveldb.nz))
                (ixloc,iyloc,izloc)=num.unravel_index(num.argmax(corrmatrix),(traveldb.nx,traveldb.ny,traveldb.nz))
                #xloc=ixloc*traveldb.dx; yloc=iyloc*traveldb.dy; zloc=izloc*traveldb.dz
                xloc=traveldb.x[ixloc]; yloc=traveldb.y[iyloc]; zloc=traveldb.x[izloc]
                out_file = open(self.output_path+'/'+event+'/'+event+'.loc','a')
                out_file.write(str(i)+' '+str(xloc)+' '+str(yloc)+' '+str(zloc)+' '+str(cmax)+' '+str(nshort)+' '+str(nlong)+'\n')
                out_file.close()
                num.save(self.output_path+'/'+event+'/'+'corrmatrix_trial_'+str(ntrial),corrmatrix)
                self.coherence_plot(self.output_path+'/'+event, corrmatrix, traveldb.x, traveldb.y, traveldb.z, i)
            self.catalogue_creation(event, traveldb.lat0, traveldb.lon0, ntrial)
        print('Ho finito!!!')


    def detection(self, extension='*', comp=['E','N','Z'], precision='single', *inputs):
        nshort_p=inputs[0]; nshort_s=inputs[1]; slrat=inputs[2]
        npr=inputs[3]
        tstep=inputs[4]
        toverlap=inputs[5]
        traveldb=traveltimes(self.db_path, self.hdr_filename)
        tp=traveldb.load_traveltimes('P', precision)
        ts=traveldb.load_traveltimes('S', precision)
        for event_path in self.data_tree:
            print(event_path)
            loc=waveforms(event_path, extension, comp)
            #event=loc.evid
            event=datetime.datetime.strptime(loc.evid,"%Y-%m-%dT%H:%M:%S.%fZ")
            tp_mod, ts_mod=self.time_extractor(tp, ts, loc.data_stations, traveldb.db_stations)
            loc.process_data(traveldb.db_stations, epsilon=0.001)
            obs_dataP, obs_dataS=loc.recstalta(nshort_p, nshort_s, slrat)
            nsamples=num.size(obs_dataP[0,:])
            noverlap=int(toverlap/loc.deltat)#200
            step=int(tstep/loc.deltat)#500
            tp_mod, ts_mod=tt_processing.tt_f2i(loc.deltat,tp_mod,ts_mod, npr)
            if not os.path.isdir(self.output_path):
               os.mkdir(self.output_path)
            kl=0
            #t0i=int(0*3600*50)
            #t0f=int(3*3600*50)
            detections=[[event,0,0,0,0]]
            for k in range(0,nsamples,step):
                initwin=k-noverlap*kl
                kl=1
                endwin=initwin+step+noverlap
                pstalta=num.zeros([loc.nstation,step+noverlap])
                sstalta=num.zeros([loc.nstation,step+noverlap])
                for i in range(loc.nstation):
                    #normfactorP=num.trapz(obs_dataP[i,initwin:endwin], dx=loc.deltat)
                    #normfactorS=num.trapz(obs_dataS[i,initwin:endwin], dx=loc.deltat)
                    normfactorP=num.max(obs_dataP[i,initwin:endwin])
                    normfactorS=num.max(obs_dataS[i,initwin:endwin])
                    sstalta[i,:]=(obs_dataS[i,initwin:endwin]/normfactorS)
                    pstalta[i,:]=(obs_dataP[i,initwin:endwin]/normfactorP)
                print('range%d-%d'%(initwin,endwin))
                #skip=loc.skip_chunk(pstalta, 5., 3)
                #print('skip ', skip)
                #if skip<15:
                #    continue
                #loc=waveforms(event_path, extension, comp)
                #event=datetime.datetime.strptime(loc.evid,"%Y-%m-%dT%H:%M:%S.%fZ")
                print('accessing to the event folder: ', event_path)
                tini=event+datetime.timedelta(seconds=initwin*loc.deltat)
                tfin=event+datetime.timedelta(seconds=(initwin+step)*loc.deltat)
                print(tini,tfin)
                tdelta=(tfin-tini).total_seconds()
                corrmatrix=ps_detection.stacking(tp_mod, ts_mod, pstalta, sstalta, step, npr)
                #corrmatrix=p_detection.stacking(tp_mod, pstalta, step, npr)
                maxloc=num.max(corrmatrix,axis=0)
                avgloc=num.mean(corrmatrix,axis=0)
                thres=0.0; coherthres=5.0;
                coheratio=maxloc/avgloc
                imaxratio=num.argmax(maxloc)
                print('coherence level', maxloc[imaxratio], 'coherence ratio', coheratio[imaxratio])
                if maxloc[imaxratio]>thres and coheratio[imaxratio]>coherthres:
                   #index=num.argmax(corrmatrix)
                   index=num.argmax(corrmatrix[:,imaxratio])
                   #itime=num.argmax(detections)
                   itime=imaxratio
                   (ixloc,iyloc,izloc)=num.unravel_index(index,(traveldb.nx,traveldb.ny,traveldb.nz))
                   xloc=traveldb.x[ixloc]; yloc=traveldb.y[iyloc]; zloc=traveldb.z[izloc]
                   late, lone = self.coordinate_conversion(xloc, yloc, traveldb.lat0, traveldb.lon0, refell=23)
                   tev=event+datetime.timedelta(seconds=(initwin+itime)*loc.deltat)
                   if num.abs((detections[-1][0]-tev).total_seconds())>120:
                      detections.append([event+datetime.timedelta(seconds=(initwin+itime)*loc.deltat),coheratio[imaxratio],late,lone,zloc])
                      f=open(self.output_path+'/raw_event_catalogue.txt','a')
                      f.write(str(event+datetime.timedelta(seconds=(initwin+itime)*loc.deltat)+datetime.timedelta(seconds=0.))+' '+str(maxloc[imaxratio])+' '+str(coheratio[imaxratio])+' '+ str(late)+' '+str(lone)+' '+str(zloc)+ '\n')
                      print('event :: ',str(event+datetime.timedelta(seconds=(initwin+itime)*loc.deltat)+datetime.timedelta(seconds=0.))+' '+str(coheratio[imaxratio])+' '+ str(late)+' '+str(lone)+' '+str(zloc)+'\n')
                      f.close()
                      #plt.imshow(corrmatrix.reshape(traveldb.nx,traveldb.ny,traveldb.nz)[:,:,5])
                      #plt.show()
            f=open(self.output_path+'/final_event_catalogue.txt','w')
            for event in detections:
                f.write(str(event[0])+' '+str(event[1])+' '+ str(event[2])+' '+str(event[3])+' '+str(event[4])+ '\n')
            f.close()
        print('Ho finito!!!')
