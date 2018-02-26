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
