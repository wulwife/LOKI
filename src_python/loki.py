import numpy as num
import traveltimes
import waveforms
import stacktraces
import tt_processing
import LatLongUTMconversion
import location
import os
import matplotlib.pyplot as plt

class Loki:
    """docstring for Loki"""

    def __init__(self, data_path, output_path, db_path, hdr_filename, mode='locator'):
        self.data_path=data_path
        self.output_path=output_path
        self.db_path=db_path
        self.hdr_filename=hdr_filename
        if mode=='locator':
            self.data_tree, self.events=self.location_data_struct(self.data_path, self.output_path)
        elif mode=='detector':
            self.data_tree, self.events=self.detection_data_struct(self.data_path, self.output_path)
        else:
            raise ValueError('mode must be "detector" or "locator"')


    def location_data_struct(self, data_path, output_path):
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


    def detection_data_struct(self, data_path, output_path):
        events=[]
        data_tree=[]
        for root,dirs,files in os.walk(data_path):
           if not dirs:
              data_tree.append(root)
              events.append(root.split('/')[-1])
        if not os.path.isdir(output_path):
              os.mkdir(output_path)
        return data_tree, events


    def location(self, extension='*', comp=['E','N','Z'], precision='single', **inputs):
        tshortp_min=inputs['tshortp_min']; tshortp_max=inputs['tshortp_max'];
        tshorts_min=inputs['tshorts_min']; tshorts_max=inputs['tshorts_max'];
        slrat=inputs['slrat']
        npr=inputs['npr']
        ntrial=inputs['ntrial']
        derivative=inputs['derivative']
        model=inputs['model']
        vfunc=inputs['vfunc']
        hfunc=inputs['hfunc']
        epsilon=inputs['epsilon']

        tshortp=num.linspace(tshortp_min,tshortp_max,ntrial)
        tshorts=num.linspace(tshorts_min,tshorts_max,ntrial)


        tobj=traveltimes.Traveltimes(self.db_path, self.hdr_filename)
        tp=tobj.load_traveltimes('P', model, precision)
        ts=tobj.load_traveltimes('S', model, precision)


        for event_path in self.data_tree:

            wobj=waveforms.Waveforms(event_path, extension, comp)
            sobj=stacktraces.Stacktraces(tobj, wobj, **inputs)
            event=sobj.evid
            sobj.cfunc_erg(ergz=False)


            print('Processing to the event folder: ', event_path, event)
            if os.path.isdir(self.output_path+'/'+event):
                continue
            else:
                os.mkdir(self.output_path+'/'+event)


            tp_mod, ts_mod=sobj.time_extractor(tp, ts)
            tp_mod, ts_mod=tt_processing.tt_f2i(sobj.deltat,tp_mod,ts_mod, npr)


            for i in range(ntrial):
                nshort_p=int(tshortp[i]//sobj.deltat)
                nshort_s=int(tshorts[i]//sobj.deltat)
                nlong_s=int(nshort_p*slrat)
                nlong_s=int(nshort_s*slrat)
                obs_dataP, obs_dataS=sobj.loc_stalta(nshort_p, nshort_s, slrat, thres=2)
                corrmatrix=location.stacking(tp_mod, ts_mod, obs_dataP, obs_dataS, npr)
                cmax=num.max(corrmatrix)
                corrmatrix=num.reshape(corrmatrix,(tobj.nx,tobj.ny,tobj.nz))
                (ixloc,iyloc,izloc)=num.unravel_index(num.argmax(corrmatrix),(tobj.nx,tobj.ny,tobj.nz))
                xloc=tobj.x[ixloc]; yloc=tobj.y[iyloc]; zloc=tobj.x[izloc]
                out_file = open(self.output_path+'/'+event+'/'+event+'.loc','a')
                out_file.write(str(i)+' '+str(xloc)+' '+str(yloc)+' '+str(zloc)+' '+str(cmax)+' '+str(nshort_p)+' '+str(nshort_s)+' '+str(slrat)+'\n')
                out_file.close()
                num.save(self.output_path+'/'+event+'/'+'corrmatrix_trial_'+str(i),corrmatrix)
                self.coherence_plot(self.output_path+'/'+event, corrmatrix, tobj.x, tobj.y, tobj.z, i)
            self.catalogue_creation(event, tobj.lat0, tobj.lon0, ntrial)
        print('Location process completed!!!')


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


    def coherence_plot(self, event_path, corrmatrix, xax, yax, zax, itrial, normalization=False):
        nx,ny,nz=num.shape(corrmatrix)
        CXY=num.zeros([ny, nx])
        for i in range(ny):
            for j in range(nx):
                CXY[i,j]=num.max(corrmatrix[j,i,:])

        CXZ=num.zeros([nz, nx])
        for i in range(nz):
            for j in range(nx):
			    CXZ[i,j]=num.max(corrmatrix[j,:,i])

        CYZ=num.zeros([nz, ny])
        for i in range(nz):
            for j in range(ny):
                CYZ[i,j]=num.max(corrmatrix[:,j,i])

        if normalization:
           nrm=Normalize(vmin=0., vmax=1.)
        else:
           nrm=None


        fig = plt.figure()
        fig.suptitle('Coherence matrix X-Y', fontsize=14, fontweight='bold')
        ax = fig.gca()
        cmap = plt.cm.get_cmap('jet', 100)
        cs = plt.contourf(xax, yax, CXY, 20, cmap=cmap, interpolation='bilinear', norm=nrm)
        ax.set_xlabel('X (km)')
        ax.set_ylabel('Y (km)')
        cbar = plt.colorbar(cs)
        plt.axes().set_aspect('equal')
        plt.savefig(event_path+'/'+'Coherence_matrix_xy'+str(itrial)+'.eps')


        fig = plt.figure()
        fig.suptitle('Coherence matrix X-Z', fontsize=14, fontweight='bold')
        ax = fig.gca()
        cmap = plt.cm.get_cmap('jet', 100)
        cs = plt.contourf(xax, zax, CXZ, 20, cmap=cmap, interpolation='bilinear', norm=nrm)
        ax.set_xlabel('X (km)')
        ax.set_ylabel('Z (km)')
        cbar = plt.colorbar(cs)
        ax.invert_yaxis()
        plt.axes().set_aspect('equal')
        plt.savefig(event_path+'/'+'Coherence_matrix_xz'+str(itrial)+'.eps')


        fig = plt.figure()
        fig.suptitle('Coherence matrix Y-Z', fontsize=14, fontweight='bold')
        ax = fig.gca()
        cmap = plt.cm.get_cmap('jet', 100)
        cs = plt.contourf(yax, zax, CYZ, 20, cmap=cmap, interpolation='bilinear', norm=nrm)
        ax.set_xlabel('Y(km)')
        ax.set_ylabel('Z (km)')
        ax.invert_yaxis()
        cbar = plt.colorbar(cs)
        plt.axes().set_aspect('equal')
        plt.savefig(event_path+'/'+'Coherence_matrix_yz'+str(itrial)+'.eps')
        plt.close("all")
