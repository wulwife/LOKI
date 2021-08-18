import os
import numpy as num
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import datetime
import copy
#
from loki import traveltimes
from loki import waveforms
from loki import stacktraces
from loki import LatLongUTMconversion
import tt_processing                    # C
import location_t0                         # C
#import location_t0_plus                         # C


class Loki:
    """docstring for Loki"""

    def __init__(self, data_path, output_path, db_path, hdr_filename, mode='locator'):
        self.data_path = data_path
        self.output_path = output_path
        self.db_path = db_path
        self.hdr_filename = hdr_filename
        if mode == 'locator':
            self.data_tree, self.events = self.location_data_struct(self.data_path, self.output_path)
        elif mode == 'detector':
            self.data_tree, self.events = self.detection_data_struct(self.data_path, self.output_path)
        else:
            raise ValueError('mode must be "detector" or "locator"')

    def location_data_struct(self, data_path, output_path):
        events=[]
        data_tree=[]
        for root, dirs, files in os.walk(data_path):
            if not dirs:
                data_tree.append(root)
        data_tree.sort()
        events = [idtree.split('/')[-1] for idtree in data_tree]
        if not os.path.isdir(output_path):
            os.mkdir(output_path)
        return data_tree, events

    def detection_data_struct(self, data_path, output_path):
        events = []
        data_tree = []
        for root, dirs, files in os.walk(data_path):
            if not dirs:
                data_tree.append(root)
        data_tree.sort()
        events = [idtree.split('/')[-1] for idtree in data_tree]
        if not os.path.isdir(output_path):
            os.mkdir(output_path)
        return data_tree, events

    def location(self, extension='*', comp=['E', 'N', 'Z'], precision='single', **inputs):
        if 'tshortp_min' in inputs:
            # need to calculate STA/LTA for stacking
            STALTA = True
            tshortp_min = inputs['tshortp_min']
            tshortp_max = inputs['tshortp_max']
            tshorts_min = inputs['tshorts_min']
            tshorts_max = inputs['tshorts_max']
            slrat = inputs['slrat']
            ntrial = inputs['ntrial']
        
            tshortp = num.linspace(tshortp_min, tshortp_max, ntrial)
            tshorts = num.linspace(tshorts_min, tshorts_max, ntrial)
        else:
            # no need to calculate STA/LTA ratio for stacking
            STALTA = False
            ntrial = 1
            
        npr = inputs['npr']
        model = inputs['model']
        if 'freq' in inputs:
            freq = inputs['freq']
        else:
            freq = None  # no filtering
        
        # load traveltime data set
        tobj = traveltimes.Traveltimes(self.db_path, self.hdr_filename)
        tp = tobj.load_traveltimes('P', model, precision)
        ts = tobj.load_traveltimes('S', model, precision)

        for event_path in self.data_tree:
            wobj = waveforms.Waveforms(event_path, extension, comp, freq)
            sobj = stacktraces.Stacktraces(tobj, wobj, **inputs)
            event = sobj.evid

            print('Processing to the event folder: ', event_path, event)
            if os.path.isdir(self.output_path+'/'+event):
                continue
            else:
                os.mkdir(self.output_path+'/'+event)

            tp_modse, ts_modse = sobj.time_extractor(tp, ts)  # traveltime table in second
            tp_mod, ts_mod = tt_processing.tt_f2i(sobj.deltat, tp_modse, ts_modse, npr)  # traveltime table in time point, for each imaging point traveltimes have substracted the minimal P traveltime

            cmax_pre = -1.0
            for i in range(ntrial):
                if STALTA:
                    # need to calculate STA/LTA from the characteristic funtion
                    # then stack the STA/LTA for imaging
                    nshort_p = int(tshortp[i]//sobj.deltat)
                    nshort_s = int(tshorts[i]//sobj.deltat)
                    obs_dataP, obs_dataS = sobj.loc_stalta(nshort_p, nshort_s, slrat, norm=1)
                else:
                    # no need to calculate STA/LTA 
                    # directly stack the characteristic function for imaging
                    obs_dataP = sobj.obs_dataV  # vertical -> P
                    obs_dataS = sobj.obs_dataH  # horizontal -> S

                iloctime, corrmatrix = location_t0.stacking(tp_mod, ts_mod, obs_dataP, obs_dataS, npr)
                evtpmin = num.amin(tp_modse[iloctime[0],:])
                event_t0 = sobj.dtime_max + datetime.timedelta(seconds=iloctime[1]*sobj.deltat) - datetime.timedelta(seconds=evtpmin)  # event origin time
                event_t0s = (event_t0).isoformat()
                # corrmatrix is the stacking matrix, in 1D format but can be 
                # reformat to 3D format, each point saves the maximum stacking 
                # value during this calculation time period
                cmax = num.max(corrmatrix)
                corrmatrix = num.reshape(corrmatrix,(tobj.nx,tobj.ny,tobj.nz))
                (ixloc, iyloc, izloc) = num.unravel_index(num.argmax(corrmatrix),(tobj.nx,tobj.ny,tobj.nz))
                xloc = tobj.x[ixloc]
                yloc = tobj.y[iyloc]
                zloc = tobj.z[izloc]
                if ntrial > 1:
                    out_file = open(self.output_path+'/'+event+'/'+event+'.loc', 'a')
                else:
                    out_file = open(self.output_path+'/'+event+'/'+event_t0s+'.loc', 'a')
                if STALTA:
                    out_file.write(str(i)+' '+str(xloc)+' '+str(yloc)+' '+str(zloc)+' '+str(cmax)+' '+str(nshort_p)+' '+str(nshort_s)+' '+str(slrat)+'\n')
                else:
                    out_file.write(str(i)+' '+str(xloc)+' '+str(yloc)+' '+str(zloc)+' '+str(cmax)+'\n')
                out_file.close()
                num.save(self.output_path+'/'+event+'/'+'corrmatrix_trial_'+str(i),corrmatrix)
                self.coherence_plot(self.output_path+'/'+event, corrmatrix, tobj.x, tobj.y, tobj.z, i)
                if cmax > cmax_pre:
                    event_t0s_final = copy.deepcopy(event_t0s)
                    cmax_pre = copy.deepcopy(cmax)
            
            self.catalogue_creation(event, event_t0s_final, tobj.lat0, tobj.lon0, ntrial, corrmatrix)
        print('Location process completed!!!')

    def catalogue_creation(self, event, event_t0s, lat0, lon0, ntrial, corrmatrix, refell=23):
        (zorig, eorig, norig) = LatLongUTMconversion.LLtoUTM(refell, lat0, lon0) #da adeguare in python 3
        if (ntrial > 1):
            ev_file = self.output_path+'/'+event+'/'+event+'.loc'
            data = num.loadtxt(ev_file)
            w = num.sum(data[:, 4])
            xb = ((num.dot(data[:, 1], data[:, 4])/w)*1000)+eorig
            yb = ((num.dot(data[:, 2], data[:, 4])/w)*1000)+norig
            late, lone = LatLongUTMconversion.UTMtoLL(refell, yb, xb, zorig)
            zb = num.dot(data[:, 3], data[:, 4])/w  # depth in km
            cb = num.mean(data[:, 4])  # the mean coherence over the ntrial realizations
            cmax = num.max(data[:, 4])  # the maximum coherence over the ntrial realizations
            merr = num.vstack((data[:, 1], data[:, 2], data[:, 3]))
            err = num.cov(merr)
            errmax = num.sqrt(num.max(num.linalg.eigvals(err)))
        else:
            ev_file = self.output_path+'/'+event+'/'+event_t0s+'.loc'
            data = num.loadtxt(ev_file)
            late, lone = LatLongUTMconversion.UTMtoLL(refell, (data[2]*1000)+norig, (data[1]*1000)+eorig, zorig)  # latitude, longitude
            zb = data[3]  # depth in km
            cmax = data[4]  # the maximum coherence over the 3D corrmatrix
            
            # nomalize corrmatrix first, let minimal->1, maximum->2
            n1 = 1.0  # minimal limit
            n2 = 2.0  # maximum limit
            dmax = num.amax(corrmatrix, axis=None, keepdims=True)
            dmin = num.amin(corrmatrix, axis=None, keepdims=True)
            k = (n2-n1)/(dmax-dmin)
            b = (dmax*n1-dmin*n2)/(dmax-dmin)
            corrmatrix = k*corrmatrix + b
            
            errmax = num.std(corrmatrix, axis=None)  # the coherence standard deviation over the 3D corrmatrix
            cb = num.median(corrmatrix, axis=None)  # the median coherence over the 3D corrmatrix
            
        f = open(self.output_path+'/'+'catalogue', 'a')
        f.write(event_t0s+'    '+str(late)+'   '+str(lone)+'   '+str(zb)+'   '+str(errmax)+'   '+str(cb)+'   '+str(cmax)+'\n')
        f.close()

    def coherence_plot(self, event_path, corrmatrix, xax, yax, zax, itrial, normalization=False):
        nx, ny, nz = num.shape(corrmatrix)
        CXY = num.zeros([ny, nx])
        for i in range(ny):
            for j in range(nx):
                CXY[i,j]=num.max(corrmatrix[j,i,:])

        CXZ = num.zeros([nz, nx])
        for i in range(nz):
            for j in range(nx):
                CXZ[i, j] = num.max(corrmatrix[j,:,i])

        CYZ = num.zeros([nz, ny])
        for i in range(nz):
            for j in range(ny):
                CYZ[i, j] = num.max(corrmatrix[:, j, i])

        if normalization:
            nrm = Normalize(vmin=0., vmax=1.)
        else:
            nrm = None

        fig = plt.figure()
        fig.suptitle('Coherence matrix X-Y', fontsize=14, fontweight='bold')
        ax = fig.gca()
        cmap = plt.cm.get_cmap('jet', 100)
        cs = plt.contourf(xax, yax, CXY, 20, cmap=cmap, interpolation='bilinear', norm=nrm)
        ax.set_xlabel('X (km)')
        ax.set_ylabel('Y (km)')
        cbar = plt.colorbar(cs)
        ax.set_aspect('equal')
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
        ax.set_aspect('equal')
        plt.savefig(event_path+'/'+'Coherence_matrix_xz'+str(itrial)+'.eps')

        fig = plt.figure()
        fig.suptitle('Coherence matrix Y-Z', fontsize=14, fontweight='bold')
        ax = fig.gca()
        cmap = plt.cm.get_cmap('jet', 100)
        cs = plt.contourf(yax, zax, CYZ, 20, cmap=cmap, interpolation='bilinear', norm=nrm)
        ax.set_xlabel('Y (km)')
        ax.set_ylabel('Z (km)')
        cbar = plt.colorbar(cs)
        ax.invert_yaxis()
        ax.set_aspect('equal')
        plt.savefig(event_path+'/'+'Coherence_matrix_yz'+str(itrial)+'.eps')
        plt.close("all")
