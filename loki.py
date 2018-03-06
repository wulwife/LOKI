import traveltimes
import waveforms
from datetime import datetime
import numpy as num
import matplotlib.pyplot as plt
import DET_STALTA
import LOC_STALTA

class Loki:

    def __init__(self, tobj, wobj, derivative=True):
        #check input objects
        self.check_sampling_rate(wobj)
        self.check_starting_time(wobj)
        self.loki_input(wobj, tobj, derivative)

    def check_sampling_rate(self,wobj):
        intsamp=1E6
        deltas=[]
        for comp in (wobj.stream).keys():
            for sta in (wobj.stream[comp]).keys():
                deltas.append(wobj.stream[comp][sta][1])
        deltas=num.array(deltas)
        ideltas=num.unique((deltas*intsamp).astype(int))
        if num.size(ideltas)==1:
           self.deltat=deltas[0]
        else:
           raise ValueError('Error!! All trace must have the same sampling rate')

    def check_starting_time(self,wobj):
        intsamp=1E6
        dtimes=[]
        self.ns=0
        for comp in (wobj.stream).keys():
            for sta in (wobj.stream[comp]).keys():
                dtimes.append(wobj.stream[comp][sta][0])
                if self.ns<num.size(wobj.stream[comp][sta][2]):
                   self.ns=num.size(wobj.stream[comp][sta][2])
        self.dtime_max=max(dtimes)
        self.evid=(self.dtime_max).isoformat()

    def loki_input(self, wobj, tobj, derivative):
        self.comp=tuple((wobj.stream).keys())
        if len(self.comp)==3:
            self.xtr=self.select_data(self.comp[0], wobj, tobj.db_stations, derivative)
            self.ytr=self.select_data(self.comp[1], wobj, tobj.db_stations, derivative)
            self.ztr=self.select_data(self.comp[2], wobj, tobj.db_stations, derivative)
        elif len(self.comp)==1:
            self.ztr=self.select_data(comp, wobj.data_stations, db_stations, derivative)
        else:
            raise ValueError('Traces must have 1 or 3 components!')

    def select_data(self, comp, wobj, db_stations, derivative):
        self.stations=tuple(wobj.data_stations & db_stations)
        self.nstation=num.size(self.stations)
        tr=num.zeros([self.nstation,self.ns])
        stream=wobj.stream[comp]
        for i,sta in enumerate(self.stations):
            nstr=num.size(stream[sta][2])
            idt=num.int((self.dtime_max-stream[sta][0]).total_seconds()/self.deltat)
            tr[i,0:nstr-idt]=stream[sta][2][idt:]
            if derivative:
               tr[i,1:self.ns]=((tr[i,1:]-tr[i,0:self.ns-1])/self.deltat)
               tr[i,0]=0.
            else:
               tr[i,:]=1.
        return tr

    def time_extractor(self, tp, ts):
        nxyz= num.size(tp[self.stations[0]])
        tp_mod=num.zeros([nxyz,self.nstation]) #optimize with pointers to tp[sta] and ts[sta]
        ts_mod=num.zeros([nxyz,self.nstation])
        for i,sta in enumerate(self.stations):
            tp_mod[:,i]=tp[sta]
            ts_mod[:,i]=ts[sta]
        return (tp_mod, ts_mod)

    def cfunc_erg(self, ergz=True):
        if ergz:
           self.obs_dataV=(self.ztr**2)
        else:
           self.obs_dataV=(self.ztr**2)
           self.obs_dataH=(self.xtr**2)*(self.ytr**2)

    def cfunc_pcafull(self, epsilon=0.001):
        obs_dataH=num.zeros([self.nstation,self.ns]); obs_dataV=num.zeros([self.nstation,self.ns])
        obs_dataH1=hilbert(self.xtr); obs_dataH2=hilbert(self.ytr); obs_dataH3=hilbert(self.ztr)
        obs_dataH1C=num.conjugate(obs_dataH1); obs_dataH2C=num.conjugate(obs_dataH2); obs_dataH3C=num.conjugate(obs_dataH3)
        xx=obs_dataH1*obs_dataH1C; xy=obs_dataH1*obs_dataH2C; xz=obs_dataH1*obs_dataH3C
        yx=obs_dataH2*obs_dataH1C; yy=obs_dataH2*obs_dataH2C; yz=obs_dataH2*obs_dataH2C
        zx=obs_dataH3*obs_dataH1C; zy=obs_dataH3*obs_dataH2C; zz=obs_dataH3*obs_dataH3C
        for i in range(self.nstation):
            for j in range(self.ns):
                cov3d=num.array([[xx[i,j], xy[i,j], xz[i,j]],[yx[i,j], yy[i,j], yz[i,j]],[zx[i,j],zy[i,j],zz[i,j]]])
                cov2d=num.array([[xx[i,j], xy[i,j]],[yx[i,j], yy[i,j]]])
                U2d, s2d, V2d = num.linalg.svd(cov2d, full_matrices=True)
                U3d, s3d, V3d = num.linalg.svd(cov3d, full_matrices=True)
                obs_dataV[i,j]=(s3d[0]**2)
                obs_dataH[i,j]=(s2d[0]**2)
            obs_dataH[i,:]=(obs_dataH[i,:]/num.max(obs_dataH[i,:]))+epsilon
            obs_dataV[i,:]=(obs_dataV[i,:]/num.max(obs_dataV[i,:]))
        self.obs_dataH=obs_dataH
        self.obs_dataV=obs_dataV

    def cfunc_pca(self, epsilon=0.001):
        obs_dataH=num.zeros([self.nstation,self.ns])
        obs_dataH1=hilbert(self.xtr)
        obs_dataH2=hilbert(self.ytr)
        obs_dataH1C=num.conjugate(obs_dataH1)
        obs_dataH2C=num.conjugate(obs_dataH2)
        xx=obs_dataH1*obs_dataH1C; xy=obs_dataH1*obs_dataH2C
        yx=obs_dataH2*obs_dataH1C; yy=obs_dataH2*obs_dataH2C
        for i in range(self.nstation):
            for j in range(self.ns):
                cov=num.array([[xx[i,j], xy[i,j]],[yx[i,j], yy[i,j]]])
                U, s, V = num.linalg.svd(cov, full_matrices=True)
                obs_dataH[i,j]=(s[0]**2)
            obs_dataH[i,:]=(obs_dataH[i,:]/num.max(obs_dataH[i,:]))+epsilon
        self.obs_dataH=obs_dataH


#locations class?
    def loc_stalta(self, nshort_p, nshort_s, slrat, thres=2):
        tshort_p=nshort_p*self.deltat; tshort_s=nshort_s*self.deltat
        tlong_p=tshort_p*slrat; tlong_s=tshort_s*slrat
        ks_p=self.deltat/tshort_p; kl_p=self.deltat/tlong_p;
        ks_s=self.deltat/tshort_s; kl_s=self.deltat/tlong_s;
        obs_dataP=LOC_STALTA.recursive_stalta(tshort_p, tlong_p, self.deltat, self.obs_dataV, kl_p, ks_p, thres)
        obs_dataS=LOC_STALTA.recursive_stalta(tshort_s, tlong_s, self.deltat, self.obs_dataH, kl_s, ks_s, thres)
        return obs_dataP, obs_dataS


#detection class?
    def det_stalta(self, nshort_p, nshort_s, slrat, staltap0, staltas0, thres=0.):
        tshort_p=nshort_p*self.deltat; tshort_s=nshort_s*self.deltat
        tlong_p=tshort_p*slrat; tlong_s=tshort_s*slrat
        obs_dataP=num.zeros([self.nstation,self.ns])
        obs_dataS=num.zeros([self.nstation,self.ns])
        for i in range(self.nstation):
            obs_dataP[i,:],stltp0=DET_STALTA.recursive_stalta(staltap0[i][0], staltap0[i][1], tshort_p, tlong_p, self.deltat, self.obs_dataV[i,:], thres)
            obs_dataS[i,:],stlts0=DET_STALTA.recursive_stalta(staltas0[i][0], staltas0[i][1], tshort_s, tlong_s, self.deltat, self.obs_dataH[i,:], thres)
            staltap0[i][0]=stltp0[0]; staltap0[i][1]=stltp0[1]
            staltas0[i][0]=stlts0[0]; staltas0[i][1]=stlts0[1]
        return obs_dataP, obs_dataS, staltap0, staltas0
