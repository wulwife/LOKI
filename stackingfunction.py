import numpy as num
import lokidata
import DET_STALTA
import LOC_STALTA
from scipy.signal import hilbert


class StackingFunction:

    def __init__(self, lobj, cfunc='ergpca', epsilon=0.001):
       self.nsta=lobj.nstation; self.nsamp=lobj.ns
       if cfunc=='erg':
           self.obs_dataV, self.obs_dataH=self.cfunc_erg(lobj.traces, False)
       elif cfunc=='ergpca':
           self.obs_dataV=self.cfunc_erg(lobj.traces)
           self.obs_dataH=self.cfunc_pca(lobj, epsilon)
       elif cfunc=='pcafull':
           self.obs_dataV, self.obs_dataH=self.cfunc_pcafull(lobj, epsilon)
       else:
           raise ValueError('wrong selection for characteristic function')

    def cfunc_erg(self, traces, ergz=True):
        if ergz:
           self.obs_dataV=(traces[2]**2)
        else:
           self.obs_dataV=(traces[2]**2)
           self.obs_dataH=(traces[0]**2)*(traces[1]**2)

    def cfunc_pcafull(self, lobj, epsilon):
        xtr=lobj.traces[0]; ytr=lobj.traces[1]; ztr=lobj.traces[2]
        nsta=lobj.nstation; nsamp=lobj.ns
        obs_dataH=num.zeros([nsta,nsamp]); obs_dataV=num.zeros([nsta,nsamp])
        obs_dataH1=hilbert(xtr); obs_dataH2=hilbert(ytr); obs_dataH3=hilbert(ztr)
        obs_dataH1C=num.conjugate(obs_dataH1); obs_dataH2C=num.conjugate(obs_dataH2); obs_dataH3C=num.conjugate(obs_dataH3)
        xx=obs_dataH1*obs_dataH1C; xy=obs_dataH1*obs_dataH2C; xz=obs_dataH1*obs_dataH3C
        yx=obs_dataH2*obs_dataH1C; yy=obs_dataH2*obs_dataH2C; yz=obs_dataH2*obs_dataH2C
        zx=obs_dataH3*obs_dataH1C; zy=obs_dataH3*obs_dataH2C; zz=obs_dataH3*obs_dataH3C
        for i in range(nsta):
            for j in range(nsamp):
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

    def cfunc_pca(self, lobj, epsilon=0.001):
        xtr=lobj.traces[0]; ytr=lobj.traces[1]
        nsta=lobj.nstation; nsamp=lobj.ns
        obs_dataH=num.zeros([nsta,nsamp])
        obs_dataH1=hilbert(xtr)
        obs_dataH2=hilbert(ytr)
        obs_dataH1C=num.conjugate(obs_dataH1)
        obs_dataH2C=num.conjugate(obs_dataH2)
        xx=obs_dataH1*obs_dataH1C; xy=obs_dataH1*obs_dataH2C
        yx=obs_dataH2*obs_dataH1C; yy=obs_dataH2*obs_dataH2C
        for i in range(nsta):
            for j in range(nsamp):
                cov=num.array([[xx[i,j], xy[i,j]],[yx[i,j], yy[i,j]]])
                U, s, V = num.linalg.svd(cov, full_matrices=True)
                obs_dataH[i,j]=(s[0]**2)
            obs_dataH[i,:]=(obs_dataH[i,:]/num.max(obs_dataH[i,:]))+epsilon
            plt.plot(obs_dataH[i,:])
            plt.show()
        self.obs_dataH=obs_dataH

    def loc_stalta(self, nshort_p, nshort_s, slrat, thres=2):
        tshort_p=nshort_p*self.deltat; tshort_s=nshort_s*self.deltat
        tlong_p=tshort_p*slrat; tlong_s=tshort_s*slrat
        ks_p=self.deltat/tshort_p; kl_p=self.deltat/tlong_p;
        ks_s=self.deltat/tshort_s; kl_s=self.deltat/tlong_s;
        obs_dataP=LOC_STALTA.recursive_stalta(tshort_p, tlong_p, self.deltat, self.obs_dataV, kl_p, ks_p, thres)
        obs_dataS=LOC_STALTA.recursive_stalta(tshort_s, tlong_s, self.deltat, self.obs_dataH, kl_s, ks_s, thres)
        return obs_dataP, obs_dataS

    def det_stalta(self, nshort_p, nshort_s, slrat, staltap0, staltas0, thres=0.):
        tshort_p=nshort_p*self.deltat; tshort_s=nshort_s*self.deltat
        tlong_p=tshort_p*slrat; tlong_s=tshort_s*slrat
        obs_dataP=num.zeros([self.nsta,self.nsamp])
        obs_dataS=num.zeros([self.nsta,self.nsamp])
        for i in range(self.nsta):
            obs_dataP[i,:],stltp0=DET_STALTA.recstalta(staltap0[i][0], staltap0[i][1], tshort_p, tlong_p, self.deltat, self.obs_dataV[i,:], thres)
            obs_dataS[i,:],stlts0=DET_STALTA.recstalta(staltas0[i][0], staltas0[i][1], tshort_s, tlong_s, self.deltat, self.obs_dataH[i,:], thres)
            staltap0[i][0]=stltp0[0]; staltap0[i][1]=stltp0[1]
            staltas0[i][0]=stlts0[0]; staltas0[i][1]=stlts0[1]
            plt.plot(obs_dataH[i,:])
            plt.show()
        return obs_dataP, obs_dataS, staltap0, staltas0
