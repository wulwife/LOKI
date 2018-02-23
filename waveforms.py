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
import C_STALTA
import tt_processing
from obspy.core   import read
from obspy.signal import filter
from scipy.signal import hilbert
import datetime


class waveforms:

    def __init__(self, event_path, extension='*', comps=['E','N','Z']):
        if not os.path.isdir(event_path):
            print('Error: data path does not exist')
        try:
            self.load_waveforms(event_path, extension, comps)
            self.station_list()
            self.check_data_consistency()
        except:
            print('Error: data not read for the event: ', event_path)

    def check_data_consistency(self):
        intsamp=1E6
        for comp in streams.keys():
            for sta in st.keys():
                if (int(self.deltat*intsamp)!=int(stream[comp][sta][0]*intsamp)):
                    raise ValueError('Error!! All trace must have the same sampling rate')

    def station_list(self):
        data_stalist=[]
        for comp in self.streams.keys():
            for sta in self.streams[comp].keys():
                if sta not in data_stalist:
                    data_stalist.append(sta)
        self.data_stations=set(data_stalist)

    def load_waveforms(self, event_path, extension, comps):
        files=os.path.join(event_path,extension)
        stream=read(files)
        self.ns=0
        self.evid=str(stream[0].stats.starttime)
        self.deltat=stream[0].stats.delta
        streams={}
        for comp in comps:
            streams[comp]={}
            for tr in stream:
                if tr.stats.channel[-1]==comp:
                    self.ns=max(num.size(tr.data),self.ns)
                    streams[comp][tr.stats.station]=[tr.stats.delta, tr.data]
        self.streams=streams





#add class



    def select_data(self, comp, db_stations, derivative):
        check_data_consistency(self.streams)
        self.stations=(station_list(self.streams) & db_stations)
        self.nstation=len(self.stations)
        tr=num.zeros([self.nstation,self.ns])
        stream=self.streams[comp]
        for i,sta in zip(range(self.nstation),self.stations):
            if sta in stream.keys():
               nstr=num.size(stream[sta])
               tr[i,0:nstr]=stream[sta][1]
               if derivative:
                  tr[i,0:nstr-1]=((tr[i,1:]-tr[i,0:nstr-1])/self.deltat)
                  tr[i,nstr-1]=0.
            else:
               tr[i,:]=1.
        return tr


    def loki_input(self, comp, db_stations, derivative=True):
        if len(comp)==3:
            xtr=self.select_data(comp[0], db_stations, derivative)
            ytr=self.select_data(comp[1], db_stations, derivative)
            ztr=self.select_data(comp[2], db_stations, derivative)
            return (xtr, ytr, ztr)
        elif len(comp)==1:
            ztr=self.select_data(comp, db_stations, derivative)
            return (ztr,)
        else:
            raise ValueError('Traces must have 1 or 3 components!')



    def stacking_function(self):
        traces=self.loki_input()


    def charfunc_erg(self, ztr, ytr, xtr):
        obs_dataV=(ztr**2)
        obs_dataH1=(xtr**2)
        obs_dataH2=(ytr**2)
        return obs_dataV, obs_dataH1, obs_dataH2

    def charfunc_ps(self, xtr, ytr, ztr, epsilon=0.001):
        nsta,nsamp=num.shape(xtr)
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
            obs_dataV[i,:]=obs_dataV[i,:]/num.max(obs_dataV[i,:])
        return obs_dataV, obs_dataH

    def charfunc_s(self, xtr, ytr, epsilon=0.001):
        nsta,nsamp=num.shape(xtr)
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
        return obs_dataH

    def loc_stalta(self, nshort_p, nshort_s, slrat, thres=2):
        tshort_p=nshort_p*self.deltat; tshort_s=nshort_s*self.deltat
        tlong_p=tshort_p*slrat; tlong_s=tshort_s*slrat
        ks_p=self.deltat/tshort_p; kl_p=self.deltat/tlong_p;
        ks_s=self.deltat/tshort_s; kl_s=self.deltat/tlong_s;
        obs_dataP=LOC_STALTA.recursive_stalta(tshort_p, tlong_p, self.deltat, self.obs_dataV, kl_p, ks_p, thres)
        obs_dataS=LOC_STALTA.recursive_stalta(tshort_s, tlong_s, self.deltat, self.obs_dataH, kl_s, ks_s, thres)
        return obs_dataP, obs_dataS

    def det_stalta(self, nshort_p, nshort_s, slrat, thres=2):
        tshort_p=nshort_p*self.deltat; tshort_s=nshort_s*self.deltat
        tlong_p=tshort_p*slrat; tlong_s=tshort_s*slrat
        ks_p=self.deltat/tshort_p; kl_p=self.deltat/tlong_p;
        ks_s=self.deltat/tshort_s; kl_s=self.deltat/tlong_s;
        obs_dataP=C_STALTA.recursive_stalta(tshort_p, tlong_p, self.deltat, self.obs_dataV, kl_p, ks_p, thres)
        obs_dataS=C_STALTA.recursive_stalta(tshort_s, tlong_s, self.deltat, self.obs_dataH, kl_s, ks_s, thres)
        sstalta10,sstalta10=DET_STALTA.recursive_stalta(stas10[i],ltas10[i],tshort_s,tlong_s,self.deltat,obs_dataS1[i,:],0.)
        return obs_dataP, obs_dataS
