import os
import numpy as num
import pickle as pkl
from loki.traveltimes import Traveltimes
from loki.waveforms import Waveforms
from loki.stacktraces import Stacktraces
import loki.latlon2cart as ll2c
import matplotlib.pyplot as plt

# define inputs

model_path = '/home/francesco/Seistools/LOKI/tests/Test_dataset/Traveltimes'
data_path = '/home/francesco/Seistools/LOKI/tests/Test_dataset/Data/e0060.188.07'
hdr_filename = 'header.hdr'

tshortp=0.1
tshorts=0.1
slrat=2
model='layer'

inputs = {}
inputs['derivative'] = True
inputs['vfunc'] = 'erg'
inputs['hfunc'] = 'pca'
inputs['epsilon'] = 0.001

evlat,evlon,depth=51.65482947253372,7.731559863077049,1.35


# extract traveltimes

times=Traveltimes(model_path,hdr_filename)
ix,iy,iz=times.event_indexes(evlat,evlon,depth)
tp=times.load_traveltimes('P','layer')
ts=times.load_traveltimes('S','layer')

# extract waveforms

traces=Waveforms(data_path, extension='*', comps=['E','N','Z'], freq=None)
straces=Stacktraces(times, traces, **inputs)
stations=straces.stations

# calculate stacking functions

nshort_p = int(tshortp//straces.deltat)
nshort_s = int(tshorts//straces.deltat)
staltaP, staltaS=straces.loc_stalta(nshort_p, nshort_s, slrat, norm=1)


# estimate correcttions

tpcalc={}; tscalc={}
tpobs={}; tsobs={}
tpminobs=1E10
tpmincalc=1E10

for i,sta in enumerate(stations):
    tpcalc[sta]=tp[sta].reshape(times.nx,times.ny,times.nz)[ix,iy,iz]
    if tpcalc[sta]<tpmincalc:
         tpmincalc=tpcalc[sta]
         stamincalc=sta
    tpobs[sta]=num.argmax(staltaP[i])*straces.deltat
    if tpobs[sta]<tpminobs:
         tpminobs=tpobs[sta]
         staminobs=sta
    tsobs[sta]=num.argmax(staltaS[i])*straces.deltat
    tscalc[sta]=ts[sta].reshape(times.nx,times.ny,times.nz)[ix,iy,iz]

corrections_p={}
corrections_s={}
for sta in stations:
    corrections_p[sta]=(tpobs[sta]-tpminobs)-(tpcalc[sta]-tpmincalc)
    corrections_s[sta]=(tsobs[sta]-tpminobs)-(tscalc[sta]-tpmincalc)


#plot corrections

for i,sta in enumerate(stations):
    tpi=num.round((tpcalc[sta]+tpminobs-tpmincalc)/straces.deltat)
    tsi=num.round((tscalc[sta]+tpminobs-tpmincalc)/straces.deltat)
    plt.plot(staltaP[i],'b')
    plt.plot(tpi,1,'y*')
    plt.plot(staltaS[i],'r')
    plt.plot(tsi,1,'g*')
    plt.show()

#save traveltimes with corrections

tp_master=times.apply_master_event_correction('P', corrections_p, label='layer', precision='single')
times.save_ttdb(tp_master,'P','master')

ts_master=times.apply_master_event_correction('S', corrections_s, label='layer', precision='single')
times.save_ttdb(ts_master,'S','master')
