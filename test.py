import waveforms
import traveltimes
import loki
import numpy as num

tt='/Users/francesco/Desktop/KOREA/korea_time/korea_herr'
ev='/Users/francesco/Desktop/KOREA/2016-09-12-mww54-south-korea'

t1=traveltimes.Traveltimes(tt,'header.hdr')
w1=waveforms.Waveforms(ev,extension='*.SAC', comps=['E','N','Z'])
l1=loki.Loki(t1,w1)
l1.cfunc_erg(False)
nshort_p=0.1; nshort_s=0.2; slrat=2.0
staltap0=[1.0,1.0]; staltas0=[1.0,1.0]
l1.det_stalta(nshort_p, nshort_s, slrat, staltap0, staltas0, 0.)
