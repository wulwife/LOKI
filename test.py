import waveforms
import traveltimes
import lokidata
import stackingfunction
import numpy as num

tt='/Users/francesco/Desktop/KOREA/korea_time/korea_herr'
ev='/Users/francesco/Desktop/KOREA/2016-09-12-mww54-south-korea'

t1=traveltimes.Traveltimes(tt,'header.hdr')
w1=waveforms.Waveforms(ev,extension='*.SAC', comps=['E','N','Z'])
l1=lokidata.LokiData(t1,w1)
stackingfunction.StackingFunction(l1, cfunc='pcafull', epsilon=0.001)
nshort_p=0.1; nshort_s=0.2; slrat=2.0
staltap0=[1.0,1.0]; staltas0=[1.0,1.0]
stackingfunction.det_stalta(nshort_p, nshort_s, slrat, staltap0, staltas0, 0.)
