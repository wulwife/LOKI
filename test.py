import waveforms
import traveltimes
import lokidata

tt='/Users/francesco/Desktop/KOREA/korea_time/korea_herr'
ev='/Users/francesco/Desktop/KOREA/2016-09-12-mww54-south-korea'

t1=traveltimes.Traveltimes(tt,'header.hdr')
w1=waveforms.Waveforms(ev,extension='*.SAC', comps=['E','N','Z'])
l1=lokidata.LokiData(t1,w1)
print('cucu')
