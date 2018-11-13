import .src_python.loki
import numpy as num

db_path='/Users/francesco/Data/KOREA/korea_time/korea_herr'
data_path='/Users/francesco/Data/KOREA/2016-09-12-mww54-south-korea'
output_path='/Users/francesco/Desktop/test'
hdr_filename='interpolated.header.hdr'
inputs={}
inputs['tshortp_min']=0.1
inputs['tshortp_max']=0.2
inputs['tshorts_min']=0.2
inputs['tshorts_max']=0.4
inputs['slrat']=2
inputs['npr']=2
inputs['ntrial']=1
inputs['derivative']=True
inputs['model']='interpolated'
precision='single'
comp=['E','N','Z']
extension='*.SAC'
#test

l1=loki.Loki(data_path, output_path, db_path, hdr_filename, mode='locator')
l1.location(extension, comp, precision, **inputs)
