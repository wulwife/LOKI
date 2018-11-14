import src_python.loki
import numpy as num

db_path='./Test_dataset/Traveltimes'
data_path='./Test_dataset/Data'
output_path='./Test_dataset/output'
hdr_filename='header.hdr'
inputs={}
inputs['tshortp_min']=0.1
inputs['tshortp_max']=0.2
inputs['tshorts_min']=0.2
inputs['tshorts_max']=0.4
inputs['slrat']=2
inputs['npr']=2
inputs['ntrial']=1
inputs['derivative']=True
inputs['vfunc']='erg'
inputs['hfunc']='pca'
inputs['model']='interpolated'
precision='single'
comp=['E','N','Z']
extension='*'
#test

l1=loki.Loki(data_path, output_path, db_path, hdr_filename, mode='locator')
l1.location(extension, comp, precision, **inputs)
