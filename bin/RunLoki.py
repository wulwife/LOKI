#!/usr/bin/env python

"""
To call LOKI you need this file and make it executable (or call python)

"""

from loki.loki import Loki

db_path = 'TRAVELTIMES-PATH'
data_path = 'WAVEFORMS-PATH'
output_path = 'OUTPUT-DIR-PATH'
hdr_filename = 'HEADER-FILENAME'
#
inputs = {}
inputs['tshortp_min'] = 0.1
inputs['tshortp_max'] = 0.1
inputs['tshorts_min'] = 0.15
inputs['tshorts_max'] = 0.15
inputs['slrat'] = 2
inputs['npr'] = 2
inputs['ntrial'] = 1
inputs['derivative'] = True
inputs['vfunc'] = 'erg'
inputs['hfunc'] = 'pca'
inputs['model'] = 'layer'
inputs['epsilon'] = 0.001
#
precision = 'single'
comp = ['E', 'N', 'Z']
extension = '*'

# =========  Call

l1 = Loki(data_path, output_path, db_path, hdr_filename, mode='locator')
l1.location(extension, comp, precision, **inputs)
