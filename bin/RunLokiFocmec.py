#!/usr/bin/env python

"""
LOKI location + focal mechanism by polarity/radiation-pattern corrected
waveform stacking. Same carbon-print as RunLoki.py, plus the 'focmec' keys.

Requirements for the focal mechanism step:
  - 3-component data (E, N, Z), traveltime tables for P and S
  - zero-phase processing is done internally (fm_freq); give velocity
    ('vel', default) or displacement ('disp') data
"""

from loki.loki import Loki

db_path = 'TRAVELTIMES-PATH'
data_path = 'WAVEFORMS-PATH'
output_path = 'OUTPUT-DIR-PATH'
hdr_filename = 'HEADER-FILENAME'
#
inputs = {}
inputs['tshortp_min'] = 0.2
inputs['tshortp_max'] = 0.2
inputs['tshorts_min'] = 0.3
inputs['tshorts_max'] = 0.3
inputs['slrat'] = 2
inputs['npr'] = 2
inputs['ntrial'] = 1
inputs['derivative'] = True
inputs['vfunc'] = 'erg'
inputs['hfunc'] = 'pca'
inputs['model'] = 'layer'
inputs['epsilon'] = 0.001
inputs['freq'] = [1.0, 8.0]          # band used for the location
#
# ----- focal mechanism (all keys optional except 'focmec')
inputs['focmec'] = True
inputs['fm_method'] = 'amp'          # 'amp' (radiation-amplitude weighted), 'pol' (P polarity only)
                                     # or 'area' (simplified waveforms: signed areas in the windows)
inputs['fm_area'] = 'ratio'          # 'area' only: 'pol' (sign of P area), 'amp' (signed areas) or
                                     # 'ratio' (per-station polarities + S/P, SV/SH amplitude ratios)
inputs['fm_awin'] = [-0.15, 0.15]    # 'area' only: integration window around the arrival (s)
inputs['fm_freesurf'] = True         # free-surface P-SV conversion at the receiver (uses fm_vrec)
inputs['fm_phases'] = 'PS'           # 'PS' (P+SV+SH) or 'P'   ('amp' only)
inputs['fm_freq'] = [0.5, 3.0]       # zero-phase band on displacement (below the corner frequency)
inputs['fm_dtype'] = 'vel'           # input data type: 'vel' or 'disp'
inputs['fm_twinp'] = [0.2, 0.6]      # P window, s before/after the theoretical arrival
inputs['fm_twins'] = [0.2, 0.8]      # S window
inputs['fm_tshift'] = 0.2            # common time shift searched (s)
inputs['fm_nneigh'] = 1              # joint search over +-1 grid nodes around the LOKI location
inputs['fm_dsdr'] = [5., 5., 5.]     # strike, dip, rake steps (deg)
inputs['fm_vrec'] = [4.5, 2.5]       # Vp, Vs below the stations (Snell + shear-wave window); None = straight rays
inputs['fm_nboot'] = 100             # station bootstrap for uncertainty
#
precision = 'single'
comp = ['E', 'N', 'Z']
extension = '*'

# =========  Call

l1 = Loki(data_path, output_path, db_path, hdr_filename, mode='locator')
l1.location(extension, comp, precision, **inputs)
# outputs: <output>/<event>/<event>.fm, <event>_fm.png and <output>/catalogue_focmec
