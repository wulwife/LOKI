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
from obspy.core   import read
import datetime

class WaveformLoadingError(Exception):
    pass

class Waveforms:

    def __init__(self, event_path, extension='*', comps=['E','N','Z']):
        if not os.path.isdir(event_path):
            raise ValueError('Error: data path does not exist')
        try:
            self.load_waveforms(event_path, extension, comps)
        except:
            raise WaveformLoadingError('Error: data not read for the event: %s' %(event_path))
        self.station_list()
        self.check_data_consistency()

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
                    self.ns=num.maximum(num.size(tr.data),self.ns)
                    streams[comp][tr.stats.station]=[tr.stats.delta, tr.data]
        self.streams=streams
