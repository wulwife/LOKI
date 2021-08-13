import os
from obspy.core import read
from datetime import datetime


class Waveforms:

    def __init__(self, event_path, extension='*', comps=['E','N','Z'], freq=None):
        if not os.path.isdir(event_path):
            raise ValueError('Error: data path does not exist')
        try:
            self.load_waveforms(event_path, extension, comps, freq)
        except:
            raise WaveformLoadingError('Error: data not read for the event: %s' %(event_path))
        self.station_list()

    def station_list(self):
        data_stalist=[]
        for comp in (self.stream).keys():
            for sta in (self.stream[comp]).keys():
                if sta not in data_stalist:
                    data_stalist.append(sta)
        self.data_stations=set(data_stalist)

    def load_waveforms(self, event_path, extension, comps, freq):
        files=os.path.join(event_path,extension)
        traces=read(files)
        
        if freq:
            traces.detrend('demean')
            traces.detrend('linear')
            if len(freq) == 1:
                traces.filter("highpass", freq=freq[0])
            elif len(freq) == 2:
                traces.filter("bandpass", freqmin=freq[0], freqmax=freq[1])
        
        self.stream={}
        for comp in comps:
            self.stream[comp]={}
            for tr in traces:
                if tr.stats.channel[-1]==comp:
                    dtime=datetime.strptime(str(tr.stats.starttime),"%Y-%m-%dT%H:%M:%S.%fZ")
                    self.stream[comp][tr.stats.station]=[dtime, tr.stats.delta, tr.data]


class WaveformLoadingError(Exception):
    pass
