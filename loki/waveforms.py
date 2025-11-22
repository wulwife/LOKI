import os
from obspy.core import read
from datetime import datetime


class Waveforms:

    def __init__(self, data_path, extension='*', comps=['E','N','Z'], freq=None, sds=False, tini=None, window=None, overlap=0.0, network='*'):
        if not os.path.isdir(data_path):
            raise ValueError('Error: data path does not exist')
        if not sds:
            try:
                self.load_waveforms(data_path, extension, comps, freq)
            except:
                raise WaveformLoadingError('Error: data not read for the event: %s' %(data_path))
        else:
            if tini is None or window is None:
                raise ValueError("For sds=True you must provide 'tini' and 'window'")
            try:
                self.load_waveforms_from_sds(data_path, comps, freq, tini, window, overlap, network)
            except:
                raise WaveformLoadingError('Error: data from sds not read in the path: %s' %(data_path))
        self.station_list()


    def station_list(self):
        data_stalist=[]
        for comp in (self.stream).keys():
            for sta in (self.stream[comp]).keys():
                if sta not in data_stalist:
                    data_stalist.append(sta)
        self.data_stations=set(data_stalist)


    def load_waveforms(self, data_path, extension, comps, freq):
        files=os.path.join(data_path,extension)
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


    def load_waveforms_from_sds(self, data_path, comps, freq, tini, window, overlap, network):
        """
        Legge dati continui da struttura SDS, unisce segmenti, riempie i gap con zeri
        e restituisce finestre pulite per ciascuna componente/stazione.

        - Gaps -> riempiti a 0
        - Segmenti sovrapposti -> fusi (merge standard ObsPy)
        - Trim finale esatto a [tini, tini+window] con pad=True (zeri)
        - Facoltativi: resampling e filtro
        """
        from obspy.clients.filesystem.sds import Client
        from obspy import UTCDateTime

        # finestra di lettura estesa per includere overlap a sinistra
        t0 = tini - overlap
        t1 = tini + window

        client = Client(data_path, sds_type='D', format='MSEED', fileborder_seconds=30, fileborder_samples=5000)

        # Leggi tutte le HH/EH/BH componenti Z/N/E (usa network fornito, stazione e location wildcard)
        st = client.get_waveforms(network, "*", "*", "[BHE]H[ENZ]", t0, t1, sds_type='D')

        # Unisci segmenti e RIEMPI GAP con zeri; gestisci anche overlap (merge standard)
        # method=1 = 'merge'; fill_value=0.0 => gap riempiti con 0
        st.merge(method=1, fill_value=0.0)

        # Detrend di base
        if len(st):
            st.detrend("demean")
            st.detrend("linear")

        # (Opzionale) filtro
        if freq:
            if len(freq) == 1:
                st.filter("lowpass", freq=freq[0], corners=4, zerophase=True)
                resampling_freq=2*freq[0]
            elif len(freq) == 2:
                resampling_freq=2*freq[1]
                st.filter("bandpass", freqmin=freq[0], freqmax=freq[1], corners=4, zerophase=True)
            st.resample(float(resampling_freq), no_filter=True)
            st.merge(method=1, fill_value=0.0)  # assicura coerenza dopo il resample

        # Trim finale ESATTO alla finestra richiesta, con padding=zeri se mancano campioni
        st.trim(tini, tini + window, pad=True, fill_value=0.0)
        st.merge(method=1, fill_value=0.0)

        # Costruisci struttura self.stream: per componente, per stazione
        self.stream = {}
        for comp in comps:
            self.stream[comp] = {}
            # seleziona canali che terminano con la componente (Z/N/E)
            for tr in st:
                if tr.stats.channel.endswith(comp):
                    sta = tr.stats.station
                    # starttime come oggetto datetime/UTCDateTime coerente
                    dtime = tr.stats.starttime.datetime  # o tr.stats.starttime se preferisci UTCDateTime
                    delta = tr.stats.delta
                    data = tr.data  # numpy array
                    self.stream[comp][sta] = [dtime, delta, data]


class WaveformLoadingError(Exception):
    pass
