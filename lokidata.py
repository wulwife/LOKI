import traveltimes
import waveforms
from datetime import datetime
import numpy as num
import matplotlib.pyplot as plt

class LokiData:

    def __init__(self, tobj, wobj, derivative=True):
        self.check_sampling_rate(wobj)
        self.check_starting_time(wobj)
        self.traces=self.loki_input(wobj, tobj, derivative)

    def check_sampling_rate(self,wobj):
        intsamp=1E6
        deltas=[]
        for comp in (wobj.stream).keys():
            for sta in (wobj.stream[comp]).keys():
                deltas.append(wobj.stream[comp][sta][1])
        deltas=num.array(deltas)
        ideltas=num.unique((deltas*intsamp).astype(int))
        if num.size(ideltas)==1:
           self.deltat=deltas[0]
        else:
           raise ValueError('Error!! All trace must have the same sampling rate')

    def check_starting_time(self,wobj):
        intsamp=1E6
        dtimes=[]
        self.ns=0
        for comp in (wobj.stream).keys():
            for sta in (wobj.stream[comp]).keys():
                dtimes.append(wobj.stream[comp][sta][0])
                if self.ns<num.size(wobj.stream[comp][sta][2]):
                   self.ns=num.size(wobj.stream[comp][sta][2])
        self.dtime_max=max(dtimes)
        self.evid=(self.dtime_max).isoformat()

    def loki_input(self, wobj, tobj, derivative):
        comp=tuple((wobj.stream).keys())
        if len(comp)==3:
            xtr=self.select_data(comp[0], wobj, tobj.db_stations, derivative)
            ytr=self.select_data(comp[1], wobj, tobj.db_stations, derivative)
            ztr=self.select_data(comp[2], wobj, tobj.db_stations, derivative)
            return (xtr, ytr, ztr)
        elif len(comp)==1:
            ztr=self.select_data(comp, wobj.data_stations, db_stations, derivative)
            return (ztr,)
        else:
            raise ValueError('Traces must have 1 or 3 components!')

    def select_data(self, comp, wobj, db_stations, derivative):
        self.stations=tuple(wobj.data_stations & db_stations)
        self.nstation=num.size(self.stations)
        tr=num.zeros([self.nstation,self.ns])
        stream=wobj.stream[comp]
        for i,sta in enumerate(self.stations):
            nstr=num.size(stream[sta][2])
            idt=num.int((self.dtime_max-stream[sta][0]).total_seconds()/self.deltat)
            tr[i,0:nstr-idt]=stream[sta][2][idt:]
            if derivative:
               tr[i,1:self.ns]=((tr[i,1:]-tr[i,0:self.ns-1])/self.deltat)
               tr[i,0]=0.
            else:
               tr[i,:]=1.
        return tr

    def time_extractor(self, tp, ts):
        nsta=len(self.stations)
        nxyz= num.size(tp[self.stations[0]])
        tp_mod=num.zeros([nxyz,nsta]) #optimize with pointers to tp[sta] and ts[sta]
        ts_mod=num.zeros([nxyz,nsta])
        for i,sta in enumerate(self.stations):
            tp_mod[:,i]=tp[sta]
            ts_mod[:,i]=ts[sta]
        return (tp_mod, ts_mod)
