import traveltimes
import waveforms
from datetime import datetime

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
                deltas.append(int(wobj.stream[comp][sta][1]))
        deltas=num.array(deltas)
        ideltas=num.unique((deltas*intsamp).astype(int))
        if num.size(ideltas)==1:
           self.deltat=deltas[0]
        else:
           raise ValueError('Error!! All trace must have the same sampling rate')

    def check_starting_time(self,wobj):
        intsamp=1E6
        dtimes=[]
        for comp in (wobj.stream).keys():
            for sta in (wobj.stream[comp]).keys():
                dtimes.append(wobj.stream[comp][sta][0])
        dtime0=dtime[0]
        for dtime in dtimes:
            if num.abs((dtime0-dtime).total_seconds())<2*wobj.deltat:
               pass
            else:
               raise ValueError('Error!! All trace must have the same starting time')
        self.evid=dtime0.isoformat()

    def loki_input(self, wobj, tobj, derivative):
        comp=tuple((wobj.stream).keys())
        if len(comp)==3:
            xtr=self.select_data(comp[0], wobj.data_stations, tobj.db_stations, derivative)
            ytr=self.select_data(comp[1], wobj.data_stations, tobj.db_stations, derivative)
            ztr=self.select_data(comp[2], wobj.data_stations, tobj.db_stations, derivative)
            return (xtr, ytr, ztr)
        elif len(comp)==1:
            ztr=self.select_data(comp, wobj.data_stations, db_stations, derivative)
            return (ztr,)
        else:
            raise ValueError('Traces must have 1 or 3 components!')

    def select_data(self, comp, data_stations, db_stations, derivative):
        self.stations=tuple(data_stations & db_stations)
        self.nstation=num.size(self.stations)
        tr=num.zeros([self.nstation,self.ns])
        stream=self.streams[comp]
        for i,sta in enumerate(self.stations):
            if sta in stream.keys():
               nstr=num.size(stream[sta])
               tr[i,0:nstr]=stream[sta][1]
               if derivative:
                  tr[i,0:nstr-1]=((tr[i,1:]-tr[i,0:nstr-1])/self.deltat)
                  tr[i,nstr-1]=0.
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
