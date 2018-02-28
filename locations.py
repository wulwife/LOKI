class ClassName(object):
    """docstring for ."""
    def __init__(self, arg):
        super(, self).__init__()
        self.arg = arg



    def location(self, extension='*', comp=['E','N','Z'], precision='single', *input):
        nshortmin=input[0]; nshortmax=input[1]; slrat=input[2]
        npr=input[3]
        ntrial=input[4]
        traveldb=traveltimes(self.db_path, self.hdr_filename)
        tp=traveldb.load_traveltimes('P', precision)
        ts=traveldb.load_traveltimes('S', precision)
        for event_path in self.data_tree:
            loc=waveforms(event_path, extension, comp)
            event=loc.evid
            print('accessing to the event folder: ', event_path, event)
            if os.path.isdir(self.output_path+'/'+event):
               continue
            else:
               os.mkdir(self.output_path+'/'+event)
            loc.process_data(traveldb.db_stations, epsilon=0.001)
            tp_mod, ts_mod=self.time_extractor(tp, ts, loc.data_stations, traveldb.db_stations)
            tp_mod, ts_mod=tt_processing.tt_f2i(loc.deltat,tp_mod,ts_mod, npr)
            for i in range(ntrial):
                gamma=num.random.random_sample()
                nshort=num.round(nshortmin(1.-gamma*nshortmax)
                nlong=nshort*slrat
                obs_dataP, obs_dataS=loc.recstalta(nshort, nlong) #da cambiare
                corrmatrix=location.stacking(tp_mod, ts_mod, obs_dataP, obs_dataS, npr)
                cmax=num.max(corrmatrix)
                corrmatrix=num.reshape(corrmatrix,(traveldb.nx,traveldb.ny,traveldb.nz))
                (ixloc,iyloc,izloc)=num.unravel_index(num.argmax(corrmatrix),(traveldb.nx,traveldb.ny,traveldb.nz))
                #xloc=ixloc*traveldb.dx; yloc=iyloc*traveldb.dy; zloc=izloc*traveldb.dz
                xloc=traveldb.x[ixloc]; yloc=traveldb.y[iyloc]; zloc=traveldb.x[izloc]
                out_file = open(self.output_path+'/'+event+'/'+event+'.loc','a')
                out_file.write(str(i)+' '+str(xloc)+' '+str(yloc)+' '+str(zloc)+' '+str(cmax)+' '+str(nshort)+' '+str(nlong)+'\n')
                out_file.close()
                num.save(self.output_path+'/'+event+'/'+'corrmatrix_trial_'+str(ntrial),corrmatrix)
                self.coherence_plot(self.output_path+'/'+event, corrmatrix, traveldb.x, traveldb.y, traveldb.z, i)
            self.catalogue_creation(event, traveldb.lat0, traveldb.lon0, ntrial)
        print('Ho finito!!!')
