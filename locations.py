import Loki

class Location:
    """docstring for Locations"""

    def __init__(self, data_path, output_path, db_path, hdr_filename):
        self.data_path=data_path
        self.output_path=output_path
        self.db_path=db_path
        self.hdr_filename=hdr_filename
        self.data_tree, self.events=self.data_struct(self.data_path, self.output_path)

    def data_struct(self, data_path, output_path):
        events=[]
        data_tree=[]
        for root,dirs,files in os.walk(data_path):
           if not dirs:
              data_tree.append(root)
              events.append(root.split('/')[-1])
        for event in events:
           if not os.path.isdir(output_path):
              os.mkdir(output_path)
        return data_tree, events

    def location(self, extension='*', comp=['E','N','Z'], precision='single', **inputs):
        nshortmin=inputs['nshortmin']; nshortmax=inputs['nshortmax']; slrat=inputs['slrat']
        npr=inputs['npr']
        ntrial=inputs['ntrial']
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



    def loc_stalta(self, nshort_p, nshort_s, slrat, thres=2):
        tshort_p=nshort_p*self.deltat; tshort_s=nshort_s*self.deltat
        tlong_p=tshort_p*slrat; tlong_s=tshort_s*slrat
        ks_p=self.deltat/tshort_p; kl_p=self.deltat/tlong_p;
        ks_s=self.deltat/tshort_s; kl_s=self.deltat/tlong_s;
        obs_dataP=LOC_STALTA.recursive_stalta(tshort_p, tlong_p, self.deltat, self.obs_dataV, kl_p, ks_p, thres)
        obs_dataS=LOC_STALTA.recursive_stalta(tshort_s, tlong_s, self.deltat, self.obs_dataH, kl_s, ks_s, thres)
        return obs_dataP, obs_dataS
