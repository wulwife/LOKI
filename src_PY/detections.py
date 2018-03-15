class ClassName(object):
    """docstring for ."""
    def __init__(self, arg):
        super(, self).__init__()
        self.arg = arg

    def detection(self, extension='*', comp=['E','N','Z'], precision='single', *inputs):
        nshort_p=inputs[0]; nshort_s=inputs[1]; slrat=inputs[2]
        npr=inputs[3]
        tstep=inputs[4]
        toverlap=inputs[5]
        traveldb=traveltimes(self.db_path, self.hdr_filename)
        tp=traveldb.load_traveltimes('P', precision)
        ts=traveldb.load_traveltimes('S', precision)
        for event_path in self.data_tree:
            print(event_path)
            loc=waveforms(event_path, extension, comp)
            #event=loc.evid
            event=datetime.datetime.strptime(loc.evid,"%Y-%m-%dT%H:%M:%S.%fZ")
            tp_mod, ts_mod=self.time_extractor(tp, ts, loc.data_stations, traveldb.db_stations)
            loc.process_data(traveldb.db_stations, epsilon=0.001)
            obs_dataP, obs_dataS=loc.recstalta(nshort_p, nshort_s, slrat)
            nsamples=num.size(obs_dataP[0,:])
            noverlap=int(toverlap/loc.deltat)#200
            step=int(tstep/loc.deltat)#500
            tp_mod, ts_mod=tt_processing.tt_f2i(loc.deltat,tp_mod,ts_mod, npr)
            if not os.path.isdir(self.output_path):
               os.mkdir(self.output_path)
            kl=0
            #t0i=int(0*3600*50)
            #t0f=int(3*3600*50)
            detections=[[event,0,0,0,0]]
            for k in range(0,nsamples,step):
                initwin=k-noverlap*kl
                kl=1
                endwin=initwin+step+noverlap
                pstalta=num.zeros([loc.nstation,step+noverlap])
                sstalta=num.zeros([loc.nstation,step+noverlap])
                for i in range(loc.nstation):
                    #normfactorP=num.trapz(obs_dataP[i,initwin:endwin], dx=loc.deltat)
                    #normfactorS=num.trapz(obs_dataS[i,initwin:endwin], dx=loc.deltat)
                    normfactorP=num.max(obs_dataP[i,initwin:endwin])
                    normfactorS=num.max(obs_dataS[i,initwin:endwin])
                    sstalta[i,:]=(obs_dataS[i,initwin:endwin]/normfactorS)
                    pstalta[i,:]=(obs_dataP[i,initwin:endwin]/normfactorP)
                print('range%d-%d'%(initwin,endwin))
                #skip=loc.skip_chunk(pstalta, 5., 3)
                #print('skip ', skip)
                #if skip<15:
                #    continue
                #loc=waveforms(event_path, extension, comp)
                #event=datetime.datetime.strptime(loc.evid,"%Y-%m-%dT%H:%M:%S.%fZ")
                print('accessing to the event folder: ', event_path)
                tini=event+datetime.timedelta(seconds=initwin*loc.deltat)
                tfin=event+datetime.timedelta(seconds=(initwin+step)*loc.deltat)
                print(tini,tfin)
                tdelta=(tfin-tini).total_seconds()
                corrmatrix=ps_detection.stacking(tp_mod, ts_mod, pstalta, sstalta, step, npr)
                #corrmatrix=p_detection.stacking(tp_mod, pstalta, step, npr)
                maxloc=num.max(corrmatrix,axis=0)
                avgloc=num.mean(corrmatrix,axis=0)
                thres=0.0; coherthres=5.0;
                coheratio=maxloc/avgloc
                imaxratio=num.argmax(maxloc)
                print('coherence level', maxloc[imaxratio], 'coherence ratio', coheratio[imaxratio])
                if maxloc[imaxratio]>thres and coheratio[imaxratio]>coherthres:
                   #index=num.argmax(corrmatrix)
                   index=num.argmax(corrmatrix[:,imaxratio])
                   #itime=num.argmax(detections)
                   itime=imaxratio
                   (ixloc,iyloc,izloc)=num.unravel_index(index,(traveldb.nx,traveldb.ny,traveldb.nz))
                   xloc=traveldb.x[ixloc]; yloc=traveldb.y[iyloc]; zloc=traveldb.z[izloc]
                   late, lone = self.coordinate_conversion(xloc, yloc, traveldb.lat0, traveldb.lon0, refell=23)
                   tev=event+datetime.timedelta(seconds=(initwin+itime)*loc.deltat)
                   if num.abs((detections[-1][0]-tev).total_seconds())>120:
                      detections.append([event+datetime.timedelta(seconds=(initwin+itime)*loc.deltat),coheratio[imaxratio],late,lone,zloc])
                      f=open(self.output_path+'/raw_event_catalogue.txt','a')
                      f.write(str(event+datetime.timedelta(seconds=(initwin+itime)*loc.deltat)+datetime.timedelta(seconds=0.))+' '+str(maxloc[imaxratio])+' '+str(coheratio[imaxratio])+' '+ str(late)+' '+str(lone)+' '+str(zloc)+ '\n')
                      print('event :: ',str(event+datetime.timedelta(seconds=(initwin+itime)*loc.deltat)+datetime.timedelta(seconds=0.))+' '+str(coheratio[imaxratio])+' '+ str(late)+' '+str(lone)+' '+str(zloc)+'\n')
                      f.close()
                      #plt.imshow(corrmatrix.reshape(traveldb.nx,traveldb.ny,traveldb.nz)[:,:,5])
                      #plt.show()
            f=open(self.output_path+'/final_event_catalogue.txt','w')
            for event in detections:
                f.write(str(event[0])+' '+str(event[1])+' '+ str(event[2])+' '+str(event[3])+' '+str(event[4])+ '\n')
            f.close()
        print('Ho finito!!!')
