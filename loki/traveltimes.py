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
import os
import sys
import numpy as num
import loki.LatLongUTMconversion as ll


class Traveltimes:

    def __init__(self, db_path, hdr_filename):
        if not os.path.isdir(db_path):
            print('Error: data or database path do not exist')
            sys.exit()
        self.db_path = db_path
        if not os.path.isfile(db_path+'/'+hdr_filename):
            print('Error: header file does not exist')
            sys.exit()
        self.hdr_filename = hdr_filename
        self.load_header()

    def load_header(self):

        f = open(os.path.join(self.db_path, self.hdr_filename))
        lines = f.readlines()
        f.close()
        self.nx, self.ny, self.nz = [ int(x)   for x in lines[0].split()]
        self.x0, self.y0, self.z0 = [ float(x) for x in lines[1].split()]
        self.dx, self.dy, self.dz = [ float(x) for x in lines[2].split()]
        self.lat0, self.lon0 = [ float(x) for x in lines[3].split()]
        self.x = self.x0+(num.arange(0,self.nx)*self.dx)
        self.y = self.y0+(num.arange(0,self.ny)*self.dy)
        self.z = self.z0+(num.arange(0,self.nz)*self.dz)
        self.nxyz=self.nx*self.ny*self.nz
        db_stalist=[]
        if len(lines[4].split())>1:
            stations_coordinates={}
        else:
            stations_coordinates=None

        for line in lines[4:]:
            toks=line.split()
            db_stalist.append(toks[0])
            if len(toks)>1:
                stations_coordinates[toks[0]]=[eval(toks[1]), eval(toks[2]), eval(toks[3])]

        self.db_stations=set(db_stalist)
        self.stations_coordinates=stations_coordinates

    def load_traveltimes(self, phase, label='layer', precision='single'):
        t={}
        for sta in self.db_stations:
            try:
               fn = os.path.join(self.db_path, '%(label)s.%(phase)s.%(station)s.time.buf' %{"label":label,"phase":phase, "station":sta} )
            except:
               print('Error: reading file for station' + sta)
               sys.exit()
            if (precision=='single'):
                t[sta]= num.fromfile(fn, dtype=num.float32)
            elif (precision=='double'):
                t[sta]= num.fromfile(fn, dtype=num.float64)
            else:
                print('Error: precision must be set to "single" or "double"!!!')
                sys.exit()
        return t

    def ttdb_reduce(self,tt,l_lim,u_lim,zlim=[]):
        (zorig,x_orig,y_orig)=ll.LLtoUTM(23, self.lat0, self.lon0)
        (zorig,x_l,y_l)=ll.LLtoUTM(23, l_lim[0], l_lim[1])
        (zorig,x_u,y_u)=ll.LLtoUTM(23, u_lim[0], u_lim[1])
        x_l=(x_l-x_orig)/1000.; y_l=(y_l-y_orig)/1000.
        x_u=(x_u-x_orig)/1000.; y_u=(y_u-y_orig)/1000.
        nx_ini=int((x_l-self.x0)/self.dx); ny_ini=int((y_l-self.y0)/self.dy)
        nx_fin=int((x_u-self.x0)/self.dx); ny_fin=int((y_u-self.y0)/self.dy)
        self.x0=x_l-self.x0; self.y0=y_l-self.y0
        if zlim:
           nz_ini=int((zlim[0]-self.z0)/self.dz)
           nz_fin=int((zlim[1]-self.z0)/self.dz)
           self.z0=zlim[0]
        else:
           nz_ini=0; nz_fin=self.nz
        nx_new=nx_fin-nx_ini; ny_new=ny_fin-ny_ini; nz_new=nz_fin-nz_ini;
        tt_new={}
        for key in tt.keys():
            t_arr=tt[key].reshape(self.nx,self.ny,self.nz)
            tt_new[key]=t_arr[nx_ini:nx_fin,ny_ini:ny_fin,nz_ini:nz_fin].reshape(nx_new*ny_new*nz_new)
        self.nx=nx_new; self.ny=ny_new; self.nz=nz_new
        self.x0=0.; self.y0=0.
        self.x = self.x0+(num.arange(0,self.nx)*self.dx)
        self.y = self.y0+(num.arange(0,self.ny)*self.dy)
        self.z = self.z0+(num.arange(0,self.nz)*self.dz)
        self.nxyz=self.nx*self.ny*self.nz
        self.lat0=l_lim[0]
        self.lon0=l_lim[1]

        return tt_new

    def interpolation(self, tt, dx_new, dy_new, dz_new):
        from scipy.interpolate import RegularGridInterpolator
        xr=self.dx/dx_new; yr=self.dy/dy_new; zr=self.dz/dz_new;
        x_old=num.arange(self.nx)*self.dx; y_old=num.arange(self.ny)*self.dy; z_old=num.arange(self.nz)*self.dz
        x_new=num.arange(1,(self.nx-1)*xr)*dx_new; y_new=num.arange(1,(self.ny-1)*yr)*dy_new; z_new=num.arange(1,(self.nz-1)*zr)*dz_new
        nx_new=num.size(x_new); ny_new=num.size(y_new); nz_new=num.size(z_new)
        grid=[]
        for i in range(nx_new):
              for j in range(ny_new):
                  for k in range(nz_new):
                        grid.append([x_new[i],y_new[j],z_new[k]])
        grid=num.array(grid)
        t_interp={}
        for key in tt.keys():
            print('interpolation station: '+key+'\n')
            t_arr=tt[key].reshape(self.nx,self.ny,self.nz)
            interpolation = RegularGridInterpolator((x_old, y_old, z_old), t_arr)
            t_int=interpolation(grid)
            t_interp[key]=t_int.reshape(nx_new*ny_new*nz_new)
        self.nx=nx_new; self.ny=ny_new; self.nz=nz_new
        self.dx=dx_new; self.dy=dy_new; self.dz=dz_new
        self.x = self.x0+x_new; self.y = self.y0+y_new; self.z = self.z0+z_new
        self.nxyz=self.nx*self.ny*self.nz

        return t_interp

    def save_ttdb(self,tt,phase,label):

        f=open(self.db_path+'/'+label+'.header.hdr','w')
        f.write('%d %d %d' %(self.nx,self.ny,self.nz)+'\n'+
                '%f %f %f' %(self.x0, self.y0, self.z0)+'\n'+
                '%f %f %f' %(self.dx, self.dy, self.dz)+'\n'+
                '%f %f' %(self.lat0, self.lon0)+'\n')

        for key in tt.keys():
            fn = os.path.join(self.db_path, '%(label)s.%(phase)s.%(station)s.time.buf' %{"label":label, "phase":phase, "station":key})
            (tt[key].astype(dtype=num.float32)).tofile(fn)
            f.write(key+'\n')

        f.close()

        return None

    def event_indexes(self,evlat,evlon,evdepth):
        zonei,xi,yi=ll.LLtoUTM(23,evlat,evlon)
        zorig,xorig,yorig=ll.LLtoUTM(23,self.lat0,self.lon0)
        zi=evdepth*km
        x0=xorig+self.x0*km; y0=yorig+self.y0*km; z0=0.0+self.z0*km
        d_spac=self.dx*km
        ix=int(num.round((xi-x0)/d_spac)); iy=int(num.round((yi-y0)/d_spac)); iz=int(num.round((zi-z0)/d_spac))
        return ix,iy,iz

    def ttdb_generator(self, velocity, phase='P'):
        # Generate traveltimes for an homogenaous medium
        if self.stations_coordinates is not None:
            for sta in self.stations_coordinates:
                print('Starting to calculate traveltimes for the station : ' + sta +' \n')
                xsta,ysta,zsta=self.stations_coordinates[sta][0:3]
                fname='homo.%(phase)s.%(station)s.time.buf' %{"phase":phase, "station":sta}
                fout=open(self.db_path+'/'+fname,'wb')
                tt=num.zeros(self.nxyz)
                print(self.nxyz,self.nx,self.ny,self.nz)
                for k in range(self.nxyz):
                    ix=k//(self.ny*self.nz)
                    iy=k//self.nz-(ix*self.ny)
                    iz=k-(iy*self.nz)-(ix*self.ny*self.nz)
                    dist=num.sqrt((self.x[ix]-xsta)**2+(self.y[iy]-ysta)**2+(self.z[iz]-zsta)**2)
                    tt[k]=dist/velocity
                tt.tofile(fout)
                fout.close()
                print('Traveltimes computation for the station : ' + sta + 'completed! \n')
        return None

    def apply_master_event_correction(self, phase, dt, label='layer', precision='single'):
        t={}
        for sta in dt.keys():
            try:
              fn = os.path.join(self.db_path, '%(label)s.%(phase)s.%(station)s.time.buf' %{"label":label,"phase":phase, "station":sta} )
              if (precision=='single'):
                  tt= num.fromfile(fn, dtype=num.float32)+dt[sta]
                  tt[tt<0.]=0.
                  t[sta]=tt
              elif (precision=='double'):
                  t[sta]= num.fromfile(fn, dtype=num.float64)+dt[sta]
              else:
                  print('Error: precision must be set to "single" or "double"!!!')
                  sys.exit()
            except:
                print('Error: reading file for station' + sta)
        return t

# db_path='/Users/francesco/Desktop/KOREA/korea_time'
# tt0=traveltimes(db_path, 'header.hdr')
# tt1=traveltimes(db_path, 'header.hdr')
# tp0=tt0.load_traveltimes('P','layer')
# tp1=tt1.load_traveltimes('P','layer')
# tp_red0=tt0.ttdb_reduce(tp0,[35.5,129.00],[36.5,130.0])
# tp_red1=tt1.ttdb_reduce(tp1,[35.5,129.00],[36.5,130.0])
# tp_int1=tt1.interpolation(tp_red1,0.5,0.5,0.5)
# tt3=traveltimes(db_path, 'interpolated.header.hdr')
# tp3=tt3.load_traveltimes('P','interpolated')
# #tt1.save_ttdb(tp_int1,'P','interpolated')
#
# #400 350 65
# #-250.000000 -150.000000 -5.000000
# #31.88 31.94  130.84 130.92
# #for k in tp0.keys():
# #    #tp1=tp[k].reshape(400, 350, 65)
# #    tp2=tp_red0[k].reshape(tt0.nx, tt0.ny, tt0.nz)
# #    tp4=tp3[k]
# #    print(num.shape(tp4),tt3.nx*tt3.ny*tt3.nz)
# #    plt.figure()
# #    plt.imshow(tp2[:,:,0])
# #    plt.figure()
# #    plt.imshow(tp3[:,:,0])
# #    plt.show()
