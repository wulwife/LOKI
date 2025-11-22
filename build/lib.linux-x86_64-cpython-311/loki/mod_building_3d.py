import loki.LatLongUTMconversion as ll
from scipy.interpolate import griddata
import numpy as num
import pylab as plt


def conversion(latmin, latmax, lonmin, lonmax, border, spacing, filename):
    zone,xmin,ymin=ll.LLtoUTM(23,latmin,lonmin)
    zone,xmax,ymax=ll.LLtoUTM(23,latmax,lonmax)
    nx=int((xmax-xmin)/spacing)
    ny=int((ymax-ymin)/spacing)
    xax=num.arange(0,nx)*spacing
    yax=num.arange(0,ny)*spacing
    f=open(filename,'r')
    points=[]
    for line in f:
        toks=line.split()
        lonobs,latobs,elvobs=eval(toks[0]),eval(toks[1]),eval(toks[2])
        if (latobs>latmin-border and latobs<latmax+border) and (lonobs>lonmin-border and lonobs<lonmax+border):
           zone,xobs,yobs=ll.LLtoUTM(23,latobs,lonobs)
           points.append([xobs-xmin,yobs-ymin,elvobs])
    points=num.array(points)
    return points, xax, yax


def interpolation(points, xax, yax):
    xgrid, ygrid = num.meshgrid(xax, yax)
    zgrid = griddata(points[:,0:2],points[:,2],(xgrid, ygrid), method='nearest')
    #plt.imshow(zgrid)
    #plt.show()
    return zgrid

filename='piton_topo'
points,xax,yax=conversion(-21.25, -21.21, 55.69, 55.74, 0.1, 50, filename)
topo=interpolation(points, xax, yax)



dx=0.05 # spacing in km
nx=len(xax); ny=len(yax); nz=450 # numero grid points lungo x,y,z
pmodel=num.zeros([nx,ny,nz], dtype=num.float32)
smodel=num.zeros([nx,ny,nz], dtype=num.float32)
vp1=3.5; vp2=3.7; vp3=4.2; vp4=4.5; vp5=5.0; vp6=6.5; vp7=6.8; vp8=8.0;
top1=2700; top2=800; top3=-790; top4=-1200; top5=-2200; top6=-4200; top7=-17000; top8=-20000;
vpsratio=1.75
vs1=0.3/vpsratio; vs2=3.5/vpsratio; vs3=3.7/vpsratio; vs4=4.2/vpsratio; vs5=5.0/vpsratio; vs6=6.5/vpsratio; vs7=6.8/vpsratio; vs8=8.0/vpsratio;
iz1=int((top1-top2)/50.)+1;
iz2=iz1+int((top2-top3)/50.)+1;
iz3=iz2+int((top3-top4)/50.)+1;
iz4=iz3+int((top5-top5)/50.)+1;
iz5=iz4+int((top5-top6)/50.)+1;
iz6=iz5+int((top6-top7)/50.)+1;
iz7=iz6+int((top7-top8)/50.)+1;
print(iz1,iz2)
pmodel[:,:,0:iz1]=vp1
pmodel[:,:,iz1:iz2]=vp2
pmodel[:,:,iz2:iz3]=vp3
pmodel[:,:,iz3:iz4]=vp4
pmodel[:,:,iz4:iz5]=vp5
pmodel[:,:,iz5:iz6]=vp6
pmodel[:,:,iz6:iz7]=vp7
pmodel[:,:,iz7:]=vp8

vptopo=0.3;
top1=2700;
vstopo=0.3/vpsratio;
for i in range(nx):
	for j in range(ny):
		iztopo=int((num.round(top1-topo[i,j])/50.))+1
		pmodel[i,j,0:iztopo]=vptopo#(1./vptopo)*dx #podvin e lecompte richiede slowness*spacing
		smodel[i,j,0:iztopo]=vstopo#(1./vstopo)*dx


print(pmodel[0,0,:])
plt.figure()
plt.imshow(num.transpose(pmodel[:,50,:]))
plt.colorbar()
plt.show()

input()

f=open('layer.P.mod.buf','wb')
for i in range(nx):
	for j in range(ny):
		for k in range(nz):
			f.write(pmodel[i,j,k])
f.close()
modelp = num.fromfile('layer.P.mod.buf', dtype=num.float32).reshape((nx,ny,nz))
plt.figure()
plt.imshow(num.transpose(modelp[75,:,:]))
plt.colorbar()
#print modelp
f=open('layer.S.mod.buf','wb')
for i in range(nx):
	for j in range(ny):
		for k in range(nz):
			f.write(smodel[i,j,k])
f.close()
models = num.fromfile('layer.S.mod.buf', dtype=num.float32).reshape((nx,ny,nz))
plt.figure()
plt.imshow(num.transpose(models[75,:,:]))
plt.colorbar()
plt.show()
