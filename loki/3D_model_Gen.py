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
    return zgrid, xgrid, ygrid

filename='piton_topo'
spacing=100.
latmin=-21.35; latmax=-21.10; lonmin=55.60; lonmax=55.85;
points,xax,yax=conversion(-21.35, -21.10, 55.60, 55.85, 0.1, 100, filename)
topo, xgrid, ygrid=interpolation(points, xax, yax)
topo=topo-100.
zone,xmin,ymin=ll.LLtoUTM(23,latmin,lonmin)
latsta=-21.33533333; lonsta=55.794
zone,xsta,ysta=ll.LLtoUTM(23,latsta,lonsta)
ix=int((xsta-xmin)/spacing)
iy=int((ysta-ymin)/spacing)







#fig, ax = plt.subplots(subplot_kw=dict(projection='3d'))

#ls = LightSource(270, 45)
# To use a custom hillshading mode, override the built-in shading and pass
# in the rgb colors of the shaded surface calculated from "shade".
#rgb = ls.shade(topo, cmap=cm.gist_earth, vert_exag=0.1, blend_mode='soft')
#surf = ax.plot_surface(xgrid, ygrid, topo, rstride=1, cstride=1, facecolors=rgb,
#                      linewidth=0, antialiased=False, shade=False)

#plt.show()







dx=0.1 # spacing in km
nx=len(xax); ny=len(yax); nz=250 # numero grid points lungo x,y,z
print('nx ny nz', nx, ny ,nz)
pmodel=num.zeros([nx,ny,nz], dtype=num.float32)
smodel=num.zeros([nx,ny,nz], dtype=num.float32)
vp1=3.5; vp2=3.7; vp3=4.2; vp4=4.5; vp5=5.0; vp6=6.5; vp7=6.8; vp8=8.0;
top1=2700; top2=800; top3=-790; top4=-1200; top5=-2200; top6=-4200; top7=-17000; top8=-20000;
vpsratio=1.75
vs1=vp1/vpsratio; vs2=vp2/vpsratio; vs3=vp3/vpsratio; vs4=vp4/vpsratio; vs5=vp5/vpsratio; vs6=vp6/vpsratio; vs7=vp7/vpsratio; vs8=vp8/vpsratio;
iz1=int((top1-top2)/spacing)+1;
iz2=iz1+int((top2-top3)/spacing)+1;
iz3=iz2+int((top3-top4)/spacing)+1;
iz4=iz3+int((top5-top5)/spacing)+1;
iz5=iz4+int((top5-top6)/spacing)+1;
iz6=iz5+int((top6-top7)/spacing)+1;
iz7=iz6+int((top7-top8)/spacing)+1;

pmodel[:,:,0:iz1]=(1./vp1)*dx
pmodel[:,:,iz1:iz2]=(1./vp2)*dx
pmodel[:,:,iz2:iz3]=(1./vp3)*dx
pmodel[:,:,iz3:iz4]=(1./vp4)*dx
pmodel[:,:,iz4:iz5]=(1./vp5)*dx
pmodel[:,:,iz5:iz6]=(1./vp6)*dx
pmodel[:,:,iz6:iz7]=(1./vp7)*dx
pmodel[:,:,iz7:]=(1./vp8)*dx


smodel[:,:,0:iz1]=(1./vs1)*dx
smodel[:,:,iz1:iz2]=(1./vs2)*dx
smodel[:,:,iz2:iz3]=(1./vs3)*dx
smodel[:,:,iz3:iz4]=(1./vs4)*dx
smodel[:,:,iz4:iz5]=(1./vs5)*dx
smodel[:,:,iz5:iz6]=(1./vs6)*dx
smodel[:,:,iz6:iz7]=(1./vs7)*dx
smodel[:,:,iz7:]=(1./vs8)*dx


vptopo=0.3;
top1=2700;
vstopo=vptopo/vpsratio;
print(vstopo, vs1, vs2, vs3, vs4, vs5, vs5)
for i in range(nx):
	for j in range(ny):
		iztopo=int((num.round(top1-topo[j,i])/spacing))+1
		pmodel[i,j,0:iztopo]=(1./vptopo)*dx#(1./vptopo)*dx #podvin e lecompte richiede slowness*spacing
		smodel[i,j,0:iztopo]=(1./vstopo)*dx#(1./vstopo)*dx
		#iz1=iztopo+int((top1-top2)/spacing)+1
		#pmodel[i,j,iztopo:iz1]=(1./vp1)*dx#(1./vptopo)*dx #podvin e lecompte richiede slowness*spacing
		#smodel[i,j,iztopo:iz1]=(1./vs1)*dx#(1./vstopo)*dx
		#iz2=iz1+int((top2-top3)/spacing)+1
		#pmodel[i,j,iz1:iz2]=(1./vp2)*dx#(1./vptopo)*dx #podvin e lecompte richiede slowness*spacing
		#smodel[i,j,iz1:iz2]=(1./vs2)*dx#(1./vstopo)*dx

print(num.shape(pmodel))

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
