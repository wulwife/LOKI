import numpy as num

_deg2rad = num.pi / 180.0
_rad2deg = 180.0 / num.pi
_km2m = 1000
_semi_major_axis = 6378137.000000
_semi_minor_axis = 6356752.314245
_eccentricity_squared = 1 - (_semi_minor_axis/_semi_major_axis)**2

class Coordinates:
    ''' 
    With this class you can convert geographical coordinates (Latitude, Longitude, Elevation)
    to generic cartesian coordinates East, North, Up (in meters). In ordert to use it you first need to define
    a reference point (Latref,Lonref,Eleref). All the calculations are performed using WGS84 ellipsoid.

    IMPORTANT: 
    For Geographical Coordinates
    Latitude is in degrees e.g. 36.117
    Longitude is in degrees e.g. -117.854
    Elevation is kilometers   
    For E,N,U Cartesian Coordinates
    East, North and Up are in meters

    USAGE:
    You generate an object coordinates as: 
    
    region=latlon2cart.Coordinates(latref,lonref,eleref)
    
    Then any other point (Lat_i, Lon_i, Ele_i) can be converted in the cartesian E,N,U frame as:

    E_i,N_i,U_i = region.geo2cart(Lat_i,Lon_i,Ele_i)

    To go back to the geographical coordinate system from the cartesian fram (E_i, N_i, U_i) you can use:

    Lat_I,Lon_I, Ele_I=region.cart2geo(E_i,N_i,U_i)

    AUTHOR: 
    Francesco Grigoli, 
    Department of Earth Sciences, 
    University of Pisa,
    Italy
    '''

    def __init__(self,lat0,lon0,ele0=0):
        X0,Y0,Z0=self.geo2cart(lat0,lon0,ele0,geo2enu=False)
        self.lat0=lat0
        self.lon0=lon0
        self.ele0=ele0
        self.X0=X0
        self.Y0=Y0
        self.Z0=Z0
        
        
    def geo2cart(self,lat,lon,ele=0,geo2enu=True):
        '''Conversion from Geographical LAT,LON,ELE(km) to Cartesian E,N,U (output in meters) frame'''

        lat*=_deg2rad
        lon*=_deg2rad
        ele*=_km2m

        N=_semi_major_axis/num.sqrt(1-_eccentricity_squared*(num.sin(lat)**2))
        
        X=(N+ele)*num.cos(lat)*num.cos(lon)
        Y=(N+ele)*num.cos(lat)*num.sin(lon)
        Z=((1-_eccentricity_squared)*N+ele)*num.sin(lat)

        if geo2enu:
            E,N,U=self.__conv2enu(X,Y,Z)
            return E,N,U
        else:
            return X,Y,Z


    def __conv2enu(self,X,Y,Z):

        lon0rad=self.lon0*_deg2rad
        lat0rad=self.lat0*_deg2rad

        DX=X-self.X0
        DY=Y-self.Y0
        DZ=Z-self.Z0

        E=DY*num.cos(lon0rad)-DX*num.sin(lon0rad)
        N=DZ*num.cos(lat0rad)-DY*num.sin(lat0rad)*num.sin(lon0rad)-DX*num.sin(lat0rad)*num.cos(lon0rad)
        U=DZ*num.sin(lat0rad)+DY*num.cos(lat0rad)*num.sin(lon0rad)+DX*num.cos(lat0rad)*num.cos(lon0rad)

        return E,N,U

    def __enu2geo(self,E,N,U):

        lon0rad=self.lon0*_deg2rad
        lat0rad=self.lat0*_deg2rad

        X=(U*num.cos(lat0rad)*num.cos(lon0rad)-E*num.sin(lon0rad)-N*num.sin(lat0rad)*num.cos(lon0rad))+self.X0
        Y=(E*num.cos(lon0rad)-N*num.sin(lat0rad)*num.sin(lon0rad)+U*num.cos(lat0rad)*num.sin(lon0rad))+self.Y0
        Z=(N*num.cos(lat0rad)+U*num.sin(lat0rad))+self.Z0

        return X,Y,Z

    def cart2geo(self,E,N,U):
        '''Conversion from Cartesian E,N,U (input in meters) to Geographical LAT,LON,ELE(km) frame'''
        X,Y,Z=self.__enu2geo(E,N,U)
        e=(_semi_major_axis**2-_semi_minor_axis**2)/_semi_minor_axis**2
        p=num.sqrt(X**2+Y**2)
        F=54*(_semi_minor_axis**2)*Z**2
        G=p**2+(1-_eccentricity_squared)*Z**2-_eccentricity_squared*(_semi_major_axis**2-_semi_minor_axis**2)
        c=((_eccentricity_squared**2)*F*p**2)/G**3
        s=num.cbrt(1+c+num.sqrt(c**2+2*c))
        k=s+1+(1/s)
        P=F/(3*k**2*G**2)
        Q=num.sqrt(1+2*(_eccentricity_squared**2)*P)
        r0_1=-((P*_eccentricity_squared*p)/(1+Q))
        r0_2=0.5*_semi_major_axis**2*(1+(1/Q))
        r0_3=(P*(1-_eccentricity_squared)*Z**2)/(Q*(1+Q))
        r0_4=0.5*P*p**2
        r0=r0_1+num.sqrt(r0_2-r0_3-r0_4)
        U=num.sqrt((p-_eccentricity_squared*r0)**2+Z**2)
        V=num.sqrt((p-_eccentricity_squared*r0)**2+(1-_eccentricity_squared)*Z**2)
        z0=(Z*_semi_minor_axis**2)/(V*_semi_major_axis)

        lat=num.arctan((Z+e*z0)/p)*_rad2deg
        lon=num.arctan2(Y,X)*_rad2deg
        ele=(U*(1-(_semi_minor_axis**2)/(V*_semi_major_axis)))/_km2m

        return lat,lon,ele



if __name__ == '__main__':
    latref=36.117; lonref=-117.854; eleref=0.
    test=Coordinates(latref,lonref,eleref)
    lat1=35.536; lon1=-118.140; ele1=0.
    E,N,U = test.geo2cart(lat1,lon1,ele1)
    lat2,lon2,ele2=test.cart2geo(E,N,U)
    print('Lat %5.3f %5.3f, Lon %5.3f %5.3f, Ele %2.1f %2.1f' %(lat1,lat2,lon1,lon2,ele1,ele2))