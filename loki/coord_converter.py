from numpy import pi, sin, cos, tan, sqrt, arange

_deg2rad = pi / 180.0
_rad2deg = 180.0 / pi

_EquatorialRadius = 6378137  # equatorial radius in wgs84
_eccentricitySquared = 0.00669438  # square of eccentricity in wgs84


def LLtoUTM(Lat, Long):
    """converts lat/long to UTM coords.
    East Longitudes are positive, West longitudes are negative.
    North latitudes are positive, South latitudes are negative
    Lat and Long are in decimal degrees
    """

    a = _EquatorialRadius
    eccSquared = _eccentricitySquared
    k0 = 0.9996

    Long = (Long+180)-int((Long+180)/360)*360-180 # -180.00 .. 179.9

    LatRad = Lat*_deg2rad
    LongRad = Long*_deg2rad

    ZoneNumber, ZoneLetter = zone_identification(Lat, Long)

    LongOrigin =(ZoneNumber - 1)*6 - 180 + 3 #+3 puts origin in middle of zone
    LongOriginRad = LongOrigin * _deg2rad

    # compute the UTM Zone from the latitude and longitude
    UTMZone = "%d%c" % (ZoneNumber, ZoneLetter)

    eccPrimeSquared = (eccSquared)/(1-eccSquared)
    N = a/sqrt(1-eccSquared*sin(LatRad)*sin(LatRad))
    T = tan(LatRad)*tan(LatRad)
    C = eccPrimeSquared*cos(LatRad)*cos(LatRad)
    A = cos(LatRad)*(LongRad-LongOriginRad)

    M = a*((1
            - eccSquared/4
            - 3*eccSquared*eccSquared/64
            - 5*eccSquared*eccSquared*eccSquared/256)*LatRad
           - (3*eccSquared/8
              + 3*eccSquared*eccSquared/32
              + 45*eccSquared*eccSquared*eccSquared/1024)*sin(2*LatRad)
           + (15*eccSquared*eccSquared/256 + 45*eccSquared*eccSquared*eccSquared/1024)*sin(4*LatRad)
           - (35*eccSquared*eccSquared*eccSquared/3072)*sin(6*LatRad))

    UTMEasting = (k0*N*(A+(1-T+C)*A*A*A/6
                        + (5-18*T+T*T+72*C-58*eccPrimeSquared)*A*A*A*A*A/120)
                  + 500000.0)

    UTMNorthing = (k0*(M+N*tan(LatRad)*(A*A/2+(5-T+9*C+4*C*C)*A*A*A*A/24
                                        + (61
                                           -58*T
                                           +T*T
                                           +600*C
                                           -330*eccPrimeSquared)*A*A*A*A*A*A/720)))

    if Lat < 0:
        UTMNorthing = UTMNorthing + 10000000.0; #10000000 meter offset for southern hemisphere
    return (UTMZone, UTMEasting, UTMNorthing)


def zone_identification(Lat, Long):

    letters=['C','D','E','F','G','H','J','K','L','M','N','P','Q','R','S','T','U','V','W','X','Z']
    lats=arange(-84,96,12)
    for l in lats:
        if Lat < l:
            # da modificare
            pass

    zone_number=None
    if Lat <= 84 and Lat >= 72:
       zone_letter='X'
       if   Long >= 0  and Long <  9: zone_number = 31
       elif Long >= 9  and Long < 21: zone_number = 33
       elif Long >= 21 and Long < 33: zone_number = 35
       elif Long >= 33 and Long < 42: zone_number = 37
    elif Lat < 72 and Lat >= 64:
       zone_letter='W'
    elif Lat < 64 and Lat >= 56:
       zone_letter='V'
       if Long > 3 and Long < 12: zone_number = 32
    elif Lat < 56 and Lat >= 48:
       zone_letter='U'
    elif Lat < 48 and Lat >= 40:
       zone_letter='T'
    elif Lat < 40 and Lat >= 32:
       zone_letter='S'
    elif Lat < 32 and Lat >= 24:
       zone_letter='R'
    elif Lat < 24 and Lat >= 16:
       zone_letter='Q'
    elif Lat < 16 and Lat >= 8:
       zone_letter='P'
    elif Lat < 8 and Lat >= 0:
       zone_letter='N'
    elif Lat < 0 and Lat >= -8:
       zone_letter='M'
    elif Lat < -8 and Lat >= -16:
       zone_letter='L'
    elif Lat < -16 and Lat >= -24:
       zone_letter='K'
    elif Lat < -24 and Lat >= -32:
       zone_letter='J'
    elif Lat < -32 and Lat >= -40:
       zone_letter='H'
    elif Lat < -40 and Lat >= -48:
       zone_letter='G'
    elif Lat <= -48 and Lat >= -56:
       zone_letter='F'
    elif Lat <= -56 and Lat >= -64:
       zone_letter='E'
    elif Lat <= -64 and Lat >= -72:
       zone_letter='D'
    elif Lat <= -72 and Lat >= -80:
       zone_letter='C'
    else:
       zone_letter='Z'

    if zone_number==None:
       zone_number = int((Long + 180)/6) + 1

    return zone_number, zone_letter


def UTMtoLL(northing, easting, zone):
    """converts UTM coords to lat/long.  Equations from USGS Bulletin 1532
    East Longitudes are positive, West longitudes are negative.
    North latitudes are positive, South latitudes are negative
    Lat and Long are in decimal degrees.
    Written by Chuck Gantz- chuck.gantz@globalstar.com
    Converted to Python by Russ Nelson <nelson@crynwr.com>
    Modified by Francesco Grigoli <fsco.grigoli@gmail.com>
    """

    k0 = 0.9996
    a = _EquatorialRadius
    eccSquared = _eccentricitySquared
    e1 = (1-sqrt(1-eccSquared))/(1+sqrt(1-eccSquared))
    #NorthernHemisphere; //1 for northern hemispher, 0 for southern

    x = easting - 500000.0 #remove 500,000 meter offset for longitude
    y = northing

    ZoneLetter = zone[-1]
    ZoneNumber = int(zone[:-1])
    if ZoneLetter >= 'N':
        NorthernHemisphere = 1  # point is in northern hemisphere
    else:
        NorthernHemisphere = 0  # point is in southern hemisphere
        y -= 10000000.0         # remove 10,000,000 meter offset used for southern hemisphere

    LongOrigin = (ZoneNumber - 1)*6 - 180 + 3  # +3 puts origin in middle of zone

    eccPrimeSquared = (eccSquared)/(1-eccSquared)

    M = y / k0
    mu = M/(a*(1-eccSquared/4-3*eccSquared*eccSquared/64-5*eccSquared*eccSquared*eccSquared/256))

    phi1Rad = (mu + (3*e1/2-27*e1*e1*e1/32)*sin(2*mu)
               + (21*e1*e1/16-55*e1*e1*e1*e1/32)*sin(4*mu)
               +(151*e1*e1*e1/96)*sin(6*mu))
    phi1 = phi1Rad*_rad2deg;

    N1 = a/sqrt(1-eccSquared*sin(phi1Rad)*sin(phi1Rad))
    T1 = tan(phi1Rad)*tan(phi1Rad)
    C1 = eccPrimeSquared*cos(phi1Rad)*cos(phi1Rad)
    R1 = a*(1-eccSquared)/pow(1-eccSquared*sin(phi1Rad)*sin(phi1Rad), 1.5)
    D = x/(N1*k0)

    Lat = phi1Rad - (N1*tan(phi1Rad)/R1)*(D*D/2-(5+3*T1+10*C1-4*C1*C1-9*eccPrimeSquared)*D*D*D*D/24
                                          +(61+90*T1+298*C1+45*T1*T1-252*eccPrimeSquared-3*C1*C1)*D*D*D*D*D*D/720)
    Lat = Lat * _rad2deg

    Long = (D-(1+2*T1+C1)*D*D*D/6+(5-2*C1+28*T1-3*C1*C1+8*eccPrimeSquared+24*T1*T1)
            *D*D*D*D*D/120)/cos(phi1Rad)
    Long = LongOrigin + Long * _rad2deg
    return (Lat, Long)


if __name__ == '__main__':
    lon0=11.2003
    lat0=45.5003
    (z, e, n) = LLtoUTM(lon0, lat0)
    (lon1,lat1)=UTMtoLL(n, e, z)
    print(lon0,' : ',lon1,' and ', lat0,' : ',lat1)
