"""
Created on Tue Oct  4 08:14:59 2022

@author: 3Q
@modified by: Irina
"""

from sgp4.api import Satrec
from datetime import datetime
from sgp4.api import jday
from astropy.coordinates import TEME, CartesianDifferential, CartesianRepresentation
from astropy import units as u
from astropy.units import Quantity
from astropy.time import Time
from astropy.coordinates import ITRS, WGS72GeodeticRepresentation
from astropy.coordinates import EarthLocation, AltAz, Angle

# TLE for satellite GPS BIIRM-2 (PRN 31) 
#s = '1 29486U 06042A   22277.55622356 -.00000022  00000+0  00000+0 0  9998'
#t = '2 29486  54.7051 204.9391 0104917  23.9319 143.2008  2.00572655117289'
#satellite = Satrec.twoline2rv(s,t)
#jd, fr = jday(2022,9,20,13,0,0)
#error, teme_position, teme_velocity = satellite.sgp4(jd,fr)

# transform teme_position vector to True Equator Mean Equinox object
#teme_position = CartesianRepresentation(teme_position*u.km)
#teme_velocity = CartesianDifferential(teme_velocity*u.km/u.s)
#teme = TEME(teme_position.with_differentials(teme_velocity),obstime=Time(jd,format='jd'))

# transform teme position to itrs geocentric coordinates
#itrs = teme.transform_to(ITRS(obstime=Time(jd,format='jd')))
#location = itrs.earth_location
#print(location.geodetic)


#lon = Angle(str(location.geodetic.lon.value) + "d")
#lat = Angle(str(location.geodetic.lat.value) + "d")
#height = (location.geodetic.height.value*1000)*u.meter
#wgs84 = WGS72GeodeticRepresentation(lon,lat,height)
#print(wgs84)
 
'''
GPS BIIR-2  (PRN 13)    
1 24876U 97035A   22276.55443887  .00000010  00000+0  00000+0 0  9993
2 24876  55.5221 151.4562 0063300  52.6528 308.0033  2.00563166184822
'''

'''
GPS BIIRM-2 (PRN 31)    
1 29486U 06042A   22277.55622356 -.00000022  00000+0  00000+0 0  9998
2 29486  54.7051 204.9391 0104917  23.9319 143.2008  2.00572655117289
'''

file=open('Galileo_sat.txt', 'r')
l=file.readline()
#n=l
#l = file.readline()
#s=l
#l= file.readline()
#t=l
c=0
#print(n)
#print(s)
#print(t)
while True:
    l=file.readline()
    if l[0]=='1':
        s=l
    elif l[0]=='2':
        t=l
    else:
        n=l
    c=c+1
    if c%3==0:
        satellite = Satrec.twoline2rv(s,t)
        jd, fr = jday(2022,9,20,13,0,0)
        error, teme_position, teme_velocity = satellite.sgp4(jd,fr)
        # transform teme_position vector to True Equator Mean Equinox object
        teme_position = CartesianRepresentation(teme_position*u.km)
        teme_velocity = CartesianDifferential(teme_velocity*u.km/u.s)
        teme = TEME(teme_position.with_differentials(teme_velocity),obstime=Time(jd,format='jd'))
        # transform teme position to itrs geocentric coordinates
        itrs = teme.transform_to(ITRS(obstime=Time(jd,format='jd')))
        location = itrs.earth_location
        print(location.geodetic)
