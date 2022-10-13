"""
Created on Tue Oct  4 08:14:59 2022

@author: 3Q
"""

from sgp4.api import Satrec
from sgp4.api import jday
from astropy.coordinates import TEME, CartesianDifferential, CartesianRepresentation
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import ITRS
from pyproj import Transformer
from tle import tle_json


def tle_to_itrs(l1,l2,day=20,month=9,year=2022,hour=12,minute=30,second=0):
    
    satellite = Satrec.twoline2rv(l1, l2)
    jd, fr = jday(year,month,day,hour,minute,second)
    error, teme_position, teme_velocity = satellite.sgp4(jd,fr)
    
    # transform teme_position vector to True Equator Mean Equinox object
    teme_position = CartesianRepresentation(teme_position*u.km)
    teme_velocity = CartesianDifferential(teme_velocity*u.km/u.s)
    teme = TEME(teme_position.with_differentials(teme_velocity),obstime=Time(jd,format='jd'))
    
    # transform teme object to itrs geocentric coordinates
    itrs = teme.transform_to(ITRS(obstime=Time(jd,format='jd')))
    position_itrs = (itrs.earth_location.geodetic.lon.value,itrs.earth_location.geodetic.lat.value,itrs.earth_location.geodetic.height.value)
    
    return position_itrs

def convert_crs(from_crs,to_crs,pos_vector):
    
    transformer =  Transformer.from_crs(from_crs,to_crs,True)
    position = transformer.transform(pos_vector[0],pos_vector[1],pos_vector[2])
    
    return position

def sat_pos():
    
    gps_path = '../data/gps.txt'
    galileo_path = '../data/galileo.txt'
    print('Reading TLE')
    satellites = [tle_json(gps_path),tle_json(galileo_path)]
    
    print('From TEME Compute Satellite Position in ETRS89/UTM32')
    for i in satellites:
        for tle in i:
            
            two_lines = list(tle.values())
            key = list(tle.keys())[0]
            
            s = two_lines[0][0]
            t = two_lines[0][1]
            #yr,mon,day,hr,minute,sec = 2022,9,20,13,0,0
        
            # get satellite position in itrf2014
            position_itrs_2k14 = tle_to_itrs(s,t)
            
            # get satellite position in itrf2000
            position_itrs_2k = convert_crs(7912, 7909, position_itrs_2k14)
            
            # get satellite position in etrf2000
            position_etrs = convert_crs(7909, 7931, position_itrs_2k)
            
            # get satellite position in rd
            position_de = convert_crs(7931,  25832, position_etrs)
            
            tle[key].append(position_de)
    
    return satellites





    
'''
IMPORTANT TO KNOW:
    " Newer realizations of WGS84 are coincident with International Terrestrial
      Refernce Frame (ITRF), see Section 9.2, at about the 10­centimeter level. 
      For these realizations there are no official transformation parameters. 
      Newer realizations are adjusted occasionally in order to update the tracking 
      station coordinates for plate
      velocity ... In general WGS84 is identical to the ITRS and its realizations
      ITRFyy at the one meter level. Therefore, in practice, when precision does 
      not really matter and the user is satisfied with coordinates at the one meter 
      level, coordinates in ITRS, or derivatives of ITRS (like the European ETRS89) 
      are sometimes simply referred to as WGS84" - (van der Marel,2020,p61)
'''

'''
IMPORTANT TO KNOW
    "Web-­based Precise Point Positioning (PPP) services, which utilize satellite 
     orbits and clocks from the International GNSS Service (IGS), allow GPS users 
     to directly compute positions in ITRF. At the same time many regional and national 
     institutions have densified the IGS network to provide dense regional and nationaL
     networks of station coordinates in the ITRF." (van der Marel,2020,p65)
'''

'''
IMPORTANT TO KNOW
    "When working with the ITRF it is typical to provide coordinates as Cartesian 
     coordinates. However, the user is free to convert these into geographic coordinates. 
     The recommended ellipsoid for ITRS is the GRS80 ellipsoid, see Table 6.1. 
     This is the same ellipsoid as used for instance by WGS84." (van der Marel,2020,p65)
'''

'''
IMPORTANT TO KNOW
    "The ETRS89 system is realized in several ways, and like with ITRS, realizations
     of a system are called reference frames. By virtue of the ETRS89 definition, 
     which ties ETRS89 to ITRS at epoch 1989.0 and the Eurasian plate, for each realization
     of the ITRS (called ITRFyy), also a corresponding frame in ETRS89 can be computed.
     These frames are labelled ETRFyy. Each realization has a new set of improved positions
     and velocities. The three most recent realizations of ETRS89 are ETRF2000, ETRF2005
     and ETRF2014. Since each realization also reflects improvements in the datum 
     definition of ITRF, which results in small jumps in the coordinate time series,
     the EUREF Governing Board (formally Technical Working Group) recommends not to
     use the ETRF2005 for practical applications, and instead to adopt ETRF2000
     as a conventional frame of the ETRS89 system. However, considering the diverse needs
     of individual countries’, it is the countries’ decision to adopt their preferred
     ETRS89 realization. Most countries adopted the recommended ETRF2000, but not every
     European country has, and considering the improved accuracy and stability of the 
     ITRF2014, some could switch to ETRF2014." (van der Marel,2020,p66)
'''

'''
IMPORTANT TO KNOW
    "The Dutch network RTK services are certified by the Dutch Kadaster and thus 
     provide a national realizationof ETRS89 linked to the ETRF2000 reference frame."
     (van der Marel,2020,p67)
'''

'''
IMPORTANT TO KNOW
    "The increasing use of GPS resulted in a redefinition of RD in 2000, whereby
     from 2000 onwards RD was linked directly to ETRS89 through a transformation
     procedure called RDNAPTRANS." (van der Marel,2020,p69)
'''

'''
IMPORTANT TO KNOW
    "This projection is conformal, which means the projection is free from angular
     distortion, and that lines intersecting at any specified angle on the ellipsoid
     project into lines intersecting at the same angle on the projection."
     (van der Marel,p70,2020)
'''


 
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
