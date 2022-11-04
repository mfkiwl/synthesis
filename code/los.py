"""
Created on Tue Oct  4 08:14:59 2022

@author: TM
"""
import os
from sgp4.api import Satrec
from sgp4.api import jday
from astropy.coordinates import TEME, CartesianDifferential, CartesianRepresentation
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import ITRS,EarthLocation
from pyproj import Transformer
from tle import tle_json


def tle_to_itrs(l1,l2,day=20,month=9,year=2022,hour=15,minute=0,second=0):
    
    satellite = Satrec.twoline2rv(l1, l2)
    jd, fr = jday(year,month,day,hour,minute,second)
    error, teme_position, teme_velocity = satellite.sgp4(jd,fr)
    
    # transform teme_position vector to True Equator Mean Equinox object
    teme_position = CartesianRepresentation(teme_position*u.km)
    teme_velocity = CartesianDifferential(teme_velocity*u.km/u.s)
    teme = TEME(teme_position.with_differentials(teme_velocity),obstime=Time(jd,format='jd'))
    #print(teme)
    
    # transform teme object to itrs geocentric coordinates
    itrs = teme.transform_to(ITRS(obstime=Time(jd,format='jd')))
    pos = [itrs.earth_location.geocentric[0].value*1000,itrs.earth_location.geocentric[1].value*1000,itrs.earth_location.geocentric[2].value*1000]
    #print(pos)
    #position_itrs = (itrs.earth_location.geodetic.lon.value,itrs.earth_location.geodetic.lat.value,itrs.earth_location.geodetic.height.value*1000)
    #print(itrs.earth_location.geodetic.lon)
    return pos#position_itrs

def convert_crs(from_crs,to_crs,pos_vector):
    
    transformer =  Transformer.from_crs(from_crs,to_crs,True)
    position = transformer.transform(pos_vector[0],pos_vector[1],pos_vector[2])
    return position





def german3d(pos_vector):
    
    transformer = Transformer.from_pipeline('''+proj=pipeline
                                               +step +proj=axisswap +order=2,1
                                               +step +proj=unitconvert +xy_in=deg +z_in=m +xy_out=rad +z_out=m
                                               +step +proj=utm +zone=32 +ellps=GRS80''')
    position = transformer.transform(pos_vector[0],pos_vector[1],pos_vector[2])
    #print(position)
    return position

def sat_pos():
    
    print(os.getcwd())
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
            print('ITRF2014:\t{}'.format(position_itrs_2k14))
            
            # get satellite position in itrf2000
            #position_itrs_2k = convert_crs(7789,7909,position_itrs_2k14)
            position_itrs_2k = convert_crs(7789,4919, position_itrs_2k14)
            print('ITRF2000:\t{}'.format(position_itrs_2k))
            
            # get satellite position in etrf2000
            position_etrs = convert_crs(4919,7930, position_itrs_2k)
            print('ETRF2000:\t{}'.format(position_etrs))

            # convert to 3d crs
            position_de = convert_crs(7930,32632,position_etrs)
            print('ETRS89/UTM32:\t{}'.format(position_de))
            
            tle[key].append(position_etrs)
    
    idx = 1
    
    point = convert_crs(32632,7930,[418521.00, 5653473.00, 346.42])
    #point = convert_crs(4979,4919,wgs84) 
    '''
    with open('../data/sat.txt','w') as sat:
        sat.write('id|wkt\n')
        for i in satellites[1]:
            v = [(i[list(i.keys())[0]][2][0] - p[0])*0.2, (i[list(i.keys())[0]][2][1] - p[1])*0.2, (i[list(i.keys())[0]][2][2] - p[2])]
            sat.write('{}|LINESTRING(418521.00 5653473.00 346.42,{} {} {})\n'.format(idx,v[0],v[1],v[2]))
            idx += 1
    '''
    
    
    with open('../data/sat.obj','w') as lines:
        for i in satellites[1]:
            v = [(i[list(i.keys())[0]][2][0]-point[0]), (i[list(i.keys())[0]][2][1]-point[1]), (i[list(i.keys())[0]][2][2]-point[2])]
            #print(v)
            lines.write('v {} {} {}\n'.format(round(v[0],2),round(v[1],2),round(v[2],2)))
        lines.write('v {} {} {}\n'.format(point[0],point[1],point[2]))
        for i in range(1,29):
            lines.write('l 29 {}\n'.format(i))
    
            
            
    return satellites

if __name__ == "__main__":
    sat_pos()


    
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
     provide a national realization of ETRS89 linked to the ETRF2000 reference frame."
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
