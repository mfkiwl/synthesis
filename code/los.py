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
from osgeo.osr import SpatialReference, CoordinateTransformation
from pyproj import Transformer
#from astropy.coordinates import EarthLocation, AltAz, Angle
#from astropy.units import Quantity
#from datetime import datetime

print('Compute Satellite Position in TEME and Transform to ITRS')

# TLE for satellite GPS BIIRM-2 (PRN 31) 
s = '1 29486U 06042A   22277.55622356 -.00000022  00000+0  00000+0 0  9998'
t = '2 29486  54.7051 204.9391 0104917  23.9319 143.2008  2.00572655117289'
satellite = Satrec.twoline2rv(s,t)
jd, fr = jday(2022,9,26,13,0,0)
error, teme_position, teme_velocity = satellite.sgp4(jd,fr)

# transform teme_position vector to True Equator Mean Equinox object
teme_position = CartesianRepresentation(teme_position*u.km)
teme_velocity = CartesianDifferential(teme_velocity*u.km/u.s)
teme = TEME(teme_position.with_differentials(teme_velocity),obstime=Time(jd,format='jd'))

# transform teme position to itrs geocentric coordinates
itrs = teme.transform_to(ITRS(obstime=Time(jd,format='jd')))
location_itrs = itrs.earth_location
print('\tITRFyy From TEME\n\tx: {}\n\ty: {}\n\tz: {}'.format(location_itrs.geodetic.lon.value,location_itrs.geodetic.lat.value,location_itrs.geodetic.height.value))

print('\nASSUMING ITRF2000:')
print('\n1. TRANSFORM USING osgeo.osr')

'''
now, assume itrs realisation == itrf2000
transform itrf2000 to etrs2000
transform etrs2000 to bessel + project rd amersfoort
using SpatialReference and CoordinateTransformation
'''

'''
ITRF2000 TO ETRF2000
'''
# define ITRF2000 system (EPSG 7909)
epsg7909 = SpatialReference()
epsg7909.ImportFromEPSG(7909)


# define itrf2014 system (3D EPSG 7912)
# itrf2014 used based on trail and error; not very engineery but no precise info on links between these different systems
# reference made to Jochem Lesparre presentation 'STOP USING WGS84'
#epsg7912 = SpatialReference()
#epsg7912.ImportFromEPSG(7912)

# define ETRF2000 system (EPSG 7931)
epsg7931 = SpatialReference()
epsg7931.ImportFromEPSG(7931)

# define ogr transformer between 2 CRSs
itrf2etrf = CoordinateTransformation(epsg7909, epsg7931)
location_etrs = itrf2etrf.TransformPoint(location_itrs.geodetic.lon.value,location_itrs.geodetic.lat.value,location_itrs.geodetic.height.value)
print('\tETRF2000 From ITRFyy\n\tx: {}\n\ty: {}\n\tz: {}'.format(location_etrs[0],location_etrs[1],location_etrs[2]))

'''
ETRF2000 TO RD
'''
# define the Rijksdriehoek projection system (EPSG 28992) 
epsg28992 = SpatialReference()
epsg28992.ImportFromEPSG(28992)

# define ogr transformer between 2 CRSs
etrf2rd = CoordinateTransformation(epsg7931, epsg28992)
location_rd = etrf2rd.TransformPoint(location_etrs[0],location_etrs[1],location_etrs[2])
print('\tRD From ETRF2000\n\tx: {}\n\ty: {}\n\tz: {}'.format(location_rd[0],location_rd[1],location_rd[2]))

print('\n2. TRANSFORM USING pyproj')

'''
now,
assume itrs realisation == itrf2000
transform itrf2000 to etrs2000
transform etrs2000 to bessel + project rd amersfoort
using SpatialReference and CoordinateTransformation
'''

'''
ITRF2000 TO ETRF2000
'''
# define pyproj Transformer between 2 CRSs
itrf2etrf_proj = Transformer.from_crs(7909, 7931)

# transform itrf2000 to etrf2000
location_etrs = itrf2etrf_proj.transform(location_itrs.geodetic.lon.value,location_itrs.geodetic.lat.value,location_itrs.geodetic.height.value)
print('\tETRF2000 From ITRFyy\n\tx: {}\n\ty: {}\n\tz: {}'.format(location_etrs[0],location_etrs[1],location_etrs[2]))


'''
ETRF2000 TO RD
'''
# define pyproj Transformer between 2 CRSs
etrf2rd_proj = Transformer.from_crs(7931, 28992)

# transform etrf2000 to rd
location_rd = etrf2rd_proj.transform(location_etrs[0],location_etrs[1],location_etrs[2])
print('\tRD From ETRF2000\n\tx: {}\n\ty: {}\n\tz: {}'.format(location_rd[0],location_rd[1],location_rd[2]))



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
      are sometimes simply referred to as WGS84" - (van der Marel,p61,2020)
'''

'''
IMPORTANT TO KNOW
    "Web-­based Precise Point Positioning (PPP) services, which utilize satellite 
     orbits and clocks from the International GNSS Service (IGS), allow GPS users 
     to directly compute positions in ITRF. At the same time many regional and national 
     institutions have densified the IGS network to provide dense regional and nationaL
     networks of station coordinates in the ITRF." (van der Marel,p65,2020)
'''

'''
IMPORTANT TO KNOW
    "When working with the ITRF it is typical to provide coordinates as Cartesian 
     coordinates. However, the user is free to convert these into geographic coordinates. 
     The recommended ellipsoid for ITRS is the GRS80 ellipsoid, see Table 6.1. 
     This is the same ellipsoid as used for instance by WGS84." (van der Marel,p65,2020)
'''

'''
IMPORTANT TO KNOW
    "The ETRS89 system is realized in several ways, and like with ITRS, realizations
     of a sysTem are called reference frames. By virtue of the ETRS89 definition, 
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
     ITRF2014, some could switch to ETRF2014." (van der Marel,p66,2020)
'''

'''
IMPORTANT TO KNOW
    "The Dutch network RTK services are certified by the Dutch Kadaster and thus 
     provide a national realizationof ETRS89 linked to the ETRF2000 reference frame."
     (van der Marel,p67,2020)
'''

'''
IMPORTANT TO KNOW
    "The increasing use of GPS resulted in a redefinition of RD in 2000, whereby
     from 2000 onwards RD was linked directly to ETRS89 through a transformation
     procedure called RDNAPTRANS." (van der Marel,p69,2020)
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
