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
from pyproj import Proj, Transformer
#from astropy.coordinates import EarthLocation, AltAz, Angle
#from astropy.units import Quantity
#from datetime import datetime

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
print(location_itrs.geodetic)


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
print(location_etrs)


'''
ETRF2000 TO RD
'''
# define the Rijksdriehoek projection system (EPSG 28992) 
epsg28992 = SpatialReference()
epsg28992.ImportFromEPSG(28992)

# define ogr transformer between 2 CRSs
etrf2rd = CoordinateTransformation(epsg7931, epsg28992)
location_rd = etrf2rd.TransformPoint(location_etrs[0],location_etrs[1],location_etrs[2])
print(location_rd)


'''
now, assume itrs realisation == itrf2000
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
print(location_etrs)

'''
ETRF2000 TO RD
'''

# define pyproj Transformer between 2 CRSs
etrf2rd_proj = Transformer.from_crs(7931, 28992)

# transform etrf2000 to rd
location_rd = etrf2rd_proj.transform(location_etrs[0],location_etrs[1],location_etrs[2])
print(location_rd)









 
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
