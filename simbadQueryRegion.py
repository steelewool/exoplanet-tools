from astroquery.simbad import Simbad
import astropy.coordinates as coord
import astropy.units as u

result_table = Simbad.query_region(coord.SkyCoord(ra=11.70, dec=10.90,unit=(u.deg, u.deg), frame='fk5'),radius=0.5 * u.deg,epoch='B1950',equinox=1950)
print(result_table)


result_table = Simbad.query_region(coord.SkyCoord("05h35m17.3s -05h23m28s", frame='icrs'),radius=0.5 * u.deg)

result_table = Simbad.query_region('m81',radius=0.1 * u.deg,epoch='B2000',equinox=2000)
