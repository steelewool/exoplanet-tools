from astroquery import open_exoplanet_catalogue as oec
from astroquery.open_exoplanet_catalogue import findvalue

cata = oec.get_catalogue()

for planet in cata.findall(".//planet"):
    print (findvalue (planet, 'name'))
           
