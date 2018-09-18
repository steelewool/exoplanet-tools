from astroquery.simbad import Simbad

table = Simbad.query_object ('M *', wildcard=True, verbose=False, get_query_payload=False)

table.pprint()

print table[49]

for x in table:
    print x


