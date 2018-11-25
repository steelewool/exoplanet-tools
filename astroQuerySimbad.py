from astroquery.simbad import Simbad


# works on linux:
# table = Simbad.query_object ('M *', wildcard=True, verbose=False, get_query_payload=False)

# works on Raspberry PI
table = Simbad.query_object ('M *', wildcard=True, verbose=False)

print Simbad.list_votable_fields()

print 'Attempt a query_object command'

x=Simbad.query_object('M 1')
x.pprint()

print(x)

# table.pprint()

print table[0]
print table[49]

# Print the name of the columns

print table.colnames

#print the RA and Declination:

print 'RA and Declination'

print table[0]['RA']
print table[0]['DEC']

# for x in table:
#    print x


