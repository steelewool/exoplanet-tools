
import commands

sourceFiles = commands.getoutput('find open_exoplanet_catalogue/ -name "*.xml"').splitlines()

count = 1
for file in sourceFiles:
    print "Line #%d: %s" % (count, file)
    count += 1
