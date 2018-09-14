# The goal of this program is to scan all of the trasiting exoplanets
# and search for ones that are visible during the night times hours
# from my location.

# I was to eventually develop a GUI for this program. But for now the
# funcionality will live in this code.


import xml.etree.ElementTree as ET

import fnmatch
import os

# astropy.time is not working on my Raspberry Pi. This is issue # 3
# and more details can be found there.

import astropy.time
from astropy.time import Time
from astropy.time import TimeDelta

import time
from datetime import date
from datetime import datetime

from astropy.coordinates import EarthLocation,SkyCoord
from astropy import units as u
from astropy.coordinates import AltAz

import cmath

import commands

# This section of code creates a directory 'xml_files' and then creates
# softlinks to the xml files. in the open_exoplanet_catalogue in the
# The goal of this program is to scan all of the trasiting exoplanets
# and search for ones that are visible during the night times hours
# from my location.

# I was to eventually develop a GUI for this program. But for now the
# funcionality will live in this code.


import xml.etree.ElementTree as ET

import fnmatch
import os

# astropy.time is not working on my Raspberry Pi. This is issue # 3
# and more details can be found there.

import astropy.time
from astropy.time import Time
from astropy.time import TimeDelta

import time
from datetime import date
from datetime import datetime

from astropy.coordinates import EarthLocation,SkyCoord
from astropy import units as u
from astropy.coordinates import AltAz

import cmath

import commands

# This section of code creates a directory 'xml_files' and then creates
# softlinks to the xml files. in the open_exoplanet_catalogue in the
# systems and systems_kepler directories. This was done because I could not
# figure out how to scan though the list of xml files from two directories.
# The other side benefit was that I could remove xml files that were
# causing my program to crash.

commands.getstatusoutput ('ls')
commands.getstatusoutput ('rm -rf xml_files')
commands.getstatusoutput ('mkdir xml_files')
commands.getstatusoutput ('cd xml_files; ln -s ../open_exoplanet_catalogue/systems/* .;cd ..')
commands.getstatusoutput ('cd xml_files; ln -s ../open_exoplanet_catalgue//systems_kepler/* .;cd ..')

commands.getstatusoutput ('rm xml_files/WISE*.xml')
commands.getstatusoutput ('rm xml_files/PSO?J318.5-22.xml')
commands.getstatusoutput ('rm xml_files/CFBDSIR2149.xml')
commands.getstatusoutput ('rm xml_files/KOI-2700.xml')
commands.getstatusoutput ('rm xml_files/KIC?12557548.xml')
commands.getstatusoutput ('rm xml_files/EPIC?204129699.xml')
commands.getstatusoutput ('rm xml_files/EPIC?201637175.xml')
commands.getstatusoutput ('rm xml_files/SIMP0136+0933.xml')
commands.getstatusoutput ('rm xml_files/SDSS?J1110+0116.xml')

# This creates a list of all of the files in systems and systems_kepler.
# If I can get this working in the 'for file' I won't need the silly softlinks

# As of 2018-08-29 the variable 'fileList' is NOT used in the program.

fileList = (os.listdir('open_exoplanet_catalogue/systems') and os.listdir('open_exoplanet_catalogue/systems_kepler'))

#worked: for file in (os.listdir('systems') and os.listdir('systems_kepler')):

#worked: for file in os.listdir('systems'):
#worked:     if fnmatch.fnmatch(file,'*xml'):
#worked:         print "file: ", file
#worked:         tree = ET.parse('systems/'+file)

# Count is the number of objects that have been identifed. The variable is
# initialized to 0 here.

count = 0

# Set up by grabbing the current date and then using the Time object from astropy.time

dateTime    = datetime.today()
nowPST      = Time (dateTime, scale='utc')

dateTimeUTC = datetime.utcnow()
now         = Time (dateTimeUTC, scale='utc')

# For testing hardwire a date/time range

print 'a'
observingRange = ['2018-08-31T01:00:00','2018-09-01T23:00:00']
rangeTime = Time(observingRange, format='isot', scale='utc')

for file in os.listdir('xml_files'):

    #    print file
    
# Because of the way I set my the xml_files directory all of the files are xml files

    if fnmatch.fnmatch(file, '*.xml'):
        
        tree = ET.parse ('xml_files/'+file)
        root = tree.getroot();

        try: 
            star = tree.find('.//star')
        except:
            print'tree.find raised an exception'
            
        for planet in star.findall('.//planet'):
            if planet.findtext ('istransiting') == '1':

                if star.findtext('magV') != None:
                    mag = star.findtext('magV')
                else:
                    if star.findtext('magB') != None:
                        mag = star.findtext('magB')
                    else:
                        if star.findtext('magJ') != None:
                            mag = star.findtext('magJ')

                planetPeriod = planet.findtext('period')

                # Look for a valid looking period, one that is not '' nore 'None'.
            
                if planetPeriod != '' and planetPeriod != None:

                    planetPeriod = float(planetPeriod)

                    if planet.findtext('transittime') != None:

                        transitTimeBJD = float(planet.findtext('transittime'))

                        transitTime = Time(transitTimeBJD, format = 'jd', scale='utc')

#                        print 'transitTimeBJD     : ', transitTimeBJD
#                        print 'transitTime        : ', transitTime
#                        print 'now.jd             : ', now.jd
#                        print 'The now.jd is different that what I saw in computation on the web.'
                        
                        delta  = now.jd - transitTimeBJD;

#                        print 'delta              : ', delta
 
                        revolutionCount = delta / planetPeriod

                        intRevolutionCount = int(revolutionCount) + 1
                        nextTransit = transitTimeBJD + (intRevolutionCount * planetPeriod)

                        nextTransitTime = Time (nextTransit, format ='jd', scale = 'utc');

                        daysToTransit = nextTransit - now.jd

#                        print 'nextTransitTime    : ', nextTransitTime
#                        print 'daysToTransit      : ', daysToTransit

#
# Change the time to PST by subtracting 8 hours from the UTC time
#

                        nextTransitTimePST = nextTransit - (1.0/24.0*7.0)
                        nTTPST = Time (nextTransitTimePST, format='jd', scale='utc')

#                        print 'nextTransitTimePST      : ', nextTransitTimePST
#                        print 'nTTPST                  : ', nTTPST
#                        print 'nTTPST.jd               : ', nTTPST.jd
#                        print 'nTTPST.fits             : ', nTTPST.fits
                        
                        starRadius   = star.findtext('radius')
                        if (starRadius == None):
                            starRadius = float(0.0)
                        else:
                            starRadius    = float(starRadius) * 1.3914 * 1000000

                        planetRadius   = planet.findtext('radius')
                        
                        if (planetRadius == None):
                            planetRadius = 0.0
                        else:
                            planetRadius = float(planetRadius) * 139822

                        if (starRadius != 0) and (planetRadius != 0):
                            starArea = cmath.pi * starRadius * starRadius
                            planetArea = cmath.pi * planetRadius * planetRadius
                            planetStarAreaRatio = planetArea / starArea
                        else:
                            planetStarAreaRatio = 0
                            
                        a = nextTransitTimePST
                        b = nowPST.jd + 1
                        c = a < b

# d start off as false and is det to true if the time is in the specifed time range

                        d = False
                        if nTTPST > rangeTime[0]:
                            if nTTPST < rangeTime[1]:
                                d = True

# e = sideral_time('apparent',longitude=None,model=None)

                        observingPosition = EarthLocation(lat=34*u.deg, lon=-118*u.deg, height=500*u.m)  

                        observingNextTransitTime = Time(nextTransitTime.fits)


# Eliminate objects based on
# a) magnitude of the star must be great the 11th magnitude,
# b) planetStarRation at least 0.01,
# c) variable 'd' (poorly named) is not true, that is the object will be eliminated if the transit is
#    not within the specified time range
# d) altitude of the object is not at least 10 degrees above the horizon
# e) transit happens during daylight hours, between 4 & 18 time.
# Still need to only output if the transit happens at night.

                        aa = AltAz(location=observingPosition, obstime=observingNextTransitTime)

                        ra = root.findtext('rightascension')
                        dec = root.findtext('declination')
                            
                        raHrMinSec = ra[0:2] + 'h' + ra[3:5] + 'm' + ra[6:8] + 's'
                        decDegMinSec = dec[0:3] + 'd' + dec[4:6] + 'm' + dec[8:10] + 's'
                            
                        skyCoord = SkyCoord (raHrMinSec + ' ' + decDegMinSec, frame='icrs')

                        altAzi = skyCoord.transform_to(AltAz(obstime=observingNextTransitTime,location=observingPosition))

# Looking for hour of transit (PST). For now day time is between 06 and 17 hours. Night would be defined
# as true if we are not in this range:

                        hour = nTTPST.fits[11:13];
                        if (hour > '04' and hour < '17'):
                            night = False
                        else:
                            night = True

                        if (float(mag) < 11) and d and (planetStarAreaRatio >= 0.01) and (altAzi.alt.degree > 20) and (night):
                            count = count + 1

                            print '------------------'
                            print 'file name                : ', file
                            print 'System name              : ', root.findtext('name')
                            print 'Planet name              : ', planet.findtext('name')
                            print 'Planet period            : ', planet.findtext('period')
                            print 'System Right Ascension   :  ', root.findtext('rightascension')
                            print 'System Declination       : ', root.findtext('declination')
                            print 'System Magnitude         : ', mag
                            print 'observingNextTransitTime : ', observingNextTransitTime
                            print 'Azimuth                  : ', altAzi.az.degree
                            print 'Altitude                 : ', altAzi.alt.degree
                            print 'Days until transit       : ', daysToTransit
                            print 'nTTPST.jd                : ', nTTPST.jd
                            print 'nTTPST.fits              : ', nTTPST.fits, 'PST'


#                            print 'transitTime.jd           : ', transitTime.jd
#                            print 'transitTime.fits         : ', transitTime.fits
#                            print 'nextTransit              : ', nextTransit
#                            print 'dateTimeUTC              : ', dateTimeUTC
#                            print 'now jd                   : ', now.jd
#                            print 'now fits                 : ', now.fits
#                            print 'delta                    : ', delta
#                            print 'revolutionCount          : ', revolutionCount
#                            print 'int revoultionCount      : ', int(revolutionCount) + 1
#                            print 'raHrMinSec               : ', raHrMinSec
#                            print 'decDegMinSec             : ', decDegMinSec
#                            print 'transitTimeBJD           : ', transitTimeBJD
#                            print 'now                      : ', now
#                            print 'nextTransitTime          : ', nextTransitTime.fits
#                            print 'nextTransitTimePST       : ', nextTransitTimePST
#                            print 'Star radius              : ', starRadius
#                            print 'Planet radius            : ', planetRadius

                            print 'Planet/Star area ratio   : ', planetStarAreaRatio
                            
                            print 'count                    : ', count
                            




        






# systems and systems_kepler directories. This was done because I could not
# figure out how to scan though the list of xml files from two directories.
# The other side benefit was that I could remove xml files that were
# causing my program to crash.

commands.getstatusoutput ('ls')
commands.getstatusoutput ('rm -rf xml_files')
commands.getstatusoutput ('mkdir xml_files')
commands.getstatusoutput ('cd xml_files; ln -s ../open_exoplanet_catalogue/systems/* .;cd ..')
commands.getstatusoutput ('cd xml_files; ln -s ../open_exoplanet_catalgue//systems_kepler/* .;cd ..')

commands.getstatusoutput ('rm xml_files/WISE*.xml')
commands.getstatusoutput ('rm xml_files/PSO?J318.5-22.xml')
commands.getstatusoutput ('rm xml_files/CFBDSIR2149.xml')
commands.getstatusoutput ('rm xml_files/KOI-2700.xml')
commands.getstatusoutput ('rm xml_files/KIC?12557548.xml')
commands.getstatusoutput ('rm xml_files/EPIC?204129699.xml')
commands.getstatusoutput ('rm xml_files/EPIC?201637175.xml')
commands.getstatusoutput ('rm xml_files/SIMP0136+0933.xml')
commands.getstatusoutput ('rm xml_files/SDSS?J1110+0116.xml')

# This creates a list of all of the files in systems and systems_kepler.
# If I can get this working in the 'for file' I won't need the silly softlinks

# As of 2018-08-29 the variable 'fileList' is NOT used in the program.

fileList = (os.listdir('open_exoplanet_catalogue/systems') and os.listdir('open_exoplanet_catalogue/systems_kepler'))

#worked: for file in (os.listdir('systems') and os.listdir('systems_kepler')):

#worked: for file in os.listdir('systems'):
#worked:     if fnmatch.fnmatch(file,'*xml'):
#worked:         print "file: ", file
#worked:         tree = ET.parse('systems/'+file)

# Count is the number of objects that have been identifed. The variable is
# initialized to 0 here.

count = 0

# Set up by grabbing the current date and then using the Time object from astropy.time

dateTime    = datetime.today()
nowPST      = Time (dateTime, scale='utc')

dateTimeUTC = datetime.utcnow()
now         = Time (dateTimeUTC, scale='utc')

# For testing hardwire a date/time range. As of right now the time
# is in the UTC time zone. The inputs SHOULD really be in the
# local time zone. For me this would either be PST or PDT.

print 'b'
observingRange = ['2018-08-31T01:00:00','2018-09-01T23:00:00']
rangeTime = Time(observingRange, format='isot', scale='utc')

for file in os.listdir('xml_files'):

    #    print file
    
# Because of the way I set my the xml_files directory all of the files are xml files

    if fnmatch.fnmatch(file, '*.xml'):
        
        tree = ET.parse ('xml_files/'+file)
        root = tree.getroot();

        try: 
            star = tree.find('.//star')
        except:
            print'tree.find raised an exception'
            
        for planet in star.findall('.//planet'):
            if planet.findtext ('istransiting') == '1':

                if star.findtext('magV') != None:
                    mag = star.findtext('magV')
                else:
                    if star.findtext('magB') != None:
                        mag = star.findtext('magB')
                    else:
                        if star.findtext('magJ') != None:
                            mag = star.findtext('magJ')

                planetPeriod = planet.findtext('period')

                # Look for a valid looking period, one that is not '' nore 'None'.
            
                if planetPeriod != '' and planetPeriod != None:

                    planetPeriod = float(planetPeriod)

                    if planet.findtext('transittime') != None:

                        transitTimeBJD = float(planet.findtext('transittime'))

                        transitTime = Time(transitTimeBJD, format = 'jd', scale='utc')

#                        print 'transitTimeBJD     : ', transitTimeBJD
#                        print 'transitTime        : ', transitTime
#                        print 'now.jd             : ', now.jd
#                        print 'The now.jd is different that what I saw in computation on the web.'
                        
                        delta  = now.jd - transitTimeBJD;

#                        print 'delta              : ', delta
 
                        revolutionCount = delta / planetPeriod

                        intRevolutionCount = int(revolutionCount) + 1
                        nextTransit = transitTimeBJD + (intRevolutionCount * planetPeriod)

                        nextTransitTime = Time (nextTransit, format ='jd', scale = 'utc');

                        daysToTransit = nextTransit - now.jd

#                        print 'nextTransitTime    : ', nextTransitTime
#                        print 'daysToTransit      : ', daysToTransit

#
# Change the time to PST by subtracting 8 hours from the UTC time
#

                        nextTransitTimePST = nextTransit - (1.0/24.0*7.0)
                        nTTPST = Time (nextTransitTimePST, format='jd', scale='utc')

#                        print 'nextTransitTimePST      : ', nextTransitTimePST
#                        print 'nTTPST                  : ', nTTPST
#                        print 'nTTPST.jd               : ', nTTPST.jd
#                        print 'nTTPST.fits             : ', nTTPST.fits
                        
                        starRadius   = star.findtext('radius')
                        if (starRadius == None):
                            starRadius = float(0.0)
                        else:
                            starRadius    = float(starRadius) * 1.3914 * 1000000

                        planetRadius   = planet.findtext('radius')
                        
                        if (planetRadius == None):
                            planetRadius = 0.0
                        else:
                            planetRadius = float(planetRadius) * 139822

                        if (starRadius != 0) and (planetRadius != 0):
                            starArea = cmath.pi * starRadius * starRadius
                            planetArea = cmath.pi * planetRadius * planetRadius
                            planetStarAreaRatio = planetArea / starArea
                        else:
                            planetStarAreaRatio = 0
                            
                        a = nextTransitTimePST
                        b = nowPST.jd + 1
                        c = a < b

# d start off as false and is det to true if the time is in the specifed time range

                        d = False
                        if nTTPST > rangeTime[0]:
                            if nTTPST < rangeTime[1]:
                                d = True

# e = sideral_time('apparent',longitude=None,model=None)

                        observingPosition = EarthLocation(lat=34*u.deg, lon=-118*u.deg, height=500*u.m)  

                        observingNextTransitTime = Time(nextTransitTime.fits)


# Eliminate objects based on
# a) magnitude of the star must be great the 11th magnitude,
# b) planetStarRation at least 0.01,
# c) variable 'd' (poorly named) is not true, that is the object will be eliminated if the transit is
#    not within the specified time range
# d) altitude of the object is not at least 10 degrees above the horizon
# e) transit happens during daylight hours, between 4 & 18 time.
# Still need to only output if the transit happens at night.

                        aa = AltAz(location=observingPosition, obstime=observingNextTransitTime)

                        ra = root.findtext('rightascension')
                        dec = root.findtext('declination')
                            
                        raHrMinSec = ra[0:2] + 'h' + ra[3:5] + 'm' + ra[6:8] + 's'
                        decDegMinSec = dec[0:3] + 'd' + dec[4:6] + 'm' + dec[8:10] + 's'
                            
                        skyCoord = SkyCoord (raHrMinSec + ' ' + decDegMinSec, frame='icrs')

                        altAzi = skyCoord.transform_to(AltAz(obstime=observingNextTransitTime,location=observingPosition))

# Looking for hour of transit (PST). For now day time is between 06 and 17 hours. Night would be defined
# as true if we are not in this range:

                        hour = nTTPST.fits[11:13];
                        if (hour > '04' and hour < '17'):
                            night = False
                        else:
                            night = True

                        if (float(mag) < 11) and d and (planetStarAreaRatio >= 0.01) and (altAzi.alt.degree > 20) and (night):
                            count = count + 1

                            print '------------------'
                            print 'file name                : ', file
                            print 'System name              : ', root.findtext('name')
                            print 'Planet name              : ', planet.findtext('name')
                            print 'Planet period            : ', planet.findtext('period')
                            print 'System Right Ascension   :  ', root.findtext('rightascension')
                            print 'System Declination       : ', root.findtext('declination')
                            print 'System Magnitude         : ', mag
                            print 'observingNextTransitTime : ', observingNextTransitTime
                            print 'Azimuth                  : ', altAzi.az.degree
                            print 'Altitude                 : ', altAzi.alt.degree
                            print 'Days until transit       : ', daysToTransit
                            print 'nTTPST.jd                : ', nTTPST.jd
                            print 'nTTPST.fits              : ', nTTPST.fits, 'PST'


#                            print 'transitTime.jd           : ', transitTime.jd
#                            print 'transitTime.fits         : ', transitTime.fits
#                            print 'nextTransit              : ', nextTransit
#                            print 'dateTimeUTC              : ', dateTimeUTC
#                            print 'now jd                   : ', now.jd
#                            print 'now fits                 : ', now.fits
#                            print 'delta                    : ', delta
#                            print 'revolutionCount          : ', revolutionCount
#                            print 'int revoultionCount      : ', int(revolutionCount) + 1
#                            print 'raHrMinSec               : ', raHrMinSec
#                            print 'decDegMinSec             : ', decDegMinSec
#                            print 'transitTimeBJD           : ', transitTimeBJD
#                            print 'now                      : ', now
#                            print 'nextTransitTime          : ', nextTransitTime.fits
#                            print 'nextTransitTimePST       : ', nextTransitTimePST
#                            print 'Star radius              : ', starRadius
#                            print 'Planet radius            : ', planetRadius

                            print 'Planet/Star area ratio   : ', planetStarAreaRatio
                            
                            print 'count                    : ', count
                            




        






