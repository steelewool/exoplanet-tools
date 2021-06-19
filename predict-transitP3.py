# The goal of this program is to scan all of the trasiting exoplanets
# and search for ones that are visible during the night times hours
# from my location.

# I was to eventually develop a GUI for this program. But for now the
# funcionality will live in this code.

from astropy import units as u

import astropy.time
from astropy.time import Time
from astropy.time import TimeDelta

from astropy.coordinates import EarthLocation,SkyCoord
from astropy.coordinates import AltAz

import cmath
import subprocess

from datetime import timedelta
from datetime import date
from datetime import datetime

import fnmatch

import os

import time

import xml.etree.ElementTree as ET

# This section of code creates a directory 'xml_files' and then creates
# softlinks to the xml files. in the open_exoplanet_catalogue in the
# The goal of this program is to scan all of the trasiting exoplanets
# and search for ones that are visible during the night times hours
# from my location.

# I was to eventually develop a GUI for this program. But for now the
# funcionality will live in this code.

# astropy.time is not working on my Raspberry Pi.
# This is issue # 3 and more details can be found there.

# This section of code creates a directory 'xml_files' and then creates
# softlinks to the xml files. in the open_exoplanet_catalogue in the
# systems and systems_kepler directories. This was done because I could not
# figure out how to scan though the list of xml files from two directories.
# The other side benefit was that I could remove xml files that were
# causing my program to crash.

subprocess.getstatusoutput ('ls')
subprocess.getstatusoutput ('rm -rf xml_files')
subprocess.getstatusoutput ('mkdir xml_files')

# subprocess.getstatusoutput ('cd xml_files; ln -s ../OpenExopolanetCatalogue/open_exoplanet_catalogue/systems/* .;cd ..')
# subprocess.getstatusoutput ('cd xml_files; ln -s ../OpenExopolanetCatalogue/open_exoplanet_catalogue/systems_kepler/* .;cd ..')

subprocess.getstatusoutput ('cd xml_files; cp ../open_exoplanet_catalogue/systems/*        .;cd ..')
subprocess.getstatusoutput ('cd xml_files; cp ../open_exoplanet_catalogue/systems_kepler/* .;cd ..')

subprocess.getstatusoutput ('rm xml_files/EPIC?211901114.xml')
subprocess.getstatusoutput ('rm xml_files/EPIC?201637175.xml')
subprocess.getstatusoutput ('rm xml_files/KIC?12557548.xml')
subprocess.getstatusoutput ('rm xml_files/SDSS?J1110+0116.xml')
subprocess.getstatusoutput ('rm xml_files/PSO?J318?5-22.xml')
subprocess.getstatusoutput ('rm xml_files/SIMP0136+0933.xml')
subprocess.getstatusoutput ('rm xml_files/CFBDSIR2149.xml')
subprocess.getstatusoutput ('rm xml_files/WISE?0855-0714.xml')
subprocess.getstatusoutput ('rm xml_files/EPIC?204129699.xml')
subprocess.getstatusoutput ('rm xml_files/EPIC?201637175.xml')

# This creates a list of all of the files in systems and systems_kepler.
# If I can get this working in the 'for file' I won't need the silly
# softlinks

# As of 2018-08-29 the variable 'fileList' is NOT used in the program.

fileList = (os.listdir('open_exoplanet_catalogue/systems') and
            os.listdir('open_exoplanet_catalogue/systems_kepler'))

# Count is the number of objects that have been identifed. The variable is
# initialized to 0 here.

count = 0

# Set up by grabbing the current date and then using the Time object
# from astropy.time

dateTime    = datetime.today()
nowPST      = Time (dateTime, scale='utc')

dateTimeUTC = datetime.utcnow()
now         = Time (dateTimeUTC, scale='utc')

# This will search for 2 days (timedelta(2))

startTime = Time(datetime.now(),              scale='utc')
endTime   = Time(datetime.now()+timedelta(2), scale='utc')

print ('startTime: ', startTime)
print ('endTime  : ', endTime)

observingMorningTime = '04'
observingEveningTime = '17'

minMagCutoff           = input ('Enter minimum magnitude   : ')
minAltCutoff           = input ('Enter minimum altitude    : ')
minPlanetStarAreaRatio = input ('Enter minimum area ration : ')

# This reads into 'file' all of the files in the xml_files directory

for file in os.listdir('xml_files'):
    
#    print ("Debugging, file: ", file)
    
# Because of the way I set my the xml_files directory all of the files
# are xml files

# This if statement may not in fact be necessary. Need to confirm this.

    if fnmatch.fnmatch(file, '*.xml'):
        
        tree = ET.parse ('xml_files/'+file)
        root = tree.getroot();

# Look for a star field in the xml - if there isn't raise an exception. This
# exception is not having in the current set of xml files.

        try: 
            star = tree.find('.//star')
        except:
            print ('tree.find raised an exception')

# Look through all of the possible planets in a system.

        for planet in star.findall('.//planet'):
            if planet.findtext ('istransiting') == '1':

# Get the magntiude of the star. Use the visual magnitude if it is available.
# If not, use the 'B' magnitude and if that isn't available use the
# 'J' magnitude.

                if star.findtext('magV') != None:
                    mag = star.findtext('magV')
                else:
                    if star.findtext('magB') != None:
                        mag = star.findtext('magB')
                    else:
                        if star.findtext('magJ') != None:
                            mag = star.findtext('magJ')
                        else:
                            mag = 0.0
                            
                planetPeriod = planet.findtext('period')

                # Look for a valid looking period, one that is not ''
                # nore 'None'.
            
                if planetPeriod != '' and planetPeriod != None:

                    planetPeriod = float(planetPeriod)

                    if planet.findtext('transittime') != None:

# These two times, transitTimeBJD and transitTime are identical times.
# Need to pick out just one for the code.

                        transitTimeBJD = float(
                                           planet.findtext('transittime'))
                        transitTime = Time(transitTimeBJD,
                                           format = 'jd',
                                           scale='utc')

# 'now' is the current time. Not sure why I'm using the current time and not
# the range of time specified in the time range.
# It seems like this should be the start of the time range.

                        delta  = now.jd - transitTimeBJD;

                        revolutionCount = delta / planetPeriod

                        intRevolutionCount = int(revolutionCount) + 1

#                        print 'delta              : ', delta
#                        print 'revolutionCount    : ', revolutionCount
#                        print 'intRevolutionCount : ', intRevolutionCount
                        
                        nextTransit = transitTimeBJD + \
                                      (intRevolutionCount * planetPeriod)

                        nextTransitTime = Time (nextTransit,
                                                format ='jd',
                                                scale = 'utc');

                        daysToTransit = nextTransit - now.jd

#
# Change the time to PST by subtracting 8 hours from the UTC time
#

                        nextTransitTimePST = nextTransit - (1.0/24.0*8.0)
                        nTTPST = Time (nextTransitTimePST,
                                       format='jd',
                                       scale='utc')

                        starRadius   = star.findtext('radius')
                        if (starRadius == None):
                            starRadius = float(0.0)
                        else:
                            starRadius    = float(starRadius) * \
                                                  1.3914      * \
                                                  1000000

                        planetRadius   = planet.findtext('radius')
                        
                        if (planetRadius == None):
                            planetRadius = 0.0
                        else:
                            planetRadius = float(planetRadius) * 139822

                        if (starRadius != 0) and (planetRadius != 0):
                            starArea            = cmath.pi   * \
                                                  starRadius * \
                                                  starRadius
                            planetArea          = cmath.pi     * \
                                                  planetRadius * \
                                                  planetRadius
                            planetStarAreaRatio = planetArea / starArea
                        else:
                            planetStarAreaRatio = 0
                            
                        a = nextTransitTimePST
                        b = nowPST.jd + 1
                        c = a < b

# d start off as false and is det to true if the time is in the specifed
# time range

                        d = False
                        
                        if nTTPST.jd > startTime.jd:
                            if nTTPST.jd < endTime.jd:
                                d = True

# e = sideral_time('apparent',longitude=None,model=None)

                        observingPosition = EarthLocation(lat=34*u.deg,
                                                          lon=-118*u.deg,
                                                          height=500*u.m)  

                        observingNextTransitTime = Time(nextTransitTime.fits)


# Eliminate objects based on
# a) magnitude of the star must be great the 11th magnitude,
# b) planetStarRation at least 0.01,
# c) variable 'd' (poorly named) is not true, that is the object will be
#    eliminated if the transit is
#    not within the specified time range
# d) altitude of the object is not at least 10 degrees above the horizon
# e) transit happens during daylight hours, between 4 & 18 time.
# Still need to only output if the transit happens at night.

                        aa = AltAz(location=observingPosition,
                                   obstime=observingNextTransitTime)

                        ra = root.findtext('rightascension')
                        dec = root.findtext('declination')
                            
                        raHrMinSec   = ra[0:2]   + 'h' + ra[3:5]  + 'm' + \
                                       ra[6:8]   + 's'
                        decDegMinSec = dec[0:3]  + 'd' + dec[4:6] + 'm' + \
                                       dec[8:10] + 's'
                            
                        skyCoord = SkyCoord (raHrMinSec + ' ' + \
                                             decDegMinSec,
                                             frame='icrs')

                        altAzi = skyCoord.transform_to(
                                      AltAz(obstime=observingNextTransitTime,
                                            location=observingPosition))

# Looking for hour of transit (PST). For now day time is between
# 06 and 17 hours. Night would be defined as true if we are not in
# this range:

                        hour = nTTPST.fits[11:13];

                        if (hour > observingMorningTime and
                            hour < observingEveningTime):
                            night = False
                        else:
                            night = True
                        
                        if (float(mag) < float(minMagCutoff))                     and \
                           d                                               and \
                           (planetStarAreaRatio >= float(minPlanetStarAreaRatio)) and \
                           (altAzi.alt.degree > float(minAltCutoff))              and \
                           night:
                            count = count + 1

                            print ('------------------')
                            print ('file name                : ', file)
                            print ('System name              : ',  \
                                   root.findtext('name'))
                            print ('Planet name              : ',  \
                                   planet.findtext('name'))
                            print ('Planet period            : ',  \
                                   planet.findtext('period'))
                            print ('System Right Ascension   :  ', \
                                   root.findtext('rightascension'))
                            print ('System Declination       : ',  \
                                   root.findtext('declination'))
                            print ('System Magnitude         : ',  \
                                   mag)
                            print ('observingNextTransitTime : ',  \
                                   observingNextTransitTime)
                            print ('Azimuth                  : ',  \
                                   altAzi.az.degree)
                            print ('Altitude                 : ',  \
                                   altAzi.alt.degree)
                            print ('Days until transit       : ',  \
                                   daysToTransit)
                            print ('nTTPST.jd                : ',  \
                                   nTTPST.jd)
                            print ('nTTPST.fits              : ',  \
                                   nTTPST.fits, 'PST')
                            print ('Planet/Star area ratio   : ',  \
                                   planetStarAreaRatio)
                            print ('count                    : ',  \
                                   count)
