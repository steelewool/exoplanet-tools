
import xml.etree.ElementTree as ET

import fnmatch
import os

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

# import subprocess

# subprocess.call ('ls')
# subprocess.call ('rm xml_files/*')
# subprocess.call ('ln -s systems/*.xml xml_files/.')

import commands
commands.getstatusoutput ('ls')
commands.getstatusoutput ('rm -rf xml_files')
commands.getstatusoutput('mkdir xml_files')
commands.getstatusoutput ('cd xml_files; ln -s ../systems/*.xml .;cd ..')
commands.getstatusoutput ('cd xml_files; ln -s ../systems_kepler/*.xml .;cd ..')
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

fileList = (os.listdir('systems') and os.listdir('systems_kepler'))

#worked: for file in (os.listdir('systems') and os.listdir('systems_kepler')):

#worked: for file in os.listdir('systems'):
#worked:     if fnmatch.fnmatch(file,'*xml'):
#worked:         print "file: ", file
#worked:         tree = ET.parse('systems/'+file)

count = 0

# Set up by grabbing the current date and then using the Time object from astropy.time

dateTime    = datetime.today()
nowPST      = Time (dateTime, scale='utc')

dateTimeUTC = datetime.utcnow()
now         = Time (dateTimeUTC, scale='utc')

# For testing hardcore a date/time range

observingRange = ['2018-08-24T18:00:00','2018-08-26T23:00:00']
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
#                    print 'star.magV           : ', star.findtext('magV')
                    mag = star.findtext('magV')
                else:
                    if star.findtext('magB') != None:
#                        print 'star.magB      : ', star.findtext('magB')
                        mag = star.findtext('magB')
                    else:
                        if star.findtext('magJ') != None:
#                            print 'star.magJ             : ', star.findtext('magJ')
                            mag = star.findtext('magJ')

                planetPeriod = planet.findtext('period')

                # Look for a valid looking period, one that is not '' nore 'None'.
            
                if planetPeriod != '' and planetPeriod != None:

                    planetPeriod = float(planetPeriod)

                    if planet.findtext('transittime') != None:

                        transitTimeBJD = float(planet.findtext('transittime'))

                        transitTime = Time(transitTimeBJD, format = 'jd', scale='utc')

                        delta  = now.jd - transitTimeBJD;
                
                        revolutionCount = delta / planetPeriod

                        intRevolutionCount = int(revolutionCount) + 1
                        nextTransit = transitTimeBJD + (intRevolutionCount * planetPeriod)

                        nextTransitTime = Time (nextTransit, format ='jd', scale = 'utc');

                        daysToTransit = nextTransit - now.jd

#
# Change the time to PST by subtracting 8 hours from the UTC time
#

                        nextTransitTimePST = nextTransit - (1.0/24.0*8.0)
                        nTTPST = Time (nextTransitTimePST, format='jd', scale='utc')

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

                        d = False
                        if nTTPST > rangeTime[0]:
                            if nTTPST < rangeTime[1]:
                                d = True

# e = sideral_time('apparent',longitude=None,model=None)

                        observingPosition = EarthLocation(lat=34*u.deg, lon=-118*u.deg, height=500*u.m)  

                        observingNextTransitTime = Time(nextTransitTime.fits)
                        
                        if (float(mag) < 11) and d and (planetStarAreaRatio >= 0.01):
                            count = count + 1

                            print '------------------'
                            print 'observingPosition     : ', observingPosition
                            print 'observingNextTransitTime: ', observingNextTransitTime
                            
                            aa = AltAz(location=observingPosition, obstime=observingNextTransitTime)
                            print 'aa    : ', aa

                            print 'ra dec:', root.findtext('rightascension')+' '+root.findtext('declination')
                            # skyCoord = SkyCoord ('05h04m20s -06d13m47s', frame='icrs')

                            ra = root.findtext('rightascension')
                            dec = root.findtext('declination')
                            
                            raHrMinSec = ra[0:2] + 'h' + ra[3:5] + 'm' + ra[6:8] + 's'
                            decDegMinSec = dec[0:3] + 'd' + dec[4:6] + 'm' + dec[8:10] + 's'
                            
                            print'raHrMinSec  : ', raHrMinSec
                            print'decDegMinSec: ', decDegMinSec
                            
                            skyCoord = SkyCoord (raHrMinSec + ' ' + decDegMinSec, frame='icrs')

                            print 'skyCoord: ', skyCoord
                            
                            altAzi = skyCoord.transform_to(AltAz(obstime=observingNextTransitTime,location=observingPosition))

                            print 'altAzi: ', altAzi
                            print 'azi   : ', altAzi.az
                            print 'alt   : ', altAzi.alt
                            print
                            print 'file name             : ', file
                            print
                            
                            print 'dateTime              : ', dateTime
                            print 'dateTimeUTC           : ', dateTimeUTC
                            print
                            print 'System name           : ', root.findtext('name')
                            print 'System rightascension : ', root.findtext('rightascension')
                            print 'System declination    : ', root.findtext('declination')
                            print 'System magnitude      : ', mag
                            print
                            print 'Planet name           : ', planet.findtext('name')
                            print'Planet period         : ', planet.findtext('period')
                            print
                            print 'transitTimeBJD        : ', transitTimeBJD
                            print 'transitTime.jd        : ', transitTime.jd
                            print 'transitTime.fits      : ', transitTime.fits
                            print 'now                   : ', now
                            print 'now jd                : ', now.jd
                            print 'now fits              : ', now.fits
                            print 'delta                 : ', delta
                            print 'revolutionCount       : ', revolutionCount
                            print 'int revoultionCount   : ', int(revolutionCount) + 1
                            print 'nextTransit           : ', nextTransit
                            print 'nextTransitTime       : ', nextTransitTime.fits
                            print 'daysToTransit         : ', daysToTransit
                            print 'nextTransitTimePST    : ', nextTransitTimePST
                            print 'nTTPST.jd             : ', nTTPST.jd
                            print 'nTTPST.fits           : ', nTTPST.fits, 'PST'


                            print 'Star radius           : ', starRadius
                            print 'Planet radius         : ', planetRadius

                            print 'Planet/Star area ratio: ', planetStarAreaRatio
                            
                            print 'count: ', count
                            print
                            




        






