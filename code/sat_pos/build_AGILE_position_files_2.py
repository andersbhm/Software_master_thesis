# MM 13/02/2017 @UiB.  Build position file with satellite coordinates (X,Y,Z,lon,lat,h) as a function of time, based on TLEs
# Edited by Anders 10.05.17

# based on:
# build_AGILE_position_files.py
# check_AGILE_position_reconstruction.py

from sgp4.earth_gravity import wgs72
from sgp4.io import twoline2rv
from astropy.time import Time
from astropy import units as u
from astropy import coordinates as coord
from datetime import datetime, timedelta
import time
import numpy as np
from calendar import timegm
#import matplotlib.pyplot as plt
from scipy import mat, cos, sin, arctan, sqrt, pi, arctan2
from math import pow, degrees, radians, sqrt
from ROOT import TFile, TTree
from array import array

from astropy.utils.iers import IERS_A
iers_a = IERS_A.open('finals2000A.all')

import rotations as rot

def computation(name, start_year, start_month, start_day, start_hour, start_minute, start_second, end_year, end_month, end_day, end_hour, end_minute, end_second):

    # set initial and final times
    aepoch = time.strptime('2004 1 1 0 0 0', '%Y %m %d %H %M %S')
    tepoch = timegm(aepoch)

    dini = datetime(start_year, start_month, start_day, start_hour, start_minute, start_second) # start date for the analysis
    dfin = datetime(end_year, end_month, end_day, end_hour, end_minute, end_second) # stop date for the analysis
    dt = timedelta(0, 1) # time interval for orbital sampling (1 s ~ 8 km ~24us max timing error)

    #  create output root file and tree

    f = TFile(name, "recreate")
    tree = TTree("position","start of 2017 to 2017-11-09T10:09:10")
    year  = array('i',[0])
    month = array('i',[0])
    day   = array('i',[0])
    hour  = array('i',[0])
    min   = array('i',[0])
    sec   = array('i',[0])
    X = array('f',[0.])
    Y = array('f',[0.])
    Z  = array('f',[0.])
    lon = array('f',[0.])
    lat = array('f',[0.])
    h  = array('f',[0.])
    obt  = array('d',[0.])
    tree.Branch('year', year, 'year/I')
    tree.Branch('month', month, 'month/I')
    tree.Branch('day', day, 'day/I')
    tree.Branch('hour', hour, 'hour/I')
    tree.Branch('min', min, 'min/I')
    tree.Branch('sec', sec, 'sec/I')
    tree.Branch('X', X, 'X/F')
    tree.Branch('Y', Y, 'Y/F')
    tree.Branch('Z', Z, 'Z/F')
    tree.Branch('lon', lon, 'lon/F')
    tree.Branch('lat', lat, 'lat/F')
    tree.Branch('h', h, 'h/F')
    tree.Branch('obt', obt, 'obt/D')

    # load AGILE TLE

    agile_tle_file = open('AGILE_TLE_2017-Feb2018.dat', 'r') #data from 22 March 2015 (newconf data)
    agile_tle_list = agile_tle_file.readlines()
    agile_tle_list.reverse()

    # track AGILE orbit

    l5 = agile_tle_list.pop()
    l6 = agile_tle_list.pop()
    agile0 = twoline2rv(l5, l6, wgs72)
    l5 = agile_tle_list.pop()
    l6 = agile_tle_list.pop()
    agile1 = twoline2rv(l5, l6, wgs72)

    # loop on time range

    d = dini

    counter = 0
    while d <= dfin :


        while d > agile1.epoch :
            print ("From time ", d, " using TLE: ", l5)
            agile0 = agile1
            l3 = agile_tle_list.pop()
            l4 = agile_tle_list.pop()
            agile1 = twoline2rv(l3, l4, wgs72)

        xagile, vagile = agile0.propagate(d.year, d.month, d.day, d.hour, d.minute, d.second + d.microsecond * 1.e-6)
        x0 = np.array(xagile)
        print xagile
        X[0] = x0[0]
        Y[0] = x0[1]
        Z[0] = x0[2]
        year[0] = d.year
        month[0] = d.month
        day[0] = d.day
        hour[0] = d.hour
        min[0] = d.minute
        sec[0] = d.second

        t = Time(d)
        t.delta_ut1_utc = iers_a.ut1_utc(t)
        s = t.sidereal_time('mean', 'greenwich');
        s.wrap_angle=180 * u.deg
        r = rot.rotation_z(s)
        x0_ecef = r.dot(x0)
        lat[0], lon[0], h[0] = rot.ecef2geodetic(x0_ecef[0], x0_ecef[1], x0_ecef[2])
        obt[0] = timegm(d.timetuple()) - tepoch

        #print (d, year[0], month[0], day[0], hour[0], min[0], sec[0], obt[0], x0[0], x0[1], x0[2], lon[0], lat[0], h[0])



        counter = counter + 1
        if counter == 5000:
            counter = 0
            print (year[0], month[0], day[0], hour[0], min[0], sec[0])

        tree.Fill()
        d = d + dt

    tree.GetCurrentFile().Write()
    tree.GetCurrentFile().Close()

    '''
    # validation of root file

    TFile f("AGILE_position.root")
    TTree *t=f.Get("position")
    t->Scan("year:month:day:hour:min:sec:obt:X:Y:Z:lon:lat:h")
    t->Draw("lon:obt")
    t->Draw("10.*lat:obt","","same")
    t->Draw("0.1*h:obt","","same")
    '''

def main():

    #computation("name.root" ,start_year, start_month, start_day, start_hour, start_minute, start_second, end_year, end_month, end_day, end_hour, end_minute, end_second)

    #computation("AGILE_position_2017.root" , 2017, 1, 1, 0, 0, 0, 2017, 12, 31, 23, 59, 59)
    #computation("AGILE_position_2017.root" , 2017, 1, 1, 0, 0, 0, 2017, 11, 9, 10, 9, 10)
    #computation("AGILE_position_2017.root" , 2017, 11, 9, 10, 9, 11, 2017, 12, 31, 23, 59, 59)


    computation("AGILE_position_2018.root" , 2018, 1, 1, 0, 0, 0, 2018, 2, 8, 0, 0, 0)


main()
