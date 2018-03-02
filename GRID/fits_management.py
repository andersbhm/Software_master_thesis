from astropy.io import fits
from astropy.table import Table
from ROOT import TFile, TTree
from array import array
import numpy as np
import os

def filesInFolderAsList(path):
    files_list = os.listdir(path)
    files_list.sort()
    return files_list

def fitsToRoot():

    # Create root file and tree
    f = TFile(path_root + "GRID_2008_2015.root", "RECREATE")
    tree = TTree("data_GRID","data_GRID")
    obt = array('d',[0])
    tree.Branch('obt', obt, 'obt/D')
    contact = array('I',[0])
    tree.Branch('contact', contact, 'contact/I')

    filenames = filesInFolderAsList(path)

    c = 0
    for name in filenames:
        c = c + 1
        print str(c) + "/" + str(len(filenames))

        # open fits file
        hdulist = fits.open(path + name)

        #hdulist.info()
        #print(hdulist['EVENTS'].columns)

        # Get table inside EVENTS
        event_data = Table(hdulist['EVENTS'].data)
        #print(event_data['TIME'])

        contact[0] = int(name[3:9])
        print contact[0]
        #Iterate through tabel and append TIME to tree
        for row in event_data:
            obt[0] = (row['TIME'])

            #if (354153660 - 60) < obt[0] and obt[0] < (362102460 + 60):
            tree.Fill()

        hdulist.close()




    #tree.Scan("","","col=20.6f")
    tree.Print()
    tree.Show(0)
    #tree.Scan("obt","", "col=20.6f")
    tree.GetCurrentFile().Write()
    tree.GetCurrentFile().Close()

path_root = "/Volumes/TOSHIBA/WWLLN_AGILE/GRID/rawdata/"
############################################
folder_name = "1.06.2008-27.02.2015"
#folder_name = "test"

#folder_name = "EVT_FM_files"
#folder_name = "EVT_FT3AB_files"
###########################################
path = path_root + folder_name + "/"
fitsToRoot()
