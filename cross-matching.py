import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from astropy.visualization import astropy_mpl_style
import urllib
import pathlib
from PIL import Image, ImageDraw, ImageFont
from alive_progress import alive_bar
import time
import pandas

#Initial variables and actions

mainpath = str(pathlib.Path(__file__).parent.resolve())


def data_coordinate_plot(table): #Ploting RA and DEC of the Galaxies
    RA = table['ra']
    DEC = table.field('dec')

    plt.scatter(RA,DEC, marker = '.')
    plt.title('Coordinate of the Galaxies in the FITS file')
    plt.xlabel('RA (deg)')
    plt.ylabel('DEC (deg)')
    plt.show()


def openFITS(path , hdu = 1): #Opening the FITS file
    image_file = get_pkg_data_filename(path)
    fits.info(image_file)
    fitsfile = fits.open(image_file)
    return fitsfile[hdu].data

def opentbl(path):
    return pandas.read_table(mainpath + path, comment='#', delim_whitespace=True)
    

ref_table = openFITS(r'sga-query.fits')
RA = '\\'
DEC = 'WXSC'
table = opentbl(r'\WXSC_Riso_1arcmin_10Jun2024.tbl')

N = len(table[DEC]) - 2
#print( type(table[DEC][2]) , len(table[DEC]) )
in_area = np.zeros(N)

with alive_bar(N) as bar:
    for i in range(N):
        ra = float(table[RA][i+2]) * np.pi/180
        dec = float(table[DEC][i+2]) * np.pi/180
        center = 22.5 * np.pi/180
        radius = 10 * np.pi/180
        dist = np.arccos(np.cos(dec)*np.cos(center)+np.sin(dec)*np.sin(center)*np.cos(np.pi+ra))
        if dist<radius :
            in_area[i]=1
        bar()

total_in_area = np.sum(in_area)
counter = 0
margin = 1/360

with alive_bar(N) as bar:
    for i in range(N):
        if in_area[i]:
            ra = float(table[RA][i+2])
            dec = float(table[DEC][i+2])
            for glx in ref_table:
                if np.abs(glx['ra']-ra) <= margin and np.abs(glx['dec']-dec) <= margin :
                    counter += 1
                    break
        bar()   

print(counter,total_in_area,100 * counter/total_in_area)
#data_coordinate_plot(table)