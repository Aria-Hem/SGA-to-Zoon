import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from astropy.visualization import astropy_mpl_style
from astropy.table import vstack, Table
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
    return np.array(pandas.read_table(mainpath + path, comment='\\',delim_whitespace=True)) , pandas.read_table(mainpath + path, comment='#',delim_whitespace=True) 
    
def area_check(table,RA,DEC,Riso):
    N = len(table[DEC]) - 2
    in_area = np.zeros(N)
    sizes = np.array([])
    with alive_bar(N) as bar:
        for i in range(N):
            ra = float(table[RA][i+2]) * np.pi/180
            dec = float(table[DEC][i+2]) * np.pi/180
            center = 22.5 * np.pi/180
            radius = 10 * np.pi/180
            dist = np.arccos(np.sin(dec)*np.sin(center)+np.cos(dec)*np.cos(center)*np.cos(np.pi+ra))
            if dist<radius:
                in_area[i]=1
                sizes = np.append(sizes,float(table[Riso][i+2]))
            bar()
    return in_area , sizes

def area_size_check(table,RA,DEC,Riso,upper,lower):
    N = len(table[DEC]) - 2
    in_area = np.zeros(N)
    sizes = np.array([])
    with alive_bar(N) as bar:
        for i in range(N):
            ra = float(table[RA][i+2]) * np.pi/180
            dec = float(table[DEC][i+2]) * np.pi/180
            center = 22.5 * np.pi/180
            radius = 10 * np.pi/180
            dist = np.arccos(np.sin(dec)*np.sin(center)+np.cos(dec)*np.cos(center)*np.cos(np.pi+ra))
            if dist<radius and float(table[Riso][i+2]) <= upper and float(table[Riso][i+2]) >= lower:
                in_area[i]=1
            bar()
    return in_area


def match_for_size(table,RA,DEC,ref_table,in_area,fornames):
    margin = 1/360
    matching = Table({'ra':[0],'dec':[0],'NAME':[''],'NED_NAME':['']},names=('ra','dec','NAME','NED_NAME'))
    unmatching = Table({'ra':[0],'dec':[0],'NAME':[''],'NED_NAME':['']},names=('ra','dec','NAME','NED_NAME'))
    N = len(table[DEC]) - 2
    counter = 0
    with alive_bar(N) as bar:
        for i in range(N):
            if in_area[i]:
                c = 0
                ra = float(table[RA][i+2])
                dec = float(table[DEC][i+2])
                temp = Table({'ra':[ra],'dec':[dec],'NAME':[fornames[i+1][1]],'NED_NAME':[fornames[i+1][0]]},names=('ra','dec','NAME','NED_NAME'))
                for glx in ref_table:
                    if np.abs(glx['ra']-ra) <= margin and np.abs(glx['dec']-dec) <= margin :
                        counter += 1
                        c = 1
                        matching = vstack([matching,temp])
                        break
                if c==0:
                    unmatching = vstack([unmatching,temp])
            bar()
    matching.remove_row(0)
    unmatching.remove_row(0)
    return matching, unmatching ,counter



def cross_matching(table, RA, DEC, ref_table, in_area):
    margin = 1/360
    N = len(table[DEC]) - 2
    counter = 0
    sizes = np.array([])
    with alive_bar(N) as bar:
        for i in range(N):
            if in_area[i]:
                ra = float(table[RA][i+2])
                dec = float(table[DEC][i+2])
                for glx in ref_table:
                    if np.abs(glx['ra']-ra) <= margin and np.abs(glx['dec']-dec) <= margin :
                        counter += 1
                        sizes = np.append(sizes,float(table[Riso][i+2]))
                        break
            bar()
    return counter, sizes

ref_table = openFITS(r'sga-query.fits')
RA = '\\'
DEC = 'WXSC'
Riso = '(jarrett'
fornames, table = opentbl(r'\WXSC_Riso_1arcmin_10Jun2024.tbl')

print(Table(openFITS('unmatching.fits')))

in_area = area_size_check(table,RA,DEC,Riso,10000,0)

ra = np.array([])
dec = np.array([])
for i in range(len(table[DEC])-2):
    if in_area[i] == 1:
        ra = np.append(ra,float(table[RA][i+2]))
        dec = np.append(dec,float(table[DEC][i+2]))
plt.scatter(ra,dec, marker = '.')
plt.title('Coordinate of the Galaxies in the FITS file')
plt.xlabel('RA (deg)')
plt.ylabel('DEC (deg)')
plt.show()

matching, unmatching , counter = match_for_size(table,RA,DEC,ref_table,in_area,fornames)
matching.write('matching.fits', format='fits',overwrite=True)
unmatching.write('unmatching.fits', format='fits',overwrite=True)


data_coordinate_plot(matching)
data_coordinate_plot(unmatching)

in_area , all_sizes = area_check(table,RA,DEC,Riso)
counter, sizes = cross_matching(table,RA,DEC,ref_table, in_area)

total_in_area = np.sum(in_area)
print(counter,total_in_area,100 * counter/total_in_area)

bins = np.linspace(30, 600, 50)

plt.hist(all_sizes, bins, alpha=0.5, label='WXSC')
plt.hist(sizes, bins, alpha=0.5, label='SGA')
plt.legend(loc='upper right')
plt.title('Histogram of Riso')
plt.xlabel('Riso (arcsec)')
plt.show()