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

def draw_ellipse_on_png(im, x0, y0, ba, pa, major_axis_diameter_arcsec,
                        pixscale, color='#3388ff', linewidth=3):  #Draws an ellipse on the galaxy (from https://github.com/moustakas/legacyhalos/)
    Image.MAX_IMAGE_PIXELS = None
    
    minor_axis_diameter_arcsec = major_axis_diameter_arcsec * ba

    overlay_height = int(major_axis_diameter_arcsec / pixscale)
    overlay_width = int(minor_axis_diameter_arcsec / pixscale)
    overlay = Image.new('RGBA', (overlay_width, overlay_height))

    draw = ImageDraw.ImageDraw(overlay)
    box_corners = (0, 0, overlay_width, overlay_height)
    draw.ellipse(box_corners, fill=None, outline=color, width=linewidth)

    rotated = overlay.rotate(pa, expand=True)
    rotated_width, rotated_height = rotated.size
    paste_shift_x = int(x0 - rotated_width / 2)
    paste_shift_y = int(y0 - rotated_height / 2)
    im.paste(rotated, (paste_shift_x, paste_shift_y), rotated)


    with alive_bar(len(table)) as bar:
        for glx in table:
            if str(glx['sga_id']) == check or begun == True:
                begun = True
            if begun == True:
                
                pixscale = 0.1
                size = int(120*glx['d26']/pixscale)
                url = 'https://www.legacysurvey.org/viewer/cutout.jpg?ra=' +str(glx['ra'])+ '&dec=' +str(glx['dec'])+ '&layer=ls-dr9&pixscale='+str(pixscale)+'&size='+str(size)
                img = Image.open(urllib.request.urlopen(url,context = ctx))

                file = open(mainpath+r'\progress.txt','w')
                file.write(str(glx['sga_id']))
                file.close()

                if ellipse:
                    draw_ellipse_on_png(img,size/2-0.5,size/2-0.5,glx['ba'],glx['pa'],glx['d26']*60,pixscale)
                
                if show:
                    plt.imshow(img)
                    plt.show()

                if save:
                    file_name = str(glx['sga_id']) + ".jpg"
                    file_path = mainpath + r"\output" + "\\"
                    img.save(file_path+file_name)    
            bar()

    file = open(mainpath+r'\progress.txt','r')
    check = file.read()
    file.close()
    if str(table[len(table)-1]['sga_id']) == check:
        file = open(mainpath+r'\progress.txt','w')
        file.write('complete')
        file.close()
    

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