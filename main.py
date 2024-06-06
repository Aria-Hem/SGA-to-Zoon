import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from astropy.visualization import astropy_mpl_style
import urllib
import pathlib
from PIL import Image
from alive_progress import alive_bar
path = str(pathlib.Path(__file__).parent.resolve())

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

def draw_ellipse_on_png(im, x0, y0, ba, pa, major_axis_diameter_arcsec,
                        pixscale, color='#3388ff', linewidth=3):  #Draws an ellipse on the galaxy (from https://github.com/moustakas/legacyhalos/)

    from PIL import Image, ImageDraw, ImageFont

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

def galaxy_image_exporter(table, ellipse = True, show = False, save = True):  #Downloading a cutout of each galaxy from www.legacysurvey.org
    for glx in table:
        pixscale = 0.1
        size = int(120*glx['d26']/pixscale)
        url = 'https://www.legacysurvey.org/viewer/cutout.jpg?ra=' +str(glx['ra'])+ '&dec=' +str(glx['dec'])+ '&layer=ls-dr9&pixscale='+str(pixscale)+'&size='+str(size)
        img = Image.open(urllib.request.urlopen(url))
        if ellipse:
            draw_ellipse_on_png(img,size/2-0.5,size/2-0.5,glx['ba'],glx['pa'],glx['d26']*60,pixscale)
        if show:
            plt.imshow(img)
            plt.show()
        if save:
            file_name = str(glx['sga_id']) + ".jpg"
            file_path = path + r"\output" + "\\"
            img.save(file_path+file_name)
        yield


table = openFITS(r'sga2020-2masx.fits')
#data_coordinate_plot(table)
with alive_bar(len(table)) as bar:
    for i in galaxy_image_exporter(table):
        bar()