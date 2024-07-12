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
import pandas as pd

#Initial variables and actions

mainpath = str(pathlib.Path(__file__).parent.resolve())


def data_coordinate_plot(table): #Ploting RA and DEC of the Galaxies
    RA = table['RA']
    DEC = table.field('DEC')

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
    return np.array(pd.read_table(mainpath + path, comment='\\',delim_whitespace=True)) , pd.read_table(mainpath + path, comment='#',delim_whitespace=True) 
    
def keep_first_numeric(data):
    temp=""
    error = 0
    minus_count = 0
    dot_count = 0
    for c in data:
        if c == '-' :
            minus_count+=1
        if c == '.' :
            dot_count+=1
        if c.isnumeric() or (c == '-' and minus_count<2) or (c == '.' and dot_count<2):
            temp+=c
        else:
            break
    if minus_count>1 or dot_count>1 or temp == '' :
        error = 1
    return str(temp) , error

def get_from_csv(name):
    df = pd.read_csv(mainpath+name,encoding='cp1252', header=0)
    colname = " pinned_globally"
    df = pd.DataFrame(df)
    print(type(df['id'][10]))
    df['sum'] = df.apply(lambda x: str(x['id'])+str(x[colname])+' &',axis=1)
    links = df['sum']
    return links


def data_extraction(links):
    RA = np.array([])
    DEC = np.array([])
    IDs = np.array([])
    Zoom = np.array([])
    Layer = np.array([])
    for row in links:
        url = str(row)
        id,error_id = keep_first_numeric(url)
        ind = row.find("legacysurvey.org")
        ra_key="ra="
        dec_key="dec="
        zoom_key = "zoom="
        layer_key = "layer="
        if ind != -1:
            ind_ra = url.find(ra_key)
            ind_ra_end = url.find("&",ind_ra)

            ind_dec = url.find(dec_key)
            ind_dec_end = url.find("&",ind_dec)

            ind_zoom = url.find(zoom_key)

            ind_layer = url.find(layer_key)
            ind_layer_end = np.min([url.find(" ",ind_layer),url.find("&",ind_layer)])

            if ind_ra != -1 and ind_dec != -1 :
                ra , error_ra = keep_first_numeric(url[ind_ra+len(ra_key):ind_ra_end])
                dec , error_dec = keep_first_numeric(url[ind_dec+len(dec_key):ind_dec_end])
                zoom , error_zoom = keep_first_numeric(url[ind_zoom+len(zoom_key):-1])
                layer = url[ind_layer+len(layer_key):ind_layer_end]
                if ind_layer == -1:
                    layer='' 
                if error_ra==0 and error_dec==0:
                    RA = np.append(RA,float(ra))
                    DEC = np.append(DEC,float(dec))
                    Zoom = np.append(Zoom,zoom)
                    IDs = np.append(IDs,id)
                    Layer = np.append(Layer,layer)
    return RA,DEC,Zoom,Layer, IDs

links_col = get_from_csv(r'\topics-publicsv.csv')
ra, dec, zoom, layer, id = data_extraction(links_col)

extracted_data = Table({'Topic id':id,'RA':ra,'DEC':dec,'zoom':zoom,'layer':layer},names=('Topic id','RA','DEC','zoom','layer'))
#extracted_data.sort(['RA','DEC'])
extracted_data.write(mainpath + '\\' +'extracted_data.fits', format='fits',overwrite=True)

print(extracted_data)
data_coordinate_plot(extracted_data)



