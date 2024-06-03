import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from astropy.visualization import astropy_mpl_style
import urllib
import pathlib
path = str(pathlib.Path(__file__).parent.resolve())

#Opening the FITS file
image_file = get_pkg_data_filename(r'sga2020-2masx.fits')

#data_list = fits.getdata(image_file, ext=0)
fits.info(image_file)
fitsfile = fits.open(image_file)
hdu = 1
table = fitsfile[hdu].data

#Ploting RA and DEC of the Galaxies
RA = table['ra']
DEC = table.field('dec')
plt.scatter(RA,DEC, marker = '.')
plt.title('Coordinate of the Galaxies in the FITS file')
plt.xlabel('RA (deg)')
plt.ylabel('DEC (deg)')
#plt.show()

#Downloading a cutout of each galaxy from www.legacysurvey.org
for i in table:
    url = 'https://www.legacysurvey.org/viewer/cutout.jpg?ra=' +str(i[1])+ '&dec=' +str(i[2])+ '&layer=ls-dr9&pixscale=0.1&size='+str(int(600*i['d26']*2))
    print(url)
    file_name = str(i[6]) + ".jpg"
    file_path = path + r"\output" + "\\"
    print(file_path+file_name)
    urllib.request.urlretrieve(url, file_path+file_name)