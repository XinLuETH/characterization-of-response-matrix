import numpy as np

#import os
from astropy.io import fits


path = "D:\Xin LU\\real response matrix\\9X9\\"
darks = "D:\Xin LU\\real response matrix\\9X9\\dark_50ms.fits"
outpath = "D:\Xin LU\\real response matrix\\9X9\\processed\\"

zernike_num = 28  # total number of used Zernike terms 11   15
zernike_coef = 5  # step of Zernike term  5 nm
phase_shift = 2   # for phase shifting method:4    for standard Zernike sensing: 2


dark_h = fits.getheader(darks)
t = dark_h['EXPTIME']
dark_data = np.median(fits.getdata(darks)[6,:,:], axis=0)  # discard the first five frames


for coef_i in range(9):
    if coef_i<4:
        file_name1 = 'linear11_N'+str(-20+coef_i*zernike_coef)[1:]
    else:
        file_name1 = 'linear11_'+str(-20+coef_i*zernike_coef)
    for zernike_i in range(2,zernike_num+1):
        file_name2 = '_Z'
        file_name2 = file_name2 + str(zernike_i)+'_'
        for phase_i in range(phase_shift):
            file_name3 = str(phase_i)+'Pi.fits'
            filename = file_name1+file_name2+file_name3
            file_name = path + filename

            h = fits.getheader(file_name)
            data = fits.getdata(file_name)[6:,:,:]
            
            t = h['EXPTIME']
            data = data - dark_data
            data = np.median(data, axis=0)  #
            
            hdu = fits.PrimaryHDU(data=data)
            hdu.header['EXPTIME'] = t
            # xx = 'standard'+filename[5:]

            # filename =xx
            file_nameS = outpath+filename
            hdu.writeto(file_nameS)
#dark = {}
#for d in [os.path.join(darks, f) for f in os.listdir(darks)]:
#    h = fits.getheader(d)
#    t = h['EXPTIME']
#    data = np.median(fits.getdata(d), axis=0)
#    dark[t] = data
#
#print("Processed darks")
#
#files = [os.path.join(path, f) for f in os.listdir(path)]

#for file in files:
#    # get integration time
#    h = fits.getheader(file)
#    data = fits.getdata(file)
#
#    t = h['EXPTIME']
#    data = data - dark[t]
#    data = np.median(data, axis=0)
#
#    filename = os.path.basename(file).split(sep='.')[0]
#    filename += 'PAL11lD_processed.fits'
#
#    hdu = fits.PrimaryHDU(data=data)
#    hdu.header['EXPTIME'] = t
#    hdu.writeto(os.path.join(outpath,filename))
#    print('Processed file %s' %filename)
