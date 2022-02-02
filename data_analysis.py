
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from aplpy import FITSFigure
import matplotlib.colors as colors
from matplotlib.lines import Line2D
from scipy.stats import norm
from scipy.stats import kde
import statsmodels.api as sm 
from subroutines_mapmaking import readfilesfinal
from subroutines_mapmaking import make_p_theta
from subroutines_mapmaking import wrapangle
from subroutines_mapmaking import polangle_diff
from subroutines_mapmaking import plot_I_polvec
import pdb
import copy

pixsize = '8.0' # '4.0' # 


### READ AND SET UP MAPS ###

### ORIGINAL IQUs ###
files_d_orig = ['N2071_HAWC+D_i+var_Jy-arcs2.fits',  # Filenames
                  'N2071_HAWC+D_q+var_Jy-arcs2.fits',
                  'N2071_HAWC+D_u+var_Jy-arcs2.fits'] 
files_e_orig = ['N2071_HAWC+E_i+var_Jy-arcs2.fits', 
                  'N2071_HAWC+E_q+var_Jy-arcs2.fits', 
                  'N2071_HAWC+E_u+var_Jy-arcs2.fits']
files_850_orig = ['maps_from_2020-07_JCMT_data/N2071_JCMT-850_i+var_Jy-arcs2.fits', 
                  'maps_from_2020-07_JCMT_data/N2071_JCMT-850_q+var_Jy-arcs2.fits', 
                  'maps_from_2020-07_JCMT_data/N2071_JCMT-850_u+var_Jy-arcs2.fits']
''' (Old pipeline)
files_d_orig = ['HAWC+_pics_old/Old_pipeline/N2071_HAWC+D_i+var_Jy-arcs2.fits',
                  'HAWC+_pics_old/Old_pipeline/N2071_HAWC+D_q+var_Jy-arcs2.fits',
                  'HAWC+_pics_old/Old_pipeline/N2071_HAWC+D_u+var_Jy-arcs2.fits'] 
files_e_orig = ['HAWC+_pics_old/Old_pipeline/N2071_HAWC+E_i+var_Jy-arcs2.fits', 
                  'HAWC+_pics_old/Old_pipeline/N2071_HAWC+E_q+var_Jy-arcs2.fits', 
                  'HAWC+_pics_old/Old_pipeline/N2071_HAWC+E_u+var_Jy-arcs2.fits']
files_850_orig = ['N2071_JCMT-850_i+var_Jy-arcs2.fits', 
                    'N2071_JCMT-850_q+var_Jy-arcs2.fits', 
                    'N2071_JCMT-850_u+var_Jy-arcs2.fits']
'''
i_d_orig, q_d_orig, u_d_orig, di_d_orig, dq_d_orig, du_d_orig = readfilesfinal(files_d_orig, units_in = 'Jy/arcsec2')  # Reading + unit conversion
i_e_orig, q_e_orig, u_e_orig, di_e_orig, dq_e_orig, du_e_orig = readfilesfinal(files_e_orig, units_in = 'Jy/arcsec2')
i_850_orig, q_850_orig, u_850_orig, di_850_orig, dq_850_orig, du_850_orig = readfilesfinal(files_850_orig, units_in = 'Jy/arcsec2')


### RESAMPLED BUT UNSMOOTHED IQUs ###
files_d_rs = ['N2071_HAWC+D_i+var_Jy-arcs2_pix' + pixsize +'-NN.fits',
                  'N2071_HAWC+D_q+var_Jy-arcs2_pix' + pixsize +'-NN.fits',
                  'N2071_HAWC+D_u+var_Jy-arcs2_pix' + pixsize +'-NN.fits'] 
files_e_rs = ['N2071_HAWC+E_i+var_Jy-arcs2_pix' + pixsize +'-NN_beam18.9.fits',      # E band: same as smoothed file
                  'N2071_HAWC+E_q+var_Jy-arcs2_pix' + pixsize +'-NN_beam18.9.fits', 
                  'N2071_HAWC+E_u+var_Jy-arcs2_pix' + pixsize +'-NN_beam18.9.fits']
if pixsize == '4.0':
    files_850_rs = ['maps_from_2020-07_JCMT_data/N2071_JCMT-850_i+var_Jy-arcs2.fits', 
                    'maps_from_2020-07_JCMT_data/N2071_JCMT-850_q+var_Jy-arcs2.fits', 
                    'maps_from_2020-07_JCMT_data/N2071_JCMT-850_u+var_Jy-arcs2.fits']
elif pixsize == '8.0':
    files_850_rs = ['N2071_JCMT-850-8as_i+var_Jy-arcs2.fits',
                    'N2071_JCMT-850-8as_q+var_Jy-arcs2.fits',
                    'N2071_JCMT-850-8as_u+var_Jy-arcs2.fits']
else:
    pass
i_d_rs, q_d_rs, u_d_rs, di_d_rs, dq_d_rs, du_d_rs = readfilesfinal(files_d_rs, units_in = 'Jy/arcsec2')
i_e_rs, q_e_rs, u_e_rs, di_e_rs, dq_e_rs, du_e_rs = readfilesfinal(files_e_rs, units_in = 'Jy/arcsec2')
i_850_rs, q_850_rs, u_850_rs, di_850_rs, dq_850_rs, du_850_rs = readfilesfinal(files_850_rs, units_in = 'Jy/arcsec2')


### RESAMPLED + SMOOTHED IQUs ###
## New pipeline
files_d_smooth = ['N2071_HAWC+D_i+var_Jy-arcs2_pix' + pixsize +'-NN_beam18.9.fits',
                  'N2071_HAWC+D_q+var_Jy-arcs2_pix' + pixsize +'-NN_beam18.9.fits',
                  'N2071_HAWC+D_u+var_Jy-arcs2_pix' + pixsize +'-NN_beam18.9.fits'] 
files_e_smooth = ['N2071_HAWC+E_i+var_Jy-arcs2_pix' + pixsize +'-NN_beam18.9.fits', 
                  'N2071_HAWC+E_q+var_Jy-arcs2_pix' + pixsize +'-NN_beam18.9.fits', 
                  'N2071_HAWC+E_u+var_Jy-arcs2_pix' + pixsize +'-NN_beam18.9.fits']
files_850_smooth = ['N2071_JCMT-850_i+var_Jy-arcs2_beam18.9_pix' + pixsize +'.fits',   # JCMT files have no regridding anyway
                    'N2071_JCMT-850_q+var_Jy-arcs2_beam18.9_pix' + pixsize +'.fits', 
                    'N2071_JCMT-850_u+var_Jy-arcs2_beam18.9_pix' + pixsize +'.fits']
''' Old pipeline (rebin rather than resample)
files_d_smooth = ['HAWC+_pics_old/Old_pipeline/N2071_HAWC+D_i+var_Jy-arcs2_beam18.2_pix4.0-sincsinc-rb.fits',
                  'HAWC+_pics_old/Old_pipeline/N2071_HAWC+D_q+var_Jy-arcs2_beam18.2_pix4.0-sincsinc-rb.fits',
                  'HAWC+_pics_old/Old_pipeline/N2071_HAWC+D_u+var_Jy-arcs2_beam18.2_pix4.0-sincsinc-rb.fits'] 
files_e_smooth = ['HAWC+_pics_old/Old_pipeline/N2071_HAWC+E_i+var_Jy-arcs2_beam18.2_pix4.0-sincsinc-rb.fits', 
                  'HAWC+_pics_old/Old_pipeline/N2071_HAWC+E_q+var_Jy-arcs2_beam18.2_pix4.0-sincsinc-rb.fits', 
                  'HAWC+_pics_old/Old_pipeline/N2071_HAWC+E_u+var_Jy-arcs2_beam18.2_pix4.0-sincsinc-rb.fits']
files_850_smooth = ['N2071_JCMT-850_i+var_Jy-arcs2_beam18.2_pix4.0.fits',   # JCMT files have no regridding anyway
                    'N2071_JCMT-850_q+var_Jy-arcs2_beam18.2_pix4.0.fits', 
                    'N2071_JCMT-850_u+var_Jy-arcs2_beam18.2_pix4.0.fits']
'''
i_d, q_d, u_d, di_d, dq_d, du_d = readfilesfinal(files_d_smooth, units_in = 'Jy/arcsec2')
i_e, q_e, u_e, di_e, dq_e, du_e = readfilesfinal(files_e_smooth, units_in = 'Jy/arcsec2')
i_850, q_850, u_850, di_850, dq_850, du_850 = readfilesfinal(files_850_smooth, units_in = 'Jy/arcsec2')


## ANCILLARY (HERSCHEL) DATA
## Temperature
fname_T = 'Ancillary_Herschel_data/N2071_T_cutout_rebin'
hT = fits.open(fname_T + pixsize + '.fits')
T_map = hT[0]
hT = fits.open(fname_T + '2rawd.fits')     # Original D-band pixels
T_map_origd = hT[0]                    
hT = fits.open(fname_T + '2rawe.fits')     # Original E-band pixels
T_map_orige = hT[0]                    

# Optical depth
fname_tau_d = 'Ancillary_Herschel_data/N2071_tau154_cutout_rebin'
htau_d = fits.open(fname_tau_d + pixsize + '.fits')
tau_map_d = htau_d[0]
htau_d = fits.open(fname_tau_d + '2rawd.fits')   # Original D-band pixels
tau_map_d_orig = htau_d[0]

fname_tau_e = 'Ancillary_Herschel_data/N2071_tau214_cutout_rebin'
htau_e = fits.open(fname_tau_e + pixsize + '.fits')
tau_map_e = htau_e[0]
htau_e = fits.open(fname_tau_e + '2rawe.fits')   # Original E-band pixels
tau_map_e_orig = htau_e[0]

# H2 back-converted from tau
m_H = 1.67e-24      # H mass
mu = 2.8            # Mol. weight for Orion B (Konyves et al. 2020)
wl_ref = 300.       # Ref. wavelength in um (Palmeirin et al. 2013)
op_ref = 0.1        # Ref. opacity in cm2/g (of dust+gas)
beta = 2.           # Beta used for Herschel fits
wl_used = 154.      # Wavelength we want
op_conv_fact = m_H * mu * op_ref * (wl_ref/wl_used)**beta
# Final grid
temp = copy.deepcopy(tau_map_d.data)
temp /= op_conv_fact
nh2_map = fits.PrimaryHDU()
nh2_map.data = temp
nh2_map.header = tau_map_d.header
# Original D-band grid
temp = copy.deepcopy(tau_map_d_orig.data)
temp /= op_conv_fact
nh2_map_origd = fits.PrimaryHDU()
nh2_map_origd.data = temp
nh2_map_origd.header = tau_map_d_orig.header

### ###

## READ 
#i_d, q_d, u_d, di_d, dq_d, du_d = readfilesfinal(files_d_smooth, units_in = 'Jy/arcsec2')
#i_e, q_e, u_e, di_e, dq_e, du_e = readfilesfinal(files_e_smooth, units_in = 'Jy/arcsec2')
#i_850, q_850, u_850, di_850, dq_850, du_850 = readfilesfinal(files_850_smooth, units_in = 'Jy/arcsec2')

## MAKE POLARIZED QUANTITIES
# Original format
Ip_d_orig, p_d_orig, theta_d_orig, dIp_d_orig, dp_d_orig, dtheta_d_orig = make_p_theta(i_d_orig, q_d_orig, u_d_orig, di_d_orig, dq_d_orig, du_d_orig, 
                                                                                       magfield = 'False')
Ip_e_orig, p_e_orig, theta_e_orig, dIp_e_orig, dp_e_orig, dtheta_e_orig = make_p_theta(i_e_orig, q_e_orig, u_e_orig, di_e_orig, dq_e_orig, du_e_orig,
                                                                                       magfield = 'False')
Ip_850_orig, p_850_orig, theta_850_orig, dIp_850_orig, dp_850_orig, dtheta_850_orig = make_p_theta(i_850_orig, q_850_orig, u_850_orig, 
                                                                                                   di_850_orig, dq_850_orig, du_850_orig, magfield = 'False')
# Resampled but unsmoothed
Ip_d_rs, p_d_rs, theta_d_rs, dIp_d_rs, dp_d_rs, dtheta_d_rs = make_p_theta(i_d_rs, q_d_rs, u_d_rs, di_d_rs, dq_d_rs, du_d_rs, magfield = 'False')
Ip_e_rs, p_e_rs, theta_e_rs, dIp_e_rs, dp_e_rs, dtheta_e_rs = make_p_theta(i_e_rs, q_e_rs, u_e_rs, di_e_rs, dq_e_rs, du_e_rs, magfield = 'False')
Ip_850_rs, p_850_rs, theta_850_rs, dIp_850_rs, dp_850_rs, dtheta_850_rs = make_p_theta(i_850_rs, q_850_rs, u_850_rs, di_850_rs, dq_850_rs, du_850_rs, magfield = 'False')
# Smoothed + resampled
Ip_d, p_d, theta_d, dIp_d, dp_d, dtheta_d = make_p_theta(i_d, q_d, u_d, di_d, dq_d, du_d, magfield = 'False')#, save='N2071_HAWC+D')
Ip_e, p_e, theta_e, dIp_e, dp_e, dtheta_e = make_p_theta(i_e, q_e, u_e, di_e, dq_e, du_e, magfield = 'False')#, save='N2071_HAWC+E')
Ip_850, p_850, theta_850, dIp_850, dp_850, dtheta_850 = make_p_theta(i_850, q_850, u_850, di_850, dq_850, du_850, magfield = 'False')

# Angle differences (smoothed + regridded)
angles_d = copy.deepcopy(theta_d)
angles_e = copy.deepcopy(theta_e)
angles_850 = copy.deepcopy(theta_850)
dangles_d = copy.deepcopy(dtheta_d)
dangles_e = copy.deepcopy(dtheta_e)
dangles_850 = copy.deepcopy(dtheta_850)
delta_theta_de = polangle_diff(angles_d.data, angles_e.data, liminf = -90., limsup = 90.)
#delta_theta_de_masked = np.ma.masked_where(mask_de, delta_theta_de)
delta_theta_SNR_de = np.abs(delta_theta_de) / np.sqrt(dangles_d.data**2 + dangles_e.data**2)
delta_theta_d850 = polangle_diff(angles_d.data, angles_850.data, liminf = -90., limsup = 90.)
#delta_theta_d850_masked = np.ma.masked_where(mask_d850, delta_theta_d850)
delta_theta_SNR_d850 = np.abs(delta_theta_d850) / np.sqrt(dangles_d.data**2 + dangles_850.data**2)
delta_theta_e850 = polangle_diff(angles_e.data, angles_850.data, liminf = -90., limsup = 90.)
delta_theta_SNR_e850 = np.abs(delta_theta_e850) / np.sqrt(dangles_e.data**2 + dangles_850.data**2)


### MASKS ###

## Masking by position
## Smoothed maps
xlen = i_850.header['NAXIS1']
ylen = i_850.header['NAXIS2']
# Distance from center (off by 1?)
xc = i_850.header['CRPIX1']
yc = i_850.header['CRPIX2']
x, y = np.meshgrid(np.arange(xlen), np.arange(ylen))
distance_array = np.sqrt(((x-xc+1) * i_850.header['CDELT1'])**2 + ((y-yc+1) * i_850.header['CDELT2'])**2)
distance_array_x = np.abs(((x-xc+1) * i_850.header['CDELT1']))
distance_array_y = np.abs(((y-yc+1) * i_850.header['CDELT2']))
if pixsize == '4.0':
    # Distance from emission peak
    #peakpos = np.argwhere(np.rot90(i_850.data.T) > 0.999*np.nanmax(i_850.data)) # Array of max position
    #peakpos = np.argwhere(i_850.data > 0.999*np.nanmax(i_850.data)) # Array of max position
    #xpeak = peakpos[0, 0]   # If more than 1 candidate, we take the 1st
    #ypeak = peakpos[0, 1]
    xlen = i_850.header['NAXIS1']
    ylen = i_850.header['NAXIS2']
    xpeak = 129
    ypeak = 128
    x, y = np.meshgrid(np.arange(xlen), np.arange(ylen))
    peak_distance_array = np.sqrt(((x-xpeak) * i_850.header['CDELT1'])**2 + ((y-ypeak) * i_850.header['CDELT2'])**2)
    # Original-pixel distance array
    distance_array_orig = distance_array
if pixsize == '8.0':
    # Distance from emission peak
    xlen = i_850.header['NAXIS1']
    ylen = i_850.header['NAXIS2']
    xpeak = 65
    ypeak = 64
    x, y = np.meshgrid(np.arange(xlen), np.arange(ylen))
    peak_distance_array = np.sqrt(((x-xpeak) * i_850.header['CDELT1'])**2 + ((y-ypeak) * i_850.header['CDELT2'])**2)
    # Original-pixel distance array
    xlen = i_850_orig.header['NAXIS1']
    ylen = i_850_orig.header['NAXIS2']
    xc = i_850_orig.header['CRPIX1']
    yc = i_850_orig.header['CRPIX2']
    x, y = np.meshgrid(np.arange(xlen), np.arange(ylen))
    distance_array_orig = np.sqrt(((x-xc) * i_850_orig.header['CDELT1'])**2 + ((y-yc) * i_850_orig.header['CDELT2'])**2)
else:
    pass
# Original D-band pixels
# Distance from emission peak
#  (TO DO)
# Distance from center
xlen = i_d_orig.header['NAXIS1']
ylen = i_d_orig.header['NAXIS2']
xc = i_d_orig.header['CRPIX1']
yc = i_d_orig.header['CRPIX2']
x, y = np.meshgrid(np.arange(xlen), np.arange(ylen))
distance_array_orig_d = np.sqrt(((x-xc+1) * i_d_orig.header['CDELT1'])**2 + ((y-yc+1) * i_d_orig.header['CDELT2'])**2)
# Original E-band pixels
xlen = i_e_orig.header['NAXIS1']
ylen = i_e_orig.header['NAXIS2']
xc = i_e_orig.header['CRPIX1']
yc = i_e_orig.header['CRPIX2']
x, y = np.meshgrid(np.arange(xlen), np.arange(ylen))
distance_array_orig_e = np.sqrt(((x-xc+1) * i_e_orig.header['CDELT1'])**2 + ((y-yc+1) * i_e_orig.header['CDELT2'])**2)

'''
## Testing for off-by-one
np.argwhere(i_850.data > 0.999*np.nanmax(i_850.data))
w = WCS(i_850.header)
# CRVAL1, CRVAL2:  86.771,      0.364361111111108
# corresponding to w.pixel_to_world(CRPIX1 - 1,  CRPIX2 - 1)
w.pixel_to_world(66,  66)  # Res: 86.76877773, 0.36658333 (for 8 arcsec pixels)
w.pixel_to_world(65,  64)  # Res: 86.771,      0.36213889 (peak, 8 arcsec)
# To invert:
from astropy.coordinates import SkyCoord
c = SkyCoord(86.771, .364, unit="deg", frame="icrs")
w.world_to_pixel(c)

distplot = fits.PrimaryHDU()
distplot.data = distance_array
distplot.header = i_850.header

i2plot = distplot # i_850 # 
fig_test = FITSFigure(i2plot, subplot=(1,1,1))
fig_test.axis_labels.set_font(size='xx-large')
fig_test.tick_labels.set_font(size='xx-large')
fig_test.show_colorscale(cmap='rainbow', vmin = 0., vmax = .03)
fig_test.show_markers(86.771, 0.36436, marker = "o", facecolor = 'white', edgecolor = 'black')
fig_test.add_colorbar()
fig_test.colorbar.set_font(size='xx-large')
fig_test.colorbar.set_axis_label_font(size='xx-large')
fig_test.show_contour(distance_array, linestyles = '--', linewidths = 2., levels = [45./3600.])
x_center = i2plot.header['CRVAL1']
y_center = i2plot.header['CRVAL2']
fig_test.recenter(x_center, y_center, radius = 90./3600.)

# Testing rotation:
#  plt.imshow(i_850.data)  # Shows x-y swap (+ rotation?)
#  irot = np.rot90(i_850.data.T)
#  plt.imshow(irot)
'''

## Masking by data quality
## Smoothed data

# All criteria (nota: if P = NaN, then I = NaN)
# P / dP > 3
mask_d = np.any([p_d.data/dp_d.data < 3, i_d.data / di_d.data < 20, i_d.data < 0.025, dp_d.data > 5., np.isnan(p_d.data), np.isnan(dp_d.data),
                 distance_array > 300. / 3600.], axis = 0)
mask_e = np.any([p_e.data/dp_e.data < 3, i_e.data / di_e.data < 20, i_e.data < 0.015, dp_e.data > 5, np.isnan(p_e.data), np.isnan(dp_e.data), 
                 distance_array > 300. / 3600.], axis = 0)
mask_850 = np.any([p_850.data/dp_850.data < 3, i_850.data / di_850.data < 20, i_850.data < 0.00025, dp_e.data > 5., np.isnan(p_850.data), 
                   np.isnan(dp_850.data), distance_array > 180. / 3600.], axis = 0)
mask_de = np.any([mask_d, mask_e], axis = 0)
mask_d850 = np.any([mask_d, mask_850], axis = 0)
mask_e850 = np.any([mask_850, mask_e], axis = 0)
mask_all = np.any([mask_d, mask_e, mask_850], axis = 0)
# Additional: focus on central 90''
mask_de_90as = np.any([mask_de, peak_distance_array > 45. / 3600.], axis = 0)    # For statistical purposes
mask_850_90as = np.any([mask_850, peak_distance_array > 45. / 3600.], axis = 0)
#mask_850_90as = np.any([p_850.data/dp_850.data < 3, i_850.data / di_850.data < 20, i_850.data < 0.00025, dp_e.data > 5., np.isnan(p_850.data),
#                   np.isnan(dp_850.data), peak_distance_array > 45. / 3600.], axis = 0)   
mask_d850_90as = np.any([mask_d, mask_850_90as], axis = 0)
mask_e850_90as = np.any([mask_850_90as, mask_e], axis = 0)
mask_all_90as = np.any([mask_d, mask_e, mask_850_90as], axis = 0)
# P / dP > 2
mask_d_psnr2 = np.any([p_d.data/dp_d.data < 2, i_d.data / di_d.data < 20, i_d.data < 0.025, dp_d.data > 5., np.isnan(p_d.data), np.isnan(dp_d.data), 
                       distance_array > 300. / 3600.], axis = 0)
mask_e_psnr2 = np.any([p_e.data/dp_e.data < 2, i_e.data / di_e.data < 20, i_e.data < 0.015, dp_e.data > 5, np.isnan(p_e.data), np.isnan(dp_e.data), 
                       distance_array > 300. / 3600.], axis = 0)
mask_850_psnr2 = np.any([p_850.data/dp_850.data < 2, i_850.data / di_850.data < 20, i_850.data < 0.00025, dp_e.data > 5., np.isnan(p_850.data), 
                   np.isnan(dp_850.data), distance_array > 180. / 3600.], axis = 0)
mask_de_psnr2 = np.any([mask_d_psnr2, mask_e_psnr2], axis = 0)
mask_d850_psnr2 = np.any([mask_d_psnr2, mask_850_psnr2], axis = 0)
mask_e850_psnr2 = np.any([mask_850_psnr2, mask_e_psnr2], axis = 0)
mask_all_psnr2 = np.any([mask_d_psnr2, mask_e_psnr2, mask_850_psnr2], axis = 0)

# No selection on P/dP
# (maintaining selection on P = NaN if debiasing)
mask_nosnrP_d = np.any([i_d.data / di_d.data < 20, i_d.data < 0.025, dp_d.data > 5, np.isnan(i_d.data), np.isnan(di_d.data), 
                        np.isnan(p_d.data), np.isnan(dp_d.data), distance_array > 300. / 3600.], axis = 0)
mask_nosnrP_e = np.any([i_e.data / di_e.data < 20, i_e.data < 0.015, dp_e.data > 5, np.isnan(i_e.data), np.isnan(di_e.data),
                        np.isnan(p_e.data), np.isnan(dp_e.data), distance_array > 300. / 3600.], axis = 0)
mask_nosnrP_850 = np.any([i_850.data / di_850.data < 20, i_850.data < 0.00025, dp_850.data > 5, np.isnan(i_850.data), np.isnan(di_850.data),
                          np.isnan(p_850.data), np.isnan(dp_850.data), distance_array > 180. / 3600.], axis = 0)
mask_nosnrP_de = np.any([mask_nosnrP_d, mask_nosnrP_e], axis = 0)
mask_nosnrP_d850 = np.any([mask_nosnrP_d, mask_nosnrP_850], axis = 0)
mask_nosnrP_e850 = np.any([mask_nosnrP_850, mask_nosnrP_e], axis = 0)
mask_nosnrP_all = np.any([mask_nosnrP_d, mask_nosnrP_e, mask_nosnrP_850], axis = 0)
# Additional: only central 90''
mask_nosnrP_850_90as = np.any([mask_nosnrP_850, peak_distance_array > 45. / 3600.], axis = 0)
#mask_nosnrP_850_90as = np.any([i_850.data / di_850.data < 20, i_850.data < 0.00025, dp_850.data > 5, np.isnan(i_850.data), np.isnan(di_850.data),
#                          np.isnan(p_850.data), np.isnan(dp_850.data), peak_distance_array > 45. / 3600.], axis = 0)  
mask_nosnrP_d850_90as = np.any([mask_nosnrP_d, mask_nosnrP_850_90as], axis = 0)
mask_nosnrP_e850_90as = np.any([mask_nosnrP_850_90as, mask_nosnrP_e], axis = 0)
mask_nosnrP_all_90as = np.any([mask_nosnrP_d, mask_nosnrP_e, mask_nosnrP_850_90as], axis = 0)

# No intensity threshold (for advanced quiver plots)
# P / dP > 3
mask_nothr_d = np.any([p_d.data/dp_d.data < 3, i_d.data / di_d.data < 20, np.isnan(p_d.data), np.isnan(dp_d.data), dp_d.data > 5.,
                       distance_array > 300. / 3600.], axis = 0)
mask_nothr_e = np.any([p_e.data/dp_e.data < 3, i_e.data / di_e.data < 20, np.isnan(p_e.data), np.isnan(dp_e.data), dp_e.data > 5.,
                       distance_array > 300. / 3600.], axis = 0)
mask_nothr_850 = np.any([p_850.data/dp_850.data < 3, i_850.data / di_850.data < 20, np.isnan(p_850.data), np.isnan(dp_850.data), dp_850.data > 5.,
                         distance_array > 180. / 3600.], axis = 0)
# P / dP > 2
mask_nothr_d_psnr2 = np.any([p_d.data/dp_d.data < 2, i_d.data / di_d.data < 20, np.isnan(p_d.data), np.isnan(dp_d.data), dp_d.data > 5.,
                             distance_array > 300. / 3600.], axis = 0)
mask_nothr_e_psnr2 = np.any([p_e.data/dp_e.data < 2, i_e.data / di_e.data < 20, np.isnan(p_e.data), np.isnan(dp_e.data), dp_e.data > 5.,
                             distance_array > 300. / 3600.], axis = 0)
mask_nothr_850_psnr2 = np.any([p_850.data/dp_850.data < 2, i_850.data / di_850.data < 20, np.isnan(p_850.data), np.isnan(dp_850.data), dp_850.data > 5.,
                               distance_array > 180. / 3600.], axis = 0)


# Selection on Q/I instead of P
#  (see V01)

## Original format data
# All criteria, P / dP > 3
mask_d_orig = np.any([p_d_orig.data/dp_d_orig.data < 3, i_d_orig.data / di_d_orig.data < 20, i_d_orig.data < 0.025, dp_d_orig.data > 5., 
                      np.isnan(i_d_orig.data), np.isnan(p_d_orig.data), np.isnan(di_d_orig.data), np.isnan(dp_d_orig.data),
                      distance_array_orig_d > 300. / 3600.], axis = 0)
mask_e_orig = np.any([p_e_orig.data/dp_e_orig.data < 3, i_e_orig.data / di_e_orig.data < 20, i_e_orig.data < 0.015, dp_e_orig.data > 5., 
                      np.isnan(i_e_orig.data), np.isnan(p_e_orig.data), np.isnan(di_e_orig.data), np.isnan(dp_e_orig.data), 
                      distance_array_orig_e > 300. / 3600.], axis = 0)
mask_850_orig = np.any([p_850_orig.data/dp_850_orig.data < 3, i_850_orig.data / di_850_orig.data < 20, i_850_orig.data < 0.00025, dp_850_orig.data > 5., 
                        np.isnan(i_850_orig.data), np.isnan(p_850_orig.data), np.isnan(di_850_orig.data), np.isnan(dp_850_orig.data), 
                        distance_array_orig > 180. / 3600.], axis = 0)
# All criteria, P / dP > 2
mask_d_orig_psnr2 = np.any([p_d_orig.data/dp_d_orig.data < 2, i_d_orig.data / di_d_orig.data < 20, i_d_orig.data < 0.025, dp_d_orig.data > 5., 
                            np.isnan(i_d_orig.data), np.isnan(p_d_orig.data), np.isnan(di_d_orig.data), np.isnan(dp_d_orig.data),
                            distance_array_orig_d > 300. / 3600.], axis = 0)
mask_e_orig_psnr2 = np.any([p_e_orig.data/dp_e_orig.data < 2, i_e_orig.data / di_e_orig.data < 20, i_e_orig.data < 0.015, dp_e_orig.data > 5., 
                            np.isnan(i_e_orig.data), np.isnan(p_e_orig.data), np.isnan(di_e_orig.data), np.isnan(dp_e_orig.data),
                            distance_array_orig_e > 300. / 3600.], axis = 0)
mask_850_orig_psnr2 = np.any([p_850_orig.data/dp_850_orig.data < 2, i_850_orig.data / di_850_orig.data < 20, i_850_orig.data < 0.00025, dp_850_orig.data > 5., 
                              np.isnan(i_850_orig.data), np.isnan(p_850_orig.data), np.isnan(di_850_orig.data), np.isnan(dp_850_orig.data),
                              distance_array_orig > 180. / 3600.], axis = 0)
# No selection on P/dP
mask_nosnrP_d_orig = np.any([i_d_orig.data / di_d_orig.data < 20, i_d_orig.data < 0.025, dp_d_orig.data > 5., np.isnan(p_d_orig.data), 
                             np.isnan(dp_d_orig.data), distance_array_orig_d > 300. / 3600.], axis = 0)
mask_nosnrP_e_orig = np.any([i_e_orig.data / di_e_orig.data < 20, i_e_orig.data < 0.015, dp_e_orig.data > 5., np.isnan(p_e_orig.data), 
                             np.isnan(dp_e_orig.data), distance_array_orig_e > 300. / 3600.], axis = 0)
mask_nosnrP_850_orig = np.any([i_850_orig.data / di_850_orig.data < 20, i_850_orig.data < 0.00025, dp_850_orig.data > 5., np.isnan(p_850_orig.data), 
                               np.isnan(dp_850_orig.data), distance_array_orig > 180. / 3600.], axis = 0)
# No selection on I threshold
mask_nothr_d_orig = np.any([p_d_orig.data / dp_d_orig.data < 3, i_d_orig.data / di_d_orig.data < 20, dp_d_orig.data > 5., np.isnan(p_d_orig.data), 
                            np.isnan(dp_d_orig.data), distance_array_orig_d > 300. / 3600.], axis = 0)
mask_nothr_e_orig = np.any([p_e_orig.data / dp_e_orig.data < 3, i_e_orig.data / di_e_orig.data < 20, dp_e_orig.data > 5., np.isnan(p_e_orig.data), 
                            np.isnan(dp_e_orig.data), distance_array_orig_e > 300. / 3600.], axis = 0)
mask_nothr_850_orig = np.any([p_850_orig.data / dp_850_orig.data < 3, i_850_orig.data / di_850_orig.data < 20, dp_850_orig.data > 5., np.isnan(p_850_orig.data), 
                              np.isnan(dp_850_orig.data), distance_array_orig > 180. / 3600.], axis = 0)
# No selection on P/dP, I threshold
mask_nothrnosnrP_d_orig = np.any([i_d_orig.data / di_d_orig.data < 20, dp_d_orig.data > 5., np.isnan(p_d_orig.data), np.isnan(dp_d_orig.data),
                                  distance_array_orig_d > 300. / 3600.], axis = 0)
mask_nothrnosnrP_e_orig = np.any([i_e_orig.data / di_e_orig.data < 20, dp_e_orig.data > 5., np.isnan(p_e_orig.data), np.isnan(dp_e_orig.data),
                                  distance_array_orig_e > 300. / 3600.], axis = 0)
mask_nothrnosnrP_850_orig = np.any([i_850_orig.data / di_850_orig.data < 20, dp_850_orig.data > 5., np.isnan(p_850_orig.data), np.isnan(dp_850_orig.data),
                                    distance_array_orig > 180. / 3600.], axis = 0)

## Resampled but usmoothed
# All criteria, P / dP > 3
mask_d_rs = np.any([p_d_rs.data/dp_d_rs.data < 3, i_d_rs.data / di_d_rs.data < 20, i_d_rs.data < 0.025, dp_d_rs.data > 5., 
                      np.isnan(i_d_rs.data), np.isnan(p_d_rs.data), np.isnan(di_d_rs.data), np.isnan(dp_d_rs.data),
                      distance_array > 300. / 3600.], axis = 0)
mask_e_rs = np.any([p_e_rs.data/dp_e_rs.data < 3, i_e_rs.data / di_e_rs.data < 20, i_e_rs.data < 0.015, dp_e_rs.data > 5., 
                      np.isnan(i_e_rs.data), np.isnan(p_e_rs.data), np.isnan(di_e_rs.data), np.isnan(dp_e_rs.data),
                      distance_array > 300. / 3600.], axis = 0)
mask_850_rs = np.any([p_850_rs.data/dp_850_rs.data < 3, i_850_rs.data / di_850_rs.data < 20, i_850_rs.data < 0.00025, dp_850_rs.data > 5., 
                        np.isnan(i_850_rs.data), np.isnan(p_850_rs.data), np.isnan(di_850_rs.data), np.isnan(dp_850_rs.data), 
                        distance_array > 180. / 3600.], axis = 0)
# No selection on P/dP
mask_nosnrP_d_rs = np.any([i_d_rs.data / di_d_rs.data < 20, i_d_rs.data < 0.025, dp_d_rs.data > 5., np.isnan(p_d_rs.data), 
                             np.isnan(dp_d_rs.data), distance_array > 300. / 3600.], axis = 0)
mask_nosnrP_e_rs = np.any([i_e_rs.data / di_e_rs.data < 20, i_e_rs.data < 0.015, dp_e_rs.data > 5., np.isnan(p_e_rs.data), 
                             np.isnan(dp_e_rs.data), distance_array > 300. / 3600.], axis = 0)
mask_nosnrP_850_rs = np.any([i_850_rs.data / di_850_rs.data < 20, i_850_rs.data < 0.00025, dp_850_rs.data > 5., np.isnan(p_850_rs.data), 
                               np.isnan(dp_850_rs.data), distance_array > 180. / 3600.], axis = 0)
# Other selections to be added as needed


# Masking angle mismatches
# Absolute difference
mask_dtheta_de = np.abs(delta_theta_de.data) > 30.
mask_dtheta_d850 = np.abs(delta_theta_d850.data) > 30.
mask_dtheta_e850 = np.abs(delta_theta_e850.data) > 30.
# Significance
mask_dtheta_SNR_de = np.abs(delta_theta_SNR_de.data) > 3.
mask_dtheta_SNR_d850 = np.abs(delta_theta_SNR_d850.data) > 3.
mask_dtheta_SNR_e850 = np.abs(delta_theta_SNR_e850.data) > 3.

## Masking by temperature
mask_T = T_map.data < 20.
mask_T_origd = T_map_origd.data < 20.
mask_T_orige = T_map_orige.data < 20.

## Masking by optical depth
mask_tau_d = tau_map_d.data > 0.1
mask_tau_d_orig = tau_map_d_orig.data > 0.1
mask_tau_e = tau_map_e.data > 0.1
mask_tau_e_orig = tau_map_e_orig.data > 0.1



### MAPS AND PLOTS ###

# Show standard plot colors:
# plt.rcParams['axes.prop_cycle'].by_key()['color']

## Obs. field or quality selection overlaid on Herschel T / col. density maps

# Preparation
hmap = T_map # nh2_map # 
minval = 12. # 1e21 # 
maxval = 26. # 1e23 # 
barlabel = r'Herschel T$_{\rm d}$ (K)' # r'Herschel column density (H$_{2}$ cm$^{-2}$)' # 
wcs = WCS(i_d.header)

# Actual plot
#plt.rcParams.update({'font.size': 16})  # Default: 10.0
plt.subplot(projection = wcs)
#plt.imshow(hmap.data, cmap = 'cividis', vmin = minval, vmax = maxval)
plt.imshow(hmap.data, cmap = 'viridis', norm = colors.LogNorm(vmin = minval, vmax = maxval))
plt.xlabel('RA (J2000)', fontsize = 'large') # 'xx-large') # 
plt.ylabel('Dec (J2000)', fontsize = 'large') # 'xx-large') # 
cbar = plt.colorbar()
cbar.set_label(barlabel, size = 'large') # 'xx-large') # 
#cbar.ax.tick_params(labelsize = 'x-large')
#plt.contour(mask_tau_d, colors = 'black', levels=[.5])
plt.contour(distance_array, colors = 'white', levels=[180./3600., 300./3600.], linestyles = ':', linewidths = 2)
IRS1_coord = wcs.world_to_pixel_values(86.7698, 0.36224)
plt.scatter(IRS1_coord[0], IRS1_coord[1], facecolor = 'white', edgecolor = 'black')  # IRS 1
#plt.savefig('Tmap_plus_distcircles')#('coldensmap_plus_distcircles')#
#plt.rcParams.update({'font.size': 10})   # ???
#mpl.rcParams.update(mpl.rcParamsDefault) # ???


## Contour levels

# I plus I contours
# 1. Preparation
I2plot = copy.deepcopy(i_e)
wcs = WCS(I2plot)
Imin = -.01 # -.02 # -.0005 # 
Imax = .25 # .5 # .01 # 
lvls = [.008, .015, .03, .06] # [.01, .025, .05, .1] # [.0001, .00025, .0005, .001] # 
labels = [r'8 mJy/arcsec$^2$', '15', '30', '60'] 
title = r'Band E (214 $\mu$m) intensity'
# 2. Plot
plt.subplot(projection=wcs)
plt.imshow(I2plot.data, cmap = 'rainbow', vmin = Imin, vmax = Imax)
plt.title(title)
plt.xlabel('RA')
plt.ylabel('Dec')
cbar = plt.colorbar()
cbar.set_label('I (Jy / arcsec2)')
contours = plt.contour(I2plot.data, colors = ['black'], levels = lvls, linewidths = [.5, .75, 1., 1.5])
#contours = plt.contour(I2plot.data, colors = ['white', 'lightgrey', 'black', 'black'], linestyles = ['-', '-', '--', '-'], 
#                       levels = lvls, linewidths = 1.5)
#plt.clabel(contours, inline = True, fontsize = 8)
for i in range(len(labels)):
    contours.collections[i].set_label(labels[i])
plt.legend(loc='upper right')
#plt.savefig('N2071_I_154um_contourlevels')

# T, P map plus selection
wcs = WCS(i_d.header)
plt.subplot(projection=wcs)
plt.imshow(T_map.data, cmap = 'rainbow')
plt.title('Dust T (Herschel)')
#plt.imshow(p_d.data, cmap = 'rainbow', vmin = 0., vmax = 25.)
#plt.title('Pol. frac. (debiased) at 154 um')
plt.xlabel('RA')
plt.ylabel('Dec')
cbar = plt.colorbar()
cbar.set_label(r'T_d (Herschel)')
#cbar.set_label('P (%), debiased')
plt.contour(mask_d850.astype(int), colors='black', levels=[.5], smooth=1, kernel='box', linewidths=1.5)       # Quality selection
#plt.contour(i_d.data, colors='white', levels=[.1], smooth=1, kernel='box', linewidths=2)                    # Turning point in P-I plt
#plt.contour(tau_map_d.data, colors='black', levels=[.05, .1, .2], smooth=1, kernel='box', linewidths=1.5)   # Tau 
#plt.savefig('Tmap_plus_d850sel_sel02')


'''
# Plots of polarized quantities
#  Options: Central 3' vs. full field
#           Smoothed+regridded vs. original
# Setup
id2plot = copy.deepcopy(i_d)#(i_d_orig)#
qd2plot = copy.deepcopy(q_d)#(q_d_orig)#
ud2plot = copy.deepcopy(u_d)#(u_d_orig)#
ie2plot = copy.deepcopy(i_e)#(i_e_orig)#
qe2plot = copy.deepcopy(q_e)#(q_e_orig)#
ue2plot = copy.deepcopy(u_e)#(u_e_orig)#
i8502plot = copy.deepcopy(i_850)#(i_850_orig)#
q8502plot = copy.deepcopy(q_850)#(q_850_orig)#
u8502plot = copy.deepcopy(u_850)#(u_850_orig)#
Ipd2plot = copy.deepcopy(Ip_d)#(Ip_d_orig)#
Ipe2plot = copy.deepcopy(Ip_e)#(Ip_e_orig)#
Ip8502plot = copy.deepcopy(Ip_850)#(Ip_850_orig)#
pd2plot = copy.deepcopy(p_d)#(p_d_orig)#
pe2plot = copy.deepcopy(p_e)#(p_e_orig)#
p8502plot = copy.deepcopy(p_850)#(p_850_orig)#
thetad2plot = copy.deepcopy(theta_d)#(theta_d_orig)#
thetae2plot = copy.deepcopy(theta_e)#(theta_e_orig)#
theta8502plot = copy.deepcopy(theta_850)#(theta_850_orig)#
zoomcentral3arcm = False
plotcontours = True
'''

# I, Q, U maps
#  (see V01)

# P, I, PI in central 3' ()
#  (see V01)


## Advanced quiver plot
i2plot = copy.deepcopy(i_d)
p2plot = copy.deepcopy(p_d)
theta2plot = theta_d
mask_std = mask_d
mask_snr2 = mask_d_psnr2
ormask_dtheta = np.any([mask_std, mask_dtheta_SNR_de, mask_dtheta_SNR_d850, mask_dtheta_SNR_e850], axis = 0)
andmask_dtheta = mask_std * mask_dtheta_SNR_de * mask_dtheta_SNR_d850 * mask_dtheta_SNR_e850
barlabel = r'I (154 $\mu$m) in Jy/arcsec$^2$' # r'I (214 $\mu$m) in Jy/arcsec$^2$' # r'I (850 $\mu$m) in Jy/arcsec$^2$' #
if pixsize == '4.0':
    stepvec = 4
    scalevec = 1.2     # (1pix = scalevec * 1% pol)
    pxscale2plot = 4.  # Pixel size in arcsec
elif pixsize == '8.0':
    stepvec = 2
    scalevec = .6
    pxscale2plot = 8.
else:
    pass
fig1 = FITSFigure(i2plot, subplot=(1,1,1))
fig1.axis_labels.set_font(size='xx-large')
fig1.tick_labels.set_font(size='xx-large')
fig1.show_colorscale(cmap='rainbow', vmin = -2.5e-2, vmax = 2.5e-1) # vmin = -2.5e-4, vmax = 2.5e-3) # 
fig1.show_markers(86.7698, 0.36224, marker = "o", facecolor = 'white', edgecolor = 'black')  # IRS 1 
fig1.add_colorbar() # (format = '%.0e')
fig1.colorbar.set_font(size='xx-large')
fig1.colorbar.set_axis_label_font(size='xx-large')
fig1.colorbar.set_axis_label_text(barlabel)
fig1.show_contour(mask_std.astype(int), levels = [.5], linewidths = .75)
fig1.show_contour(ormask_dtheta.astype(int), levels = [.5], linewidths = .75)
fig1.show_contour(distance_array, linestyles = ':', linewidths = 1.5, levels = [300./3600.]) # , levels = [180./3600.]) #
'''
# If you want scaled vectors
vectscale = scalevec * pxscale2plot/3600
p2plot.data[mask_snr2] = np.nan
fig1.show_vectors(p2plot, theta2plot, scale=scalevec, step = stepvec, color = 'white')
p2plot.data[mask_std] = np.nan
fig1.show_vectors(p2plot, theta2plot, scale=scalevec, step = stepvec, color = 'black')
fig1.add_scalebar(20 * vectscale, "p = 20%",corner='top right',frame=True)
fig1.scalebar.set_font(size='xx-large')
fig1.scalebar.set_linewidth(1.5)
'''
# If you want fixed-length vectors
p2plot.data /= p2plot.data
p2plot.data *= scalevec
p2plot.data[mask_snr2] = np.nan
fig1.show_vectors(p2plot, theta2plot, step = stepvec, scale = 3, color = 'white', linewidth = 2.)
p2plot.data[mask_std] = np.nan
fig1.show_vectors(p2plot, theta2plot, step = stepvec, scale = 3, color = 'black', linewidth = 2.)
# Zoom on central 3'
x_center = i2plot.header['CRVAL1']
y_center = i2plot.header['CRVAL2']
fig1.recenter(x_center, y_center, radius=90./3600.)
fig1.show_contour(mask_std.astype(int), levels = [.5], linewidths = 1., color = 'black')
fig1.show_markers(86.7723, 0.36390, marker = "o", facecolor = 'white', edgecolor = 'black')  # IRS 2
fig1.show_markers(86.7698, 0.36391, marker = "o", facecolor = 'white', edgecolor = 'black')  # IRS 3
fig1.add_label(86.769, 0.36224, '1', horizontalalignment = 'left', verticalalignment = 'center', weight = 'bold', size = 'xx-large')
fig1.add_label(86.773, 0.36390, '2', horizontalalignment = 'right', verticalalignment = 'center', weight = 'bold', size = 'xx-large')
fig1.add_label(86.769, 0.36391, '3', horizontalalignment = 'left', verticalalignment = 'center', weight = 'bold', size = 'xx-large')
#plt.savefig('Quiverplot_I+pol154_newpipeline_pix4.0_3arcmin.png')

# Center coord (8 arcsec pixel):
#  Pixel (66, 66)
#  Deg (86.771, 0.364361111111108)
# Max value (8 arcsec pixel):
#  154 um: 1.88
#  214 um: ...
#  850 um: ...


## Quiver plot comparison
# Updated version
frankenstein = fits.PrimaryHDU()
frankenstein.data = nh2_map.data
frankenstein.header = p_d.header
fig1 = FITSFigure(frankenstein, subplot=(1,1,1))
fig1.show_colorscale(cmap = 'viridis', vmin = 0, vmax = 1e23)
fig1.show_contour(peak_distance_array, colors = 'black', levels = [45./3600.], linestyles = '--')
#fig1.show_contour(distance_array, colors = 'black', levels = [180./3600.], linestyles = ':')
fig1.add_colorbar()
fig1.colorbar.set_axis_label_text(r'Herschel column density (H$_{2}$ cm$^{-2}$)')
fig1.axis_labels.set_font(size='xx-large')
fig1.tick_labels.set_font(size='xx-large')
fig1.colorbar.set_font(size='xx-large')
fig1.colorbar.set_axis_label_font(size='xx-large')
x_center = frankenstein.header['CRVAL1']
y_center = frankenstein.header['CRVAL2']
fig1.recenter(x_center, y_center, radius = 180./3600.)
p_d_2plot = copy.deepcopy(p_d)
p_d_2plot.data[mask_d] = np.nan
p_e_2plot = copy.deepcopy(p_e)
p_e_2plot.data[mask_e] = np.nan
p_850_2plot = copy.deepcopy(p_850)
p_850_2plot.data[mask_850] = np.nan
p_850_2plot.data /= p_850_2plot.data
p_d_2plot.data /= p_d_2plot.data
p_e_2plot.data /= p_e_2plot.data
fig1.show_vectors(p_d_2plot, theta_d, step = 2, scale = 2, color = 'white')
fig1.show_vectors(p_e_2plot, theta_e, step = 2, scale = 2, color = 'magenta')
fig1.show_vectors(p_850_2plot, theta_850, step = 2, scale = 2, color = 'black')
fig1.show_markers(86.7698, 0.36224, marker = "o", facecolor = 'white', edgecolor = 'black')  # IRS 1 
fig1.show_markers(86.7723, 0.36390, marker = "o", facecolor = 'white', edgecolor = 'black')  # IRS 2
fig1.show_markers(86.7698, 0.36391, marker = "o", facecolor = 'white', edgecolor = 'black')  # IRS 3
fig1.add_label(86.769, 0.36224, '1', horizontalalignment = 'left', verticalalignment = 'center', weight = 'bold', size = 'xx-large')
fig1.add_label(86.773, 0.36390, '2', horizontalalignment = 'right', verticalalignment = 'center', weight = 'bold', size = 'xx-large')
fig1.add_label(86.769, 0.36391, '3', horizontalalignment = 'left', verticalalignment = 'center', weight = 'bold', size = 'xx-large')
#plt.savefig('Quiverplot_compare_3bands_pix8_newpipeline_90as')


# Q, U uncertainties (need them ~uniform for Rice)
wcs = WCS(i_d_orig.header) # (i_d.header) # 
plt.subplot(projection=wcs)
plt.imshow(dq_d_orig.data * 1e3, cmap = 'rainbow',  norm=colors.LogNorm(vmin = 1e-1, vmax =1.)) # (dq_d_orig.data * 1e3, cmap = 'rainbow', vmin = 0, vmax = 1e-1) # 
plt.title('Q uncertainty map, D band') # ('Q uncertainty map (D band, smoothed)') # 
plt.xlabel('RA')
plt.ylabel('Dec')
cbar = plt.colorbar()
cbar.set_label(r'$\delta Q$ in mJy / arcsec$^2$')
# (mask_d.astype(int), colors='black', levels=[.5], linewidths=1.5) # 
plt.contour(mask_central6_orig.astype(int), colors='black', levels=[.5], linewidths=1., linestyles = 'dashed') # (mask_central6.astype(int), ... # 
#plt.savefig('N2071_Dband_dQ.png')



## Angle differences

# Parameters
whatdata = 'de'
central90as = True
absval = False

# Setup
if whatdata == 'de':
    bartext = r'$\Delta {\rm \theta}$ (154 - 214 $\mu$m)'
    selection = mask_de
    delta_theta = np.ma.masked_where(selection, delta_theta_de)
    delta_theta_SNR = np.ma.masked_where(selection, delta_theta_SNR_de)
elif whatdata == 'd850':
    bartext = r'$\Delta {\rm \theta}$ (154 - 850 $\mu$m)'
    if central90as == True:
        selection = mask_d850_90as
    else:
        selection = mask_d850
    delta_theta = np.ma.masked_where(selection, delta_theta_d850)
    delta_theta_SNR = np.ma.masked_where(selection, delta_theta_SNR_d850)
elif whatdata == 'e850':
    bartext = r'$\Delta {\rm \theta}$ (214 - 850 $\mu$m)'
    if central90as == True:
        selection = mask_e850_90as
    else:
        selection = mask_e850
    delta_theta = np.ma.masked_where(selection, delta_theta_e850)
    delta_theta_SNR = np.ma.masked_where(selection, delta_theta_SNR_e850)
else:
    pass

# Plot (central zoom version w/ FITSFigure))
delta_theta_4aplpy = fits.PrimaryHDU()                        # Need to create an HDU
if absval == False:
    delta_theta_4aplpy.data = np.ma.filled(delta_theta, np.nan)
    cmap2use = 'twilight_shifted'
    mintheta = -90.
    maxtheta = 90.
else:
    delta_theta_4aplpy.data = np.ma.filled(np.abs(delta_theta), np.nan)
    cmap2use = 'cividis_r'
    mintheta = 0.
    maxtheta = 90.
delta_theta_4aplpy.header = angles_d.header
ax_zoom = FITSFigure(delta_theta_4aplpy, subplot = (1,1,1))
ax_zoom.show_colorscale(cmap = cmap2use, vmin = mintheta, vmax = maxtheta)
#if absval == False:
#    ax_zoom.show_colorscale(cmap='twilight_shifted', vmin=-90., vmax=90.)
#else:
#    ax_zoom.show_colorscale(cmap='cividis_r', vmin=-90., vmax=90.)
ax_zoom.add_colorbar()
ax_zoom.colorbar.set_axis_label_text(bartext)
ax_zoom.axis_labels.set_font(size='xx-large')
ax_zoom.tick_labels.set_font(size='xx-large')
ax_zoom.colorbar.set_font(size='xx-large')
ax_zoom.colorbar.set_axis_label_font(size='xx-large')
x_center = angles_d.header['CRVAL1']  # Image center coordinates in deg
y_center = angles_d.header['CRVAL2']
ax_zoom.recenter(x_center, y_center, radius=90./3600.)
#contours1 = ax_zoom.show_contour(delta_theta_SNR.data, colors=['black'], levels=[5.], linewidths=[2.])
#contours2 = ax_zoom.show_contour(nh2_map.data, colors=['black'], levels=[5e22], linewidths=[3.], linestyles = ['dashed'])
ax_zoom.show_contour('Ancillary_Herschel_data/N2071_NH2_cutout-15arcm.fits', levels = [5e22], colors = 'black', linewidths = 3., linestyles = 'dashed')
ax_zoom.show_contour('CO_maps/NGC2071_12CO_hred_integ_r.fits', levels = [2.1], colors = 'black') # [2.1, 7., 21.] # 
ax_zoom.show_contour('CO_maps/NGC2071_12CO_hblue_integ_r.fits', levels = [2.1], colors = 'black')
ax_zoom.show_contour('CO_maps/NGC2071_C18O_blue_integ_r.fits', levels = [2.5], colors = 'white', linewidths = 3.) # [2.5, 2.9, 3.3, 3.7] # 
ax_zoom.show_contour('CO_maps/NGC2071_C18O_blue_integ_r.fits', levels = [2.5], colors = 'blue', linewidths = 1.)
ax_zoom.show_contour('CO_maps/NGC2071_C18O_red_integ_r.fits', levels = [2.5], colors = 'white', linewidths = 3.)
ax_zoom.show_contour('CO_maps/NGC2071_C18O_red_integ_r.fits', levels = [2.5], colors = 'red', linewidths = 1.)
ax_zoom.show_markers(86.7698, 0.36224, marker = "o", facecolor = 'white', edgecolor = 'black')  # IRS 1 
ax_zoom.show_markers(86.7723, 0.36390, marker = "o", facecolor = 'white', edgecolor = 'black')  # IRS 2
ax_zoom.show_markers(86.7698, 0.36391, marker = "o", facecolor = 'white', edgecolor = 'black')  # IRS 3
ax_zoom.add_label(86.769, 0.36224, '1', horizontalalignment = 'left', verticalalignment = 'center', weight = 'bold', size='xx-large')
ax_zoom.add_label(86.773, 0.36390, '2', horizontalalignment = 'right', verticalalignment = 'center', weight = 'bold', size='xx-large')
ax_zoom.add_label(86.769, 0.36391, '3', horizontalalignment = 'left', verticalalignment = 'center', weight = 'bold', size='xx-large')
#ax_zoom.savefig('N2071_DTheta_154-850_zoom_newpipeline_pix8.png')

'''
# Correlation between 154-214 and 154-850 angle differences
x_dangle = np.ma.masked_where(mask_all, delta_theta_de)
y_dangle = np.ma.masked_where(mask_all, delta_theta_d850)
fig, ax = plt.subplots(1, 1)
plt.scatter(x_dangle, y_dangle)
plt.xlabel('154 - 214 um angle difference')
plt.ylabel('154 - 850 um angle difference')
#plt.axis(axlims)
#plt.savefig('DTheta_comparison_noerrorbars.png')
'''

# Location of high-DTheta points on P map
# New version, new pipeline

# Data preparation
p2plot = fits.PrimaryHDU()
whatdata = 'd850' # 'e850' # 'de' # 
if whatdata == 'de':
    figtitle = r'154-214 $\mu$m pol. angle difference' 
    savename = 'High-Dtheta154-214-locations_3sigma.png' 
    dummy = np.ma.masked_where(mask_de, delta_theta_de.data) 
    sel_angle = np.any([mask_de, delta_theta_SNR_de < 3.], axis = 0) 
    sel_angle2 = np.any([mask_de, delta_theta_SNR_de < 5.], axis = 0)
    #sel_angle = np.any([mask_de, np.abs(delta_theta_de.data) < 30.], axis = 0) # (Old version)
    plot6arcmin = False
elif whatdata == 'd850':
    figtitle = r'154-850 $\mu$m pol. angle difference' 
    savename = 'High-Dtheta154-850-locations_3sigma.png' 
    dummy = np.ma.masked_where(mask_d850, delta_theta_d850.data) 
    sel_angle = np.any([mask_d850, delta_theta_SNR_d850 < 3.], axis = 0) 
    sel_angle2 = np.any([mask_d850, delta_theta_SNR_d850 < 5.], axis = 0)
    plot6arcmin = True
elif whatdata == 'e850':
    figtitle = r'214-850 $\mu$m pol. angle difference' 
    savename = 'High-Dtheta214-850-locations_3sigma.png' 
    dummy = np.ma.masked_where(mask_e850, delta_theta_e850.data) 
    sel_angle = np.any([mask_e850, delta_theta_SNR_e850 < 3.], axis = 0) 
    sel_angle2 = np.any([mask_e850, delta_theta_SNR_e850 < 5.], axis = 0)
    plot6arcmin = True
else:
    pass
p2plot.data = np.ma.filled(dummy, np.nan)
p2plot.header = i_d.header
x_center = p2plot.header['CRVAL1']
y_center = p2plot.header['CRVAL2']

fig = FITSFigure(p2plot, subplot = (1,1,1))
fig.show_colorscale(cmap='twilight_shifted', vmin = -90., vmax = 90.)
fig.add_colorbar()
fig.recenter(x_center, y_center, radius = 300./3600.)
fig.show_contour(sel_angle.astype(int), levels=[.5], colors='black', linewidths = 2.) # , linestyles = ':') 
#fig.show_contour(sel_angle2.astype(int), levels=[.5], colors='black', linewidths = 2.)
if plot6arcmin == True:
    fig.show_contour(peak_distance_array, levels=[45./3600.], colors='black', linewidths = 2., linestyles = ['--'])
    fig.show_contour(distance_array, levels=[180./3600.], colors='black', linewidths = 2., linestyles = [':'])
fig.set_title(figtitle, fontsize = 'xx-large') 
fig.axis_labels.set_font(size='xx-large')
fig.tick_labels.set_font(size='xx-large')
fig.axis_labels.set_font(size='xx-large')
fig.tick_labels.set_font(size='xx-large')
fig.colorbar.set_font(size = 'x-large')
#plt.savefig(savename)

# Old version, old pipeline:
#  (See V01)
#   (includes older version of 214-850 plot)


## SNR GAIN
wcsraw = WCS(i_d_orig.header)
wcs = WCS(i_d.header)
SNR_i_d_orig = i_d_orig.data/di_d_orig.data
SNR_i_d = i_d.data / di_d.data
SNR_i_e = i_e.data / di_e.data

plt.subplot(projection = wcsraw)
plt.imshow(i_d_orig.data, cmap = 'rainbow', vmin=-.03, vmax=.5)
plt.colorbar()
cts = plt.contour(SNR_i_d_orig, colors='black', levels=[10., 100., 1000., 10000.], smooth=1, kernel='box', linewidths=1.5)
plt.clabel(cts, inline=True, fmt='%1.1f')
plt.title('I (154 um) before smothing + S/N contours')
plt.xlabel('RA')
plt.ylabel('Dec')

plt.subplot(projection = wcs)
plt.imshow(i_d.data, cmap = 'rainbow', vmin=-.03, vmax=.5)
plt.colorbar()
cts = plt.contour(SNR_i_d, colors='black', levels=[100., 1000., 10000.], smooth=1, kernel='box', linewidths=1.5)
plt.clabel(cts, inline=True, fmt='%1.1f')
plt.xlabel('RA')
plt.ylabel('Dec')

plt.subplot(projection = wcs)
plt.imshow(i_e.data, cmap = 'rainbow', vmin=-.03, vmax=.5)
plt.colorbar()
cts = plt.contour(SNR_i_e, colors='black', levels=[10., 100., 1000., 10000.], smooth=1, kernel='box', linewidths=1.5)
plt.clabel(cts, inline=True, fmt='%1.1f')
plt.xlabel('RA')
plt.ylabel('Dec')




### SCATTER PLOTS ###

## P vs. I plots

lowess = sm.nonparametric.lowess
i_d_orig_dmasked = np.ma.masked_where(mask_nosnrP_d_orig, i_d_orig.data)
p_d_orig_dmasked = np.ma.masked_where(mask_nosnrP_d_orig, p_d_orig.data)
i_d_orig_opthin = np.ma.masked_where(np.any([mask_nosnrP_d_orig, mask_tau_d_orig], axis = 0), i_d_orig.data)    # One-band selection
p_d_orig_opthin = np.ma.masked_where(np.any([mask_nosnrP_d_orig, mask_tau_d_orig], axis = 0), p_d_orig.data)
i_e_orig_emasked = np.ma.masked_where(mask_nosnrP_e_orig, i_e_orig.data)
p_e_orig_emasked = np.ma.masked_where(mask_nosnrP_e_orig, p_e_orig.data)
i_e_orig_opthin = np.ma.masked_where(np.any([mask_nosnrP_e_orig, mask_tau_e_orig], axis = 0), i_e_orig.data)
p_e_orig_opthin = np.ma.masked_where(np.any([mask_nosnrP_e_orig, mask_tau_e_orig], axis = 0), p_e_orig.data)
if pixsize == '4.0':
    i_850_orig_850masked = np.ma.masked_where(mask_nosnrP_850_orig, i_850_orig.data)
    p_850_orig_850masked = np.ma.masked_where(mask_nosnrP_850_orig, p_850_orig.data)
    # (Nota: no need for 'opthin' version)
elif pixsize == '8.0':
    i_850_orig_850masked = np.ma.masked_where(mask_nosnrP_850_rs, i_850_rs.data)
    p_850_orig_850masked = np.ma.masked_where(mask_nosnrP_850_rs, p_850_rs.data)
else:
    pass

fig, ax = plt.subplots(1) # 154 um
plt.xscale('log')         # Using 2D histogram due to the large number of data points
plt.yscale('log')
xbins = np.logspace(np.log10(.025), np.log10(3.), 30)
ybins = np.logspace(-1., 2., 30)
h = plt.hist2d(i_d_orig_dmasked.ravel().compressed(), p_d_orig_dmasked.ravel().compressed(), bins = [xbins, ybins], cmap = 'Blues', norm = colors.LogNorm(vmin = 1., vmax = 100.))
fig.colorbar(h[3], ax=ax)
# LOWESS
data_lowess = lowess(p_d_orig_dmasked.ravel().compressed(), i_d_orig_dmasked.ravel().compressed(), frac = .5)  
i_lowess = data_lowess[:, 0]
p_lowess = data_lowess[:, 1]
data_lowess_opthin = lowess(p_d_orig_opthin.ravel().compressed(), i_d_orig_opthin.ravel().compressed(), frac = .5)
i_lowess_opthin = data_lowess_opthin[:, 0]
p_lowess_opthin = data_lowess_opthin[:, 1]
plt.plot(i_lowess, p_lowess, color = 'black')
plt.plot(i_lowess_opthin, p_lowess_opthin, color = 'black', linestyle = '--')
# Error bar example
dp_d_orig_dmasked = np.ma.masked_where(mask_nosnrP_d_orig, dp_d_orig.data)  # mask_d_nosnrP_orig ??
rel_error = np.nanmedian(dp_d_orig_dmasked.ravel().compressed() / p_d_orig_dmasked.ravel().compressed())
plt.errorbar(2., 4., yerr = 4. * rel_error, capsize = 3)
# Final touches
plt.axis([1e-2, 3, 1e-1, 1e2])
plt.text(np.sqrt(3)/10, 70., r'HAWC+ band D (154 $\mu$m)', backgroundcolor = 'white', horizontalalignment = 'center')
plt.ylabel('$P$ (%)')
plt.xlabel(r'$I$ (Jy arcsec$^{-2}$)')
plt.legend([r'LOWESS (all $\tau$)', r'LOWESS ($\tau \leq 0.1)$'], loc = (.55, .78)) 
#plt.savefig('P-vs-I-154um_newpipeline_lowess+errbar+norm')

fig, ax = plt.subplots(1)  # 214 um
plt.xscale('log')
plt.yscale('log')
xbins = np.logspace(np.log10(.015), np.log10(3.), 30)
ybins = np.logspace(-1., 2., 30)
h = plt.hist2d(i_e_orig_emasked.ravel().compressed(), p_e_orig_emasked.ravel().compressed(), bins = [xbins, ybins], cmap = 'Blues', norm = colors.LogNorm(vmin = 1., vmax = 100.))
plt.colorbar(h[3], ax=ax)
# LOWESS
data_lowess = lowess(p_e_orig_emasked.ravel().compressed(), i_e_orig_emasked.ravel().compressed(), frac = .5)  
i_lowess = data_lowess[:, 0]
p_lowess = data_lowess[:, 1]
data_lowess_opthin = lowess(p_e_orig_opthin.ravel().compressed(), i_e_orig_opthin.ravel().compressed(), frac = .5)
i_lowess_opthin = data_lowess_opthin[:, 0]
p_lowess_opthin = data_lowess_opthin[:, 1]
plt.plot(i_lowess, p_lowess, color = 'black')
plt.plot(i_lowess_opthin, p_lowess_opthin, color = 'black', linestyle = '--')
# Error bar example
dp_e_orig_emasked = np.ma.masked_where(mask_nosnrP_e_orig, dp_e_orig.data)  # mask_e_nosnrP_orig
rel_error = np.nanmedian(dp_e_orig_emasked.ravel().compressed() / p_e_orig_emasked.ravel().compressed())
plt.errorbar(2., 4., yerr = 4. * rel_error, capsize = 3)
# Finishing plot
plt.text(np.sqrt(3)/10, 70., r'HAWC+ band E (214 $\mu$m)', backgroundcolor = 'white', horizontalalignment = 'center')
plt.axis([1e-2, 3, 1e-1, 1e2])
plt.ylabel('$P$ (%)')
plt.xlabel(r'$I$ (Jy arcsec$^{-2}$)')
#plt.savefig('P-vs-I-214um_newpipeline_lowess+errbar+norm')

fig, ax = plt.subplots(1)  # 850 um
plt.xscale('log')
plt.yscale('log')
xbins = np.logspace(np.log10(2.5)-4, np.log10(3.)-2, 30)
ybins = np.logspace(-1., 2., 30)
h = plt.hist2d(i_850_orig_850masked.ravel().compressed(), p_850_orig_850masked.ravel().compressed(), bins = [xbins, ybins], cmap = 'Blues', norm = colors.LogNorm(vmin = .5, vmax = 20.))
plt.colorbar(h[3], ax=ax)
# LOWESS
data_lowess = lowess(p_850_orig_850masked.ravel().compressed(), i_850_orig_850masked.ravel().compressed(), frac = .5)  
i_lowess = data_lowess[:, 0]
p_lowess = data_lowess[:, 1]
plt.plot(i_lowess, p_lowess, color = 'black')
# Error bar
if pixsize == '4.0':
    dp_850_orig_850masked = np.ma.masked_where(mask_nosnrP_850_orig, dp_850_orig.data)
elif pixsize == '8.0':
    dp_850_orig_850masked = np.ma.masked_where(mask_nosnrP_850_rs, dp_850_rs.data)
else:
    pass
rel_error = np.nanmedian(dp_850_orig_850masked.ravel().compressed() / p_850_orig_850masked.ravel().compressed())
plt.errorbar(2e-2, 4., yerr = 4. * rel_error, capsize = 3)
# Finish
plt.text(np.sqrt(3)*1e-3, 70., r'POL-2 (850 $\mu$m)', backgroundcolor = 'white', horizontalalignment = 'center')
plt.axis([1e-4, 3e-2, 1e-1, 1e2])
plt.ylabel('$P$ (%)')
plt.xlabel(r'$I$ (Jy arcsec$^{-2}$)')
#plt.savefig('P-vs-I-850um_lowess+errbar+norm')


## VARIANT: tau instead of I
tau_d_orig_dmasked = np.ma.masked_where(mask_nosnrP_d_orig, tau_map_d_orig.data)
tau_e_orig_emasked = np.ma.masked_where(mask_nosnrP_e_orig, tau_map_e_orig.data)
tau_map_850 = fits.PrimaryHDU()
tau_map_850.data = copy.deepcopy(tau_map_d.data) * (154./850.)**2
tau_map_850.header = tau_map_d.header
if pixsize == '4.0':
    tau_850_orig_850masked = np.ma.masked_where(mask_nosnrP_850_orig, tau_map_850.data)
elif pixsize == '8.0':
    tau_850_orig_850masked = np.ma.masked_where(mask_nosnrP_850_rs, tau_map_850.data)
else:
    pass

fig, ax = plt.subplots(1)
plt.xscale('log')
plt.yscale('log')
#plt.scatter(tau_d_orig_dmasked, p_d_orig_dmasked)
xbins = np.logspace(np.log10(5)-3, 0., 50)
ybins = np.logspace(-1., 2., 30)
h = plt.hist2d(tau_d_orig_dmasked.ravel().compressed(), p_d_orig_dmasked.ravel().compressed(), bins = [xbins, ybins], cmap = 'Blues', norm = colors.LogNorm(vmin = 1., vmax = 100.))
plt.colorbar(h[3], ax=ax)
data_lowess = lowess(p_d_orig_dmasked.ravel().compressed(), tau_d_orig_dmasked.ravel().compressed(), frac = .5)
tau_lowess = data_lowess[:, 0]
p_lowess = data_lowess[:, 1]
plt.plot(tau_lowess, p_lowess, color = 'black')
plt.axis([5e-3, 1, 1e-1, 1e2])
plt.text(np.sqrt(5e-3), 70., r'HAWC+ band D (154 $\mu$m)', backgroundcolor = 'white', horizontalalignment = 'center')
plt.ylabel('$P$ (%)')
plt.xlabel(r'$\tau$ from HGBS maps') # (r'$\tau$ at 154 $\mu$m') # 
plt.plot([.1, .1], [1e-1, 1e2], color = 'black', linestyle = ':')
#plt.savefig('P-vs-tau-154um_newpipeline_lowess+norm')

fig, ax = plt.subplots(1)
#plt.scatter(tau_e_orig_emasked, p_e_orig_emasked)
xbins = np.logspace(np.log10(5)-3, 0., 50)
ybins = np.logspace(-1., 2., 30)
h = plt.hist2d(tau_e_orig_emasked.ravel().compressed(), p_e_orig_emasked.ravel().compressed(), bins = [xbins, ybins], cmap = 'Blues', norm = colors.LogNorm(vmin = 1., vmax = 100.))
plt.colorbar(h[3], ax=ax)
data_lowess = lowess(p_e_orig_emasked.ravel().compressed(), tau_e_orig_emasked.ravel().compressed(), frac = .5)
tau_lowess = data_lowess[:, 0]
p_lowess = data_lowess[:, 1]
plt.plot(tau_lowess, p_lowess, color = 'black')
plt.xscale('log')
plt.yscale('log')
plt.axis([5e-3, 1, 1e-1, 1e2])
plt.text(np.sqrt(5e-3), 70., r'HAWC+ band E (214 $\mu$m)', backgroundcolor = 'white', horizontalalignment = 'center')
plt.ylabel('$P$ (%)')
plt.xlabel(r'$\tau$ from HGBS maps')
plt.plot([.1, .1], [1e-1, 1e2], color = 'black', linestyle = ':')
#plt.savefig('P-vs-tau-214um_newpipeline_lowess+norm')

fig, ax = plt.subplots(1)
#plt.scatter(tau_850_orig_850masked, p_850_orig_850masked)
xbins = np.logspace(np.log10(2)-4, np.log10(4)-2, 50)
ybins = np.logspace(-1., 2., 30)
h = plt.hist2d(tau_850_orig_850masked.ravel().compressed(), p_850_orig_850masked.ravel().compressed(), bins = [xbins, ybins], cmap = 'Blues', norm = colors.LogNorm(vmin = .5, vmax = 20.))
plt.colorbar(h[3], ax=ax)
data_lowess = lowess(p_850_orig_850masked.ravel().compressed(), tau_850_orig_850masked.ravel().compressed(), frac = .5)
tau_lowess = data_lowess[:, 0]
p_lowess = data_lowess[:, 1]
plt.plot(tau_lowess, p_lowess, color = 'black')
plt.xscale('log')
plt.yscale('log')
plt.axis([2e-4, 4e-2, 1e-1, 1e2])
plt.text(np.sqrt(8e-6), 70., r'POL-2 (850 $\mu$m)', backgroundcolor = 'white', horizontalalignment = 'center')
plt.ylabel('$P$ (%)')
plt.xlabel(r'$\tau$ from HGBS maps')
#plt.savefig('P-vs-tau-850um_lowess+norm')

'''
# Save data as ASCII
i_850_2save = i_850_orig_850masked.ravel().compressed()
di_850_2save = di_850_orig_850masked.ravel().compressed()
p_850_2save = p_850_orig_850masked.ravel().compressed()
dp_850_2save = dp_850_orig_850masked.ravel().compressed()
tau_850_2save = tau_850_orig_850masked.ravel().compressed()
save_850 = np.array([i_850_2save, di_850_2save, p_850_2save, dp_850_2save, tau_850_2save])
f = open('PItau_850um.asc', 'w')
np.savetxt(f, save_850.T)
f.close
'''



## P vs. P plots

LegendElements = [Line2D([0], [0], marker='o', linestyle="", label="P/$\delta$P > 3"),
                  Line2D([0], [0], marker='$\u25EF$', linestyle="", label="P/$\delta$P < 3")]
                  #Line2D([0], [0], marker='$\u2298$', linestyle="", label="2 < P/$\delta$P < 3"), 
                  #Line2D([0], [0], marker='$\u25EF$', linestyle="", label="P/$\delta$P < 2")]
LegendElementsLowess = [Line2D([0], [0], color = 'black', linestyle = ":", label="All data"),
                        Line2D([0], [0], color = 'black', linestyle="--", label=r"P/$\delta$P > 3"), 
                        Line2D([0], [0], color = 'black', linestyle="-", label=r"P/$\delta$P > 3, $\Delta\theta < 3 \sigma$")]
                        #Line2D([0], [0], color = 'black', linestyle=":", label=r"P/$\delta$P > 3, $\Delta\theta < 30^\circ$")]

# Dust models
p_wls = np.array([154., 214., 850.])    # Values from model
p_mod_A = np.array([5.9, 9.3, 13.3])    #  G18 model A, G0 = 1.0
p_mod_B = np.array([6.0, 9.3, 12.9])    #  G18 model B, G0 = 1.0
p_mod_C = np.array([7.1, 9.8, 12.7])    #  G18 model C, G0 = 1.0
p_mod_D = np.array([9.5, 11.5, 13.1])   #  G18 model D, G0 = 1.0
p_modA_grid_154 = np.array([[3.2, 5.1, 7.3, 9.1, 10.8],   # for a_alig = 70 nm, G0 = .1, .3, 1., 3., 10.
                           [2.6, 4.2, 6.1, 7.8, 9.3],     # a_alig = 100 nm
                            [1.9, 3.3, 4.8, 6.2, 7.5]])   # a_alig = 150 nm 
p_modA_grid_214 = np.array([[7.6, 9.4, 11.1, 12.4, 13.4],
                           [6.4, 8.1,  9.6, 10.8, 11.8],
                            [5.0, 6.4,  7.7,  8.8,  9.6]])
p_modA_grid_850 = np.array([[16.3, 15.8, 15.3, 14.7, 14.1],
                           [14.6, 14.2, 13.7, 13.2, 12.6],
                            [12.1, 11.8, 11.4, 11.0, 10.5]])
p_modD_grid_154 = np.array([[7.8, 9.8, 11.0, 11.7, 12.1],
                           [5.9, 7.6,  8.8,  9.4,  9.8],
                            [3.8, 5.0,  5.9,  6.4,  6.7]])
p_modD_grid_214 = np.array([[12.1, 12.8, 13.2, 13.3, 13.1],
                           [ 9.5, 10.3, 10.7, 10.8, 10.7],
                            [ 6.4,  7.0,  7.3,  7.5,  7.5]])
p_modD_grid_850 = np.array([[16.3, 15.5, 14.7, 14.0, 13.2],
                           [13.6, 12.9, 12.2, 11.5, 10.9],
                            [ 9.6,  9.1,  8.6,  8.2,  7.7]])
G0_text_xpos = np.array([1., 1., 1., 1., 1.])
G0_text_ypos = np.array([1., 1., 1., 1., 1.])
alig_text_xpos = np.array([1., 1., 1.])
alig_text_ypos = np.array([.95, .95, .95])


# Set-up:
whatdata = 'e850' # 'de' # 'd850' # 
model2use = 'D' # 'none' # 'A' # 
plot_inset = False  # Default value; may be changed in the next cycle
central90as = True  # Only useful if using 850 um data
plot_lowess = False
absval = False

if whatdata == 'de':
    axlims = [0., 25., 0., 20.]
    if model2use == 'none':
        plot_inset = True
    else:
        pass
    axlims_inset = [.025, .475, .5, .5]
    first_lgd_pos = 'lower right'
    lowess_lgd_pos = 'upper right'
    ## No selection on P/dP
    mask_noPsel = mask_nosnrP_de
    x_noPsel = np.ma.masked_where(mask_noPsel, p_d.data)
    y_noPsel = np.ma.masked_where(mask_noPsel, p_e.data)
    #tau_d_noPsel = np.ma.masked_where(mask_nPsel, tau_map_d.data)  # Alternative 3rd parameters
    #T_noPsel = np.ma.masked_where(mask_noPsel, T_map.data)
    #dist_noPsel = np.ma.masked_where(mask_noPsel, distance_array*60.)
    delta_theta = delta_theta_de
    dtheta_noPsel = np.ma.masked_where(mask_noPsel, delta_theta)
    ## Selection on P/dP > 3
    mask_Psel = mask_de
    x_Psel = np.ma.masked_where(mask_Psel, p_d.data)
    y_Psel = np.ma.masked_where(mask_Psel, p_e.data)
    dtheta_Psel = np.ma.masked_where(mask_Psel, delta_theta)
    ## Selection on both P/dP > 3 and DTheta/dDtheta > 3
    mask_DthetaPsel = np.any([mask_dtheta_SNR_de, mask_de], axis = 0) # ([mask_dtheta_de, mask_de], axis = 0) # 
    x_DthetaPsel = np.ma.masked_where(mask_DthetaPsel, p_d.data)
    y_DthetaPsel = np.ma.masked_where(mask_DthetaPsel, p_e.data)
    ## Text
    xtext = r'$P$ (%) at 154 $\mu$m (HAWC+ band D)'
    ytext = r'$P$ (%) at 214 $\mu$m (HAWC+ band E)'
    titletext = r'N2071 pol. frac. at 154 $\mu$m vs. 214 $\mu$m'
    clabel = r'$\Delta \theta_{154-214}$' # r'T$_d$ (K)' # r'$\tau$ at 154 um' # 'Angular distance from center (arcmin)' #
    G0_text = np.array([r'0.1', '0.3', '1', '3', 'G$_0$ = 10'])
    alig_text = np.array([r'70 nm', '100 nm', 'a$_{alig}$ = 150 nm'])
    if model2use == 'A':
        G0_text = np.array([r'0.1', '0.3', '1', '3', 'G$_0$ = 10'])
        alig_text = np.array([r'70 nm', '100 nm', 'a$_{alig}$ = 150 nm'])
        G0_text_xpos = np.array([.5, .7, .92, .95, .95])
        G0_text_ypos = np.array([1.08, 1.07, 1.06, 1.05, 1.05])
        alig_text_xpos = np.array([1.02, 1.02, 1.02])
        alig_text_ypos = np.array([.96, .96, .98])
    elif model2use == 'D':
        G0_text = np.array([r'0.1', '0.3', '1', '3', '10'])
        alig_text = np.array([r'70 nm', '100 nm', '150 nm'])
        G0_text_xpos = np.array([.8, .85, .96, .97, 1.01])
        G0_text_ypos = np.array([1.04, 1.04, 1.04, 1.04, 1.04])
        alig_text_xpos = np.array([1.02, 1.02, 1.02])
        alig_text_ypos = np.array([.96, .96, .96])
    else:
        pass
    ## Model selection
    xind = 0
    yind = 1
    x_grid_A = p_modA_grid_154
    y_grid_A = p_modA_grid_214
    x_grid_D = p_modD_grid_154
    y_grid_D = p_modD_grid_214
elif whatdata == 'd850':
    if central90as == True:
        if model2use == 'none':
            axlims = [0., 2., 0., 5.]
        else:
            axlims = [0., 16., 0., 20.]
        mask_noPsel = mask_nosnrP_d850_90as # (no P selection)
        mask_Psel = mask_d850_90as          # (P/dP > 3)
        mask_DthetaPsel = np.any([mask_dtheta_SNR_d850, mask_d850_90as], axis = 0) # (P/dP > 3, matching theta)
    else:
        if model2use == 'none':
            axlims = [0., 8., 0., 6.]            
        else:
            axlims = [0., 16., 0., 20.]
        mask_noPsel = mask_nosnrP_d850 # (no P selection)
        mask_Psel = mask_d850          # (P/dP > 3)
        mask_DthetaPsel = np.any([mask_dtheta_SNR_d850, mask_d850], axis = 0) # (P/dP > 3, matching theta)
    #plot_inset = False
    #axlims_inset = [.025, .475, .5, .5]
    first_lgd_pos = 'upper left'
    lowess_lgd_pos = 'upper right'
    ## No selection on P/dP
    x_noPsel = np.ma.masked_where(mask_noPsel, p_d.data)
    y_noPsel = np.ma.masked_where(mask_noPsel, p_850.data)
    delta_theta = delta_theta_d850
    dtheta_noPsel = np.ma.masked_where(mask_noPsel, delta_theta)
    ## Selection on P/dP > 3
    x_Psel = np.ma.masked_where(mask_Psel, p_d.data)
    y_Psel = np.ma.masked_where(mask_Psel, p_850.data)
    dtheta_Psel = np.ma.masked_where(mask_Psel, delta_theta)
    ## Selection on both P/dP > 3 and DTheta/dDtheta > 3
    x_DthetaPsel = np.ma.masked_where(mask_DthetaPsel, p_d.data)
    y_DthetaPsel = np.ma.masked_where(mask_DthetaPsel, p_850.data)
    ## Text
    xtext = r'$P$ (%) at 154 $\mu$m (HAWC+ band D)'
    ytext = r'P (%) at 850 $\mu$m (POL-2)'
    titletext = r'N2071 pol. frac. at 154 $\mu$m vs. 850 $\mu$m'
    clabel = r'$\Delta \theta_{154-850}$'
    G0_text = np.array(['0.1', '0.3', '1', '3', '10'])
    alig_text = np.array(['70 nm', '100 nm', '150 nm'])
    if model2use == 'A':
        G0_text_xpos = np.array([.7, .8, .97, .98, .97])
        G0_text_ypos = np.array([1.04, 1.04, 1.04, 1.04, 1.04])
        alig_text_xpos = np.array([1., 1., 1.])
        alig_text_ypos = np.array([.95, .95, .95])
    elif model2use == 'D':
        G0_text_xpos = np.array([.88, .9, .98, .99, 1.])
        G0_text_ypos = np.array([1.03, 1.035, 1.04, 1.04, 1.04])
        alig_text_xpos = np.array([1., 1., 1.])
        alig_text_ypos = np.array([.97, .96, .95])
    else:
        pass
    ## Model selection
    xind = 0
    yind = 2
    x_grid_A = p_modA_grid_154
    y_grid_A = p_modA_grid_850
    x_grid_D = p_modD_grid_154
    y_grid_D = p_modD_grid_850
elif whatdata == 'e850':
    if central90as == True:
        if model2use == 'none':
            axlims = [0., 2., 0., 5.]
        else:
            axlims = [0., 16., 0., 20.]
        mask_noPsel = mask_nosnrP_e850_90as # (no P selection)
        mask_Psel = mask_e850_90as          # (P/dP > 3)
        mask_DthetaPsel = np.any([mask_dtheta_SNR_e850, mask_e850_90as], axis = 0) # (P/dP > 3, matching theta)
    else:
        if model2use == 'none':
            axlims = [0., 6., 0., 6.]
        else:
            axlims = [0., 16., 0., 20.]
        mask_noPsel = mask_nosnrP_e850 # (no P selection)
        mask_Psel = mask_e850          # (P/dP > 3)
        mask_DthetaPsel = np.any([mask_dtheta_SNR_e850, mask_e850], axis = 0) # (P/dP > 3, matching theta)
    #plot_inset = False
    #axlims_inset = [.025, .475, .5, .5]
    first_lgd_pos = 'upper right'
    lowess_lgd_pos = 'upper left'
    ## No selection on P/dP
    x_noPsel = np.ma.masked_where(mask_noPsel, p_e.data)
    y_noPsel = np.ma.masked_where(mask_noPsel, p_850.data)
    delta_theta = delta_theta_e850
    dtheta_noPsel = np.ma.masked_where(mask_noPsel, delta_theta)
    ## Selection on P/dP > 3
    x_Psel = np.ma.masked_where(mask_Psel, p_e.data)
    y_Psel = np.ma.masked_where(mask_Psel, p_850.data)
    dtheta_Psel = np.ma.masked_where(mask_Psel, delta_theta)
    ## Selection on both P/dP > 3 and Dtheta/dDtheta > 3
    x_DthetaPsel = np.ma.masked_where(mask_DthetaPsel, p_e.data)
    y_DthetaPsel = np.ma.masked_where(mask_DthetaPsel, p_850.data)
    ## Text
    xtext = r'$P$ (%) at 214 $\mu$m (HAWC+ band E)'
    ytext = r'P (%) at 850 $\mu$m (POL-2)'
    titletext = r'N2071 pol. frac. at 214 $\mu$m vs. 850 $\mu$m'
    clabel = r'$\Delta \theta_{214-850}$'
    G0_text = np.array(['0.1', '0.3', '1', '3', '10'])
    alig_text = np.array(['70 nm', '100 nm', '150 nm'])
    if model2use == 'A':
        G0_text_xpos = np.array([.87, .9, .98, .99, .99])
        G0_text_ypos = np.array([1.03, 1.035, 1.04, 1.04, 1.04])
        alig_text_xpos = np.array([1., 1., 1.])
        alig_text_ypos = np.array([.97, .96, .95])
    elif model2use == 'D':
        G0_text_xpos = np.array([.98, .99, 1.02, 1.03, 1.02])
        G0_text_ypos = np.array([1.04, 1.03, 1.02, 1., .98])
        alig_text_xpos = np.array([.98, .98, .98])
        alig_text_ypos = np.array([.95, .94, .93])
    else:
        pass
    ## Model selection
    xind = 1
    yind = 2
    x_grid_A = p_modA_grid_214
    y_grid_A = p_modA_grid_850
    x_grid_D = p_modD_grid_214
    y_grid_D = p_modD_grid_850
else:
    pass

# Nonparametric (LOWESS) models
lowess = sm.nonparametric.lowess
data_lowessmoothed_all = lowess(y_noPsel.ravel().compressed(), x_noPsel.ravel().compressed(), frac = .5)
x_lowessmoothed_all = data_lowessmoothed_all[:, 0]
y_lowessmoothed_all = data_lowessmoothed_all[:, 1]
data_lowessmoothed_Psel = lowess(y_Psel.ravel().compressed(), x_Psel.ravel().compressed(), frac = .5)
x_lowessmoothed_Psel = data_lowessmoothed_Psel[:, 0]
y_lowessmoothed_Psel = data_lowessmoothed_Psel[:, 1]
# ADD: Dtheta selection by sigma instead of abs. angle
data_lowessmoothed_DthetaPsel = lowess(y_DthetaPsel.ravel().compressed(), x_DthetaPsel.ravel().compressed(), frac = .5)
x_lowessmoothed_DthetaPsel = data_lowessmoothed_DthetaPsel[:, 0]
y_lowessmoothed_DthetaPsel = data_lowessmoothed_DthetaPsel[:, 1]
#data_lowessmoothed_DthetaPsel = lowess(y_DthetaPsel.ravel().compressed(), x_DthetaPsel.ravel().compressed(), frac = .5)
#x_lowessmoothed_DthetaPsel = data_lowessmoothed_DthetaPsel[:, 0]
#y_lowessmoothed_DthetaPsel = data_lowessmoothed_DthetaPsel[:, 1]



# Plot obs. data
fig, ax = plt.subplots(1, 1)
if absval == False:
    plt.scatter(x_noPsel, y_noPsel, marker="$\u25EF$", c = dtheta_noPsel, cmap = 'twilight_shifted', vmin = -90., vmax = 90.)
    plt.scatter(x_Psel, y_Psel, c = dtheta_Psel, cmap = 'twilight_shifted', vmin = -90., vmax = 90.)
else:
    plt.scatter(x_noPsel, y_noPsel, marker="$\u25EF$", c = np.abs(dtheta_noPsel), cmap = 'cividis_r', vmin = 0., vmax = 90.)
    plt.scatter(x_Psel, y_Psel, c = np.abs(dtheta_Psel), cmap = 'cividis_r', vmin = 0., vmax = 90.)
plt.xlabel(xtext)
plt.ylabel(ytext)
#plt.title(titletext)
plt.axis(axlims)
cb = plt.colorbar()
cb.set_label(clabel)
if whatdata == 'de' and model2use == 'none':
    first_legend = plt.legend(handles = LegendElements, loc = first_lgd_pos)
    lgd = plt.gca().add_artist(first_legend)
else:
    pass

# Tick + label size regulation
plt.xlabel(xtext, fontsize = 'xx-large')
plt.ylabel(ytext, fontsize = 'xx-large')
cb.set_label(clabel, fontsize = 'x-large')
plt.xticks(fontsize = 'x-large')
plt.yticks(fontsize = 'x-large')
plt.xticks(fontsize = 'x-large')
cb.ax.tick_params(labelsize = 'x-large')

# Plot LOWESS smoothing:
if plot_lowess == True:
    plt.plot(x_lowessmoothed_all, y_lowessmoothed_all, color = 'black', linestyle = ':')
    plt.plot(x_lowessmoothed_Psel, y_lowessmoothed_Psel, color = 'black', linestyle = '--')
    plt.plot(x_lowessmoothed_DthetaPsel, y_lowessmoothed_DthetaPsel, color = 'black')
    if whatdata == 'de':
        plt.legend(handles = LegendElementsLowess, loc = lowess_lgd_pos) # 'upper right') #

# Inset zoom plot
if plot_inset == True:
    axins = ax.inset_axes(axlims_inset)
    # XY coordinates followed by size (in fractions of original)
    # [.45, .45, .5, .5] for top right corner ([.475, .275, .5, .5] if there is a legend)
    # [.025, .475, .5, .5] for top left corner
    if absval == False:
        axins.scatter(x_noPsel, y_noPsel, marker="$\u25EF$", c = dtheta_noPsel, cmap = 'twilight_shifted', vmin = -90., vmax = 90.)
        axins.scatter(x_Psel, y_Psel, c = dtheta_Psel, cmap = 'twilight_shifted', vmin = -90., vmax = 90.)
    else:
        axins.scatter(x_noPsel, y_noPsel, marker="$\u25EF$", c = np.abs(dtheta_noPsel), cmap = 'cividis_r', vmin = 0., vmax = 90.)
        axins.scatter(x_Psel, y_Psel, c = np.abs(dtheta_Psel), cmap = 'cividis_r', vmin = 0., vmax = 90.)
    if plot_lowess == True:
        axins.plot(x_lowessmoothed_all, y_lowessmoothed_all, color = 'black', linestyle = ':')
        axins.plot(x_lowessmoothed_Psel, y_lowessmoothed_Psel, color = 'black', linestyle = '--')
        axins.plot(x_lowessmoothed_DthetaPsel, y_lowessmoothed_DthetaPsel, color = 'black')
    axins.set_xlim(0., 3.)
    axins.set_ylim(0., 3.)
    #axins.set_xticklabels('')
    axins.set_yticklabels('')
    
# Plot default dust models:
if model2use == 'A':
    p_mod = p_mod_A
    x_grid = x_grid_A
    y_grid = y_grid_A
    plotgrid = True
elif model2use == 'D':
    p_mod = p_mod_D
    x_grid = x_grid_D
    y_grid = y_grid_D
    plotgrid = True
else:
    plotgrid = False
    
if plotgrid == True:
    plt.plot(p_mod[xind], p_mod[yind], c = 'black', marker = 's')                                           # Basic model
    plt.plot(p_mod[xind]*np.array([0., 1.]), p_mod[yind]*np.array([0., 1.]), linestyle = ':', c = 'black')  # "Turbulence" line
    #plt.text(p_mod[xind], p_mod[yind] * 1.05, model2use, horizontalalignment = 'center')                    # Add model name
    for i, col in enumerate(x_grid):
        plt.plot(col, y_grid[i,:], 'k', markerfacecolor = 'none', markeredgecolor = 'k', marker = 's')      # Overplot model grid
    for j, dummy in enumerate(x_grid[0,:]):                                                                 # Annotations for G_0
        plt.text(x_grid[0,j] * G0_text_xpos[j], y_grid[0,j] * G0_text_ypos[j], G0_text[j], horizontalalignment = 'left',
                 c = 'red', size = 'x-large')
    for j, dummy in enumerate(x_grid[:,4]):                                                                 # Annotations for a_alig
        plt.text(x_grid[j,4] * alig_text_xpos[j], y_grid[j,4] * alig_text_ypos[j], alig_text[j],
                 verticalalignment = 'top', size = 'x-large')
    ModelLegend = [Line2D([0], [0], marker='s', color='k', linestyle="", label="Standard model"),           # Legend
                   Line2D([0], [0], marker='s', markerfacecolor='none', markeredgecolor='k', linestyle = '', label = 'Model variations')]
    plt.legend(handles = ModelLegend, loc = 'lower right')
#plt.savefig('p-d_vs_p-e_newpipeline_pix8')

'''
# Intermission: plot of I_wl1/I_wl2 ratio
x_Psel_1d = x_Psel.ravel().compressed()
y_Psel_1d = y_Psel.ravel().compressed()
x_DthetaPsel_1d = x_DthetaPsel.ravel().compressed()
y_DthetaPsel_1d = y_DthetaPsel.ravel().compressed()
ratio = y_Psel_1d / x_Psel_1d
ratio_Dtheta = y_DthetaPsel_1d / x_DthetaPsel_1d
plt.scatter(x_Psel_1d, ratio)
plt.scatter(x_DthetaPsel_1d, ratio_Dtheta)
'''


## Q(U)/I vs. Q(U)/I plots
#  (See V01)


## I vs. I plots
i_d_2plot = copy.deepcopy(i_d.data)
i_e_2plot = copy.deepcopy(i_e.data)
i_850_2plot = copy.deepcopy(i_850.data)
T_2plot = copy.deepcopy(T_map.data)
tau_d_2plot = copy.deepcopy(tau_map_d.data)
dist_2plot = copy.deepcopy(distance_array)*60.
dtheta_2plot = copy.deepcopy(delta_theta)  # Check which version you are using: 154-214 or 154-850
'''
i_d_2plot = np.ma.masked_where(mask_d850, i_d.data)
i_850_2plot = np.ma.masked_where(mask_d850, i_850.data)
T_2plot = np.ma.masked_where(mask_d850, T_map.data)
tau_d_2plot = np.ma.masked_where(mask_d850, tau_map_d.data)
dist_2plot = np.ma.masked_where(mask_d850, distance_array*60.)
dtheta_2plot = np.ma.masked_where(mask_d850, delta_theta)
'''
plt.scatter(i_d_2plot, i_e_2plot, c = dtheta_2plot, cmap = 'twilight_shifted')
#plt.scatter(i_d_2plot, i_850_2plot, c=T_2plot, cmap = 'rainbow')
#plt.scatter(i_d_2plot, i_850_2plot, c=tau_d_2plot, cmap = 'rainbow', norm=colors.LogNorm())
#plt.scatter(i_d_2plot, i_850_2plot, c=dist_2plot, cmap = 'rainbow')
plt.xscale('log')
plt.yscale('log')
plt.axis([1e-3, 3., 1e-3, 3.]) # ([1e-3, 3., 1e-5, 3e-2]) # 
plt.plot([2.5e-2, 2.5e-2], [1e-3, 3.], c = 'black', linestyle = '--') # [1e-5, 3e-2], c = 'black', linestyle = '--') # 
plt.plot([1e-3, 3.], [1.5e-2, 1.5e-2], c = 'black', linestyle = '--')
#plt.plot([1e-3, 3.], [2.5e-4, 2.5e-4], [1e-5, 3e-2], c = 'black', linestyle = '--')
plt.title('N2071 I at 154 um vs. 214 um') # 850 um') # 
plt.xlabel(r'I (154 $\mu$m) in Jy/arcsec$^2$')
plt.ylabel(r'I (214 $\mu$m) in Jy/arcsec$^2$') # (r'I (850 $\mu$m) in Jy/arcsec$^2$') # 
cb = plt.colorbar()
cb.set_label(r'$\Delta \theta (^\circ)$')
#cb.set_label(r'T$_d$ (K)')
#cb.set_label(r'$\tau$ at 154 um')
#cb.set_label('Angular distance from center (arcmin)')
#plt.savefig('N2071_i-d_vs_i-850_Tcoded_nosel')


# Angle vs. angle plots
theta_d_sel = np.ma.masked_where(mask_d, theta_d.data)
theta_d_3bdsel = np.ma.masked_where(mask_all, theta_d.data)
theta_e_sel = np.ma.masked_where(mask_e, theta_e.data)
theta_e_3bdsel = np.ma.masked_where(mask_all, theta_e.data)
theta_850_sel = np.ma.masked_where(mask_850, theta_850.data)
theta_850_3bdsel = np.ma.masked_where(mask_all, theta_850.data)
equline = np.linspace(-90., 90.)
###
plt.plot(equline, equline, color = 'black')
#plt.gca().set_prop_cycle(None)  # Resetting color cycle
ax = plt.gca()
ax.set_aspect(1.)
plt.xlim([-90., 90.])
plt.ylim([-90., 90.])
plt.fill_between(equline, equline-20., equline+20., alpha = 0.2, color = 'darkturquoise') # , color = 'grey') # 
plt.scatter(theta_d_sel, theta_e_sel, color ='darkturquoise') # (theta_d_sel, theta_850_sel, color ='grey') # 
plt.scatter(theta_d_3bdsel, theta_e_3bdsel, color ='orangered') # (theta_d_3bdsel, theta_850_3bdsel, color ='orangered') # 
#plt.title('Angle comparison') 
plt.xlabel(r'Pol. angle $\theta$ at 154 $\mu$m')
plt.ylabel(r'$\theta$ at 214 $\mu$m') # (r'$\theta$ at 850 $\mu$m') # 
#plt.savefig('Polangle_d-vs-850_comparison')

# Delta_theta vs. delta_theta
delta_theta_de = polangle_diff(angles_d.data, angles_e.data, liminf = -90., limsup = 90.)
delta_theta_de_masked = np.ma.masked_where(mask_all, delta_theta_de)
#delta_theta_de_masked2 = np.ma.masked_where(np.any([mask_all, mask_central2], axis = 0), delta_theta_de)
delta_theta_d850 = polangle_diff(angles_d.data, angles_850.data, liminf = -90., limsup = 90.)
delta_theta_d850_masked = np.ma.masked_where(mask_all, delta_theta_d850)
#delta_theta_d850_masked2 = np.ma.masked_where(np.any([mask_all, mask_central2], axis = 0), delta_theta_d850)

plt.scatter(delta_theta_de_masked, delta_theta_d850_masked)
plt.plot(equline, equline, color = 'black')
ax = plt.gca()
ax.set_aspect(1.)
plt.xlim([-90., 90.])
plt.ylim([-90., 90.])
plt.xlabel(r'$\Delta\theta$ (154 -- 214 $\mu$m)')
plt.ylabel(r'$\Delta\theta$ (154 -- 850 $\mu$m)')
#plt.savefig('Delta_Theta_Comparison')

# Delta_theta vs. SNR
whatmask = mask_d850 # _de # _all # 
#delta_theta_de = polangle_diff(angles_d.data, angles_e.data, liminf = -90., limsup = 90.)
#delta_theta_de_masked = np.ma.masked_where(mask_all, delta_theta_de)
delta_theta_d850 = polangle_diff(angles_d.data, angles_850.data, liminf = -90., limsup = 90.)
delta_theta_d850_masked = np.ma.masked_where(whatmask, delta_theta_d850)
#delta_theta_d850_masked = np.ma.masked_where(np.any([whatmask, p_850_SNR.data < 5], axis = 0), delta_theta_d850)
p_d_SNR = np.ma.masked_where(whatmask, p_d.data/dp_d.data)
p_e_SNR = np.ma.masked_where(whatmask, p_e.data/dp_e.data)
p_850_SNR = np.ma.masked_where(whatmask, p_850.data/dp_850.data)

# mask_PdP = np.any([whatmask, p_850_SNR.data < 30], axis = 0)
# delta_theta_d850_flat = delta_theta_d850_masked.ravel().compressed()
# p_850_SNR_flat = p_850_SNR.ravel().compressed()
# ... = np.ma.masked_where(mask_PdP, delta_theta_d850)
# np.isfinite(delta_theta_d850_flat).sum()
# (np.abs(delta_theta_d850_flat) > 30.).sum()
## Fraction of delta_theta_d850 > 30:
##   226/676 = 33.4 % (P/dP > 3)
##   139/441 = 31.5 % (P/dP > 5)
##    63/167 = 37.7 % (P/dP > 10)
##    36/ 68 = 52.9 % (P/dP > 15)
##    87/235 = 37.0 % (3 < P/dP < 5)
##    76/274 = 30.9 % (5 < P/dP < 10)
##    27/ 99 = 27.3 % (10 < P/dP < 15)
## Adding r < 6 arcmin condition:
##   

plt.scatter(p_850_SNR, np.abs(delta_theta_d850_masked))
plt.xscale('log')
plt.ylim([0., 90.])
plt.xlabel(r'$P/\delta P$ at 154 $\mu$m')
plt.ylabel(r'$\Delta\theta$ (154 - 850 $\mu$m)')
#plt.savefig('DTheta154-850_vs_SNR_P154')


# DTheta vs. col. dens.
sel1 = mask_de # mask_d850  # 
#sel2 = np.any([sel1, delta_theta_SNR_de < 3.], axis = 0) # _d850 < 3.], axis = 0) # 
sel2 = np.any([sel1, np.abs(delta_theta_de) < 30.], axis = 0) # _d850) < 30.], axis = 0) # 
x_1 = np.ma.masked_where(sel1, nh2_map.data) / 1e22
xx_1 = np.ma.masked_where(sel1, distance_array)
y_1 = np.ma.masked_where(sel1, np.abs(delta_theta_de)) # _d850)) # 
x_2 = np.ma.masked_where(sel2, nh2_map.data) / 1e22
xx_2 = np.ma.masked_where(sel2, distance_array)
y_2 = np.ma.masked_where(sel2, np.abs(delta_theta_de)) # _d850)) # 

fig, ax = plt.subplots(1, 1)
plt.scatter(x_1, y_1) # , marker="$\u25EF$")
plt.scatter(x_2, y_2)
plt.xscale('log')
#plt.title(r'Angle difference (154 - 850 $\mu$m) vs. column density')
plt.xlabel(r'$N_{H_2}$ in $10^{22}$ cm$^{-2}$') # (r'Distance from center (arcmin)') # 
plt.ylabel(r'$\Delta\theta$ (154 -- 214 $\mu$m)') # (r'$\Delta\theta$ (154 -- 850 $\mu$m)') # 
#plt.legend([r'$\Delta \theta/\delta \Delta \theta < 3$', r'$\Delta \theta/\delta \Delta \theta > 3$'])
plt.legend([r'$\Delta \theta < 30^\circ$', r'$\Delta \theta > 30^\circ$'])
#plt.savefig('DTheta_154-850_vs_coldens_30degsel')

y_1compress = y_1.ravel().compressed()
y_2compress = y_2.ravel().compressed()
sel3 = np.any([sel, delta_theta_SNR_d850 > 3.], axis = 0) # _de < 3.], axis = 0) # 
y_3 = np.ma.masked_where(sel3, np.abs(delta_theta_d850)) # _de)) # 
y_3compress = y_3.ravel().compressed()





### HISTOGRAMS ###

## Intensity
hi = plt.hist(i_d.data.ravel(), bins = 130, range = [-.1, 2.5], log = True)
plt.xlim([-.05, 2.1])
plt.title('D-band (154 um) intensities')
plt.xlabel('I (Jy / arcsec2)')
#plt.savefig('N2071_Dband_intensity_histogram')


## ANGLES
# Delta_Theta histograms

# Parameters
central90as = True  # Only relevant if using 850 um data
absval = False
whatdata = 'd850' # 'e850' # 'de' # 

# Preparation
if whatdata == 'de':
    titletext = r'154 and 214 $\mu$m bands'
    if absval == True:
        fname = 'Histogram_absDTheta_154-214_pix8_3sigma'
    else:
        fname = 'Histogram_DTheta_154-214_pix8_3sigma'
    # V01: selection by band(s)
    #delta_theta1 = np.ma.masked_where(selection1, delta_theta_de) # 
    #delta_theta2 = np.ma.masked_where(mask_all, delta_theta_de) #
    # V02: Selection by sigma threshold
    selection1 = mask_de
    dangleSNR = delta_theta_SNR_de
    selection2 = np.any([dangleSNR < 3., selection1], axis = 0)
    delta_theta1 = np.ma.masked_where(selection1, delta_theta_de)
    delta_theta2 = np.ma.masked_where(selection2, delta_theta_de)
elif whatdata == 'd850':
    titletext = r'154 and 850 $\mu$m bands'
    if absval == True:
        fname = 'Histogram_absDTheta_154-850_pix8_3sigma'
    else:
        fname = 'Histogram_DTheta_154-850_pix8_3sigma'
    # V02: Selection by sigma threshold
    if central90as == True:
        selection1 = mask_d850_90as
        fname += '_90as'
    else:
        selection1 = mask_d850
    dangleSNR = delta_theta_SNR_d850
    selection2 = np.any([dangleSNR < 3., selection1], axis = 0)
    delta_theta1 = np.ma.masked_where(selection1, delta_theta_d850)
    delta_theta2 = np.ma.masked_where(selection2, delta_theta_d850)
elif whatdata == 'e850':
    titletext = r'214 and 850 $\mu$m bands'
    if absval == True:
        fname = 'Histogram_absDTheta_214-850_pix8_3sigma'
    else:
        fname = 'Histogram_DTheta_214-850_pix8_3sigma'
    # V02: Selection by sigma threshold
    if central90as == True:
        selection1 = mask_e850_90as
        fname += '_90as'
    else:
        selection1 = mask_e850
    dangleSNR = delta_theta_SNR_e850
    selection2 = np.any([dangleSNR < 3., selection1], axis = 0)
    delta_theta1 = np.ma.masked_where(selection1, delta_theta_e850)
    delta_theta2 = np.ma.masked_where(selection2, delta_theta_e850)
else:
    pass

if absval == True:
    thetarange = (0., 90.)
    binnum = 9
    delta_theta1_1d = np.abs(delta_theta1.ravel().compressed())
    delta_theta2_1d = np.abs(delta_theta2.ravel().compressed())
else:
    thetarange = (-90., 90.)
    binnum = 18
    delta_theta1_1d = delta_theta1.ravel().compressed()
    delta_theta2_1d = delta_theta2.ravel().compressed()

# Actual plot
h = plt.hist(delta_theta1_1d, range = thetarange, bins = binnum, histtype = 'step', fill = True)
plt.hist(delta_theta2_1d, range = thetarange, bins = binnum, histtype = 'step', fill = True)
if absval == False:
    plt.plot([-30., -30.], [0., 1.05 * max(h[0])], color = 'black', linestyle = ':')
plt.plot([30., 30.], [0., 1.05 * max(h[0])], color = 'black', linestyle = ':')
plt.title(titletext, fontsize = 'xx-large')
plt.xlabel(r'$\Delta \theta \, (^\circ)$', fontsize = 'x-large')
plt.ylabel('Pixel count', fontsize = 'x-large') 
#plt.xticks(fontsize = 'x-large')
#plt.yticks(fontsize = 'x-large')
#plt.savefig(fname)




# Some accounting
# mask_toto = np.any([mask_d850, np.abs(delta_theta_d850) < 20, delta_theta_SNR_d850 < 3], axis = 0)
# thing = np.ma.masked_where(mask_toto, delta_theta_d850)
# other = thing.ravel().compressed()
# N of finite elements in array: np.isfinite(array).sum()
# N of elements > 20 in array: (np.abs(array) > 20.).sum()
# # To check if extra SNR selection is making any difference:
#   testmask = np.any([mask_d850, np.abs(delta_theta_d850.data) < 30.], axis = 0)
#   test = np.ma.masked_where(testmask, delta_theta_SNR_d850)
#   np.nanmin(test.ravel().compressed())
# # If you want an ordered list to check consecutive values: s = set(test)
#
# 4 ARCSEC PIXELS
# Selection on both DTheta > 20/30 and SNR(DTheta) > 3:
#   DTheta 154-214: 142 / 617 are > 20;  82 / 340 for 3-band selection
#                    62 / 617 are > 30;  50 / 340 for 3-band sel.
#   DTheta 154-850: 328 / 676 are > 20; 184 / 340 for 3-band sel.
#                   226 / 676 are > 30; 116 / 340 for 3-band sel.
#   DTheta 214-850: 134 / 399 are > 30; 125 / 340 for 3-band sel.
# Same, but for data with P/dP >= 15
#   DTheta 154-214:  2 / 59, min SNR (for DTheta) = 24.63
#   DTheta 154-850: 36 / 67, min SNR (for DTheta) = 15.14
#   DTheta 214-850: 23 / 34, min SNR (for DTheta) = 17.99
# Same, but for data with P/dP >= 30
#   DTheta 154-214: 0 / 20
#   DTheta 154-850: 5 /  6, min SNR (for DTheta) = 42.9
#   DTheta 214-850: 1 /  1, min SNR (for DTheta) = 18.9
## Selection on DTheta alone
#   DTheta 154-214: 163 / 617 are > 20;  95 / 340 for 3-band sel.
#                    62 / 617 are > 30;  50 / 340 for 3-band sel.
#   DTheta 154-850: 338 / 676 are > 20; 187 / 340 for 3-band sel.
#                   226 / 676 are > 30; 116 / 340 for 3-band sel.
#   DTheta 214-850: 135 / 399 are > 30; 126 / 340 for 3-band sel.,  min SNR = 2.68, 2nd-min = 3.03
# Selection on SNR alone
#   DTheta 154-214: 235 / 617; 158 / 340 for 3-band sel.
#   DTheta 154-850: 412 / 676; 237 / 340 for 3-band sel.
# 
## ADD SELECTION ON CENTRAL 6'
#  ...
#
# 12 ARCSEC PIXELS
# Selection on both DTheta > 30 and SNR > 9:
#   DTheta 154-214: 22 / 130
#   DTheta 154-850: 37 / 101
#   DTheta 214-850: 33 / 125
# Selection on DTheta > 30 alone
#   DTheta 154-214: 29 / 130, min SNR = 4.69
#   DTheta 154-850: 44 / 101, min SNR = 4.15
#   DTheta 214-850: 56 / 125, min SNR = 3.50


## Pol. Frac.
selection = mask_all #mask_de #  
pmasked = np.ma.masked_where(selection, p_850.data) #delta_theta_de) #  
pmasked_1d = pmasked.ravel().compressed()
np.min(pmasked_1d), np.median(pmasked_1d), np.max(pmasked_1d)
np.isfinite(pmasked_1d).sum(), (pmasked_1d > 20.).sum()
# No threshold on I:
#  P max:     187% (154 um),   41% (214 um), 28  % (850 um)
#  P median: 14.9% (154 um), 10.3% (214 um),  2.9% (850 um)
#  P >  20%:  36.8% = 4259 / 11562 (154 um), 14.3% = 301 / 2107 (214 um), .5% = 11 / 2280 (850 um)
#    >  30%:  20.0% = 2318 / 11562 (154 um),  1.9% =  40 / 2107 (214 um)
#    >  40%:  10.6% = 1231 / 11562 (154 um),   .1% =   2 / 2107 (214 um)
#    >  50%:   5.6% =  645 / 11562 (154 um)
#    > 100%:    .4% =   47 / 11562 (154 um)
# Threshold on I, single-band selection:
#  P max:    24.9% (154 um), 17.5% (214 um), 5.7% (850 um)
#  P median:  4.2% (154 um),  4.2% (214 um), 1.5% (850 um)
#  P > 20%: .5% = 8 / 1584 (154 um), 0.0% = 0 / 771 (214 um), 0.0% = 0 / 1183 (850 um)
# Threshold on I, 3-band selection:
#  P max:    12.59% (154 um), 7.30% (214 um), 3.84% (850 um)
#  P median:   .85% (154 um), 1.49% (214 um),  .86% (850 um)
#  P > 20%:  None





### FITS ###

# D vs. E band (for comparison with Michail+21)
#selection = mask_de
#selection = np.any([mask_dtheta_SNR_de, mask_de], axis = 0)
selection = np.any([mask_de, peak_distance_array > 45./3600.], axis = 0)
#selection = np.any([mask_dtheta_SNR_de, mask_de, peak_distance_array > 45./3600.], axis = 0)
x_all = p_d.data
y_all = p_e.data
dy_all = dp_e.data
x_sel = np.ma.masked_where(selection, x_all)
y_sel = np.ma.masked_where(selection, y_all)
dy_sel = np.ma.masked_where(selection, dy_all)
x = x_sel.ravel().compressed()
y = y_sel.ravel().compressed()
dy = dy_sel.ravel().compressed()

# Slope
from subroutines_mapmaking import line
from scipy.optimize import curve_fit
res, cov = curve_fit(line, x, y, sigma = dy)
errs = np.sqrt(np.diag(cov))

# Plot
x_recal = np.linspace(0., 25., 26)
y_recal = x_recal * res[0] + res[1]
plt.scatter(x, y)
plt.plot(x_recal, y_recal, color = 'black', linestyle = '--')

# Median stats
from scipy.stats import median_abs_deviation as MAD
np.median(y/x)
MAD(y/x)

# Results:
# No Dtheta sel.
#   Slope:              0.42 +- 0.01
#   Ratio (everywhere): 0.63 +- 0.23  (median, MAD)
#   Ratio (< 90 as):    1.21 +- 0.42
#   Ratio (> 90 as):    0.49 +- 0.10
# Dtheta < 3 sig.
#   Slope:              0.38 +- 0.01
#   Ratio (everywhere): 0.48 +- 0.10
#   Ratio (< 90 as):    1.77 +- 0.64
#   Ratio (> 90 as):    0.45 +- 0.06


# Older version (P vs. I comparison with Aran) -- see V01
