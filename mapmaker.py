
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from astropy.io import fits
#from reproject import reproject_interp
#from astropy.wcs import WCS
#import copy
from subroutines_mapmaking import convertfiles
from subroutines_mapmaking import readfilesfinal
from subroutines_mapmaking import make_p_theta
from subroutines_mapmaking import plot_levels
from subroutines_mapmaking import plot_I_polvec


folder = ''


# INITIAL FITS FILES FORMATION/REFORMATTING
# No need to run this if you already did once

# N2071 data filenames
## HAWC+ files
filed = folder + 'SOFIA_retrievals_21-09-06/data/SOFIA/HAWC_PLUS/L4/p10196/data/g5/F0615_HA_POL_0701308_HAWDHWPD_PMP_071-114.fits'
filee = folder + 'SOFIA_retrievals_21-09-06/data/SOFIA/HAWC_PLUS/L4/p10197/data/g1/F0616_HA_POL_0701309_HAWEHWPE_PMP_067-102.fits'
''' Old pipeline
filed = folder + 'SOFIA_retrievals_19-12-19/sofia_2019-10-02_HA_F615/p7514/F0615_HA_POL_0701308_HAWDHWPD_PMP_071-114.fits'    #N2071 D
filee = folder + 'SOFIA_retrievals_19-12-19/sofia_2019-10-02_HA_F615-2/p7514/F0615_HA_POL_0701307_HAWEHWPE_PMP_067-102.fits'  #N2071 E
filed = folder + 'SOFIA_retrievals_19-12-19/sofia_2019-10-02_HA_F621-2/p7526/F0621_HA_POL_0701306_HAWDHWPD_PMP_033-049.fits' #Serpens D
filee = folder + 'SOFIA_retrievals_19-12-19/sofia_2019-10-10_HA_F621/p7526/F0621_HA_POL_0701305_HAWEHWPE_PMP_000-032.fits'   #Serpens E
'''
## JCMT 850 um files
# 4 arcsec pixel version
files_850 = [folder + 'NGC_2071/NGC2071_4aspx_maps_2020-09/iext_all_subco_725.fits',  # I
             folder + 'NGC_2071/NGC2071_4aspx_maps_2020-09/qext_all_725.fits',        # Q
             folder + 'NGC_2071/NGC2071_4aspx_maps_2020-09/uext_all_725.fits']        # U
# 8 arcsec pixel version
files_850 = [folder + 'JCMT-iext_850_8pixel.fits',
             folder + 'JCMT-qext_850_8pixel.fits',
             folder + 'JCMT-uext_850_8pixel.fits']
'''
# Older 4 arcsec pixel version (until 2020-09-13)
files_850 = [folder + 'NGC_2071/NGC2071_4aspx_maps_2020-07/iext_850_4p.fits',  # I
             folder + 'NGC_2071/NGC2071_4aspx_maps_2020-07/qext_850_4p.fits',  # Q
             folder + 'NGC_2071/NGC2071_4aspx_maps_2020-07/uext_850_4p.fits']  # U
# 12 arcsec pixel version
files_850 = [folder + 'NGC_2071/NGC2071_12aspx_maps_2020-09/iext_12pixel.fits',  # I
             folder + 'NGC_2071/NGC2071_12aspx_maps_2020-09/qext_12pixel.fits',  # Q
             folder + 'NGC_2071/NGC2071_12aspx_maps_2020-09/uext_12pixel.fits']  # U
'''
# Ancillary Herschel files
file_h160 = folder + 'Ancillary_Herschel_data/N2071_160um_cutout.fits' # 'Ancillary_Herschel_data/orionB-S-160.fits' # 
file_h250 = folder + 'Ancillary_Herschel_data/N2071_250um_cutout.fits' # 'Ancillary_Herschel_data/orionB-250.fits'   # 


# Make files to be converted to NDF and smoothed by Starlink
imap_d, qmap_d, umap_d = convertfiles(filed, instru = 'HAWC+', units = 'Jy/arcsec2')#, savefile = 'N2071_HAWC+D')
imap_e, qmap_e, umap_e = convertfiles(filee, instru = 'HAWC+', units = 'Jy/arcsec2')#, savefile = 'N2071_HAWC+E')
imap_850, qmap_850, umap_850 = convertfiles(files_850, instru = 'SCUBA2', units = 'Jy/arcsec2', beam = 14.1, savefile = 'N2071_JCMT-850-12as')
'''
# Older version
imap_d, qmap_d, umap_d = convertfiles(filed, instru = 'HAWC+')#, savefile = 'N2071_HAWC+D')
imap_e, qmap_e, umap_e = convertfiles(filee, instru = 'HAWC+')#, savefile = 'N2071_HAWC+E')
imap_850, qmap_850, umap_850 = convertfiles(files_850, instru = 'SCUBA2')#, savefile = 'N2071_JCMT-850')
#imap_h160 = convertfiles(file_h160, instru = 'Herschel')  # Only needed if you are using the full image
#imap_h250 = convertfiles(file_h250, instru = 'Herschel')  # Only needed if you are using the full image
'''


### AFTER RUNNING THE CODE ABOVE, RUN THE STARLINK COMMANDS FROM 'Starlink_commands_full_sample.txt' ###
### THE CODE BELOW IS MEANT TO BE USED ONLY AFTER STARLINK HAS BEEN RUN AT LEAST ONCE ###


# Read starlink-smoothed files, convert units
## Unsmoothed version (for consistency test)
files_d_raw = [folder + 'N2071_HAWC+D_i+var_Jy-arcs2.fits',  # Filenames
                  folder + 'N2071_HAWC+D_q+var_Jy-arcs2.fits',
                  folder + 'N2071_HAWC+D_u+var_Jy-arcs2.fits'] 
files_e_raw = [folder + 'N2071_HAWC+E_i+var_Jy-arcs2.fits', 
                  folder + 'N2071_HAWC+E_q+var_Jy-arcs2.fits', 
                  folder + 'N2071_HAWC+E_u+var_Jy-arcs2.fits']
files_850_raw = [folder + 'N2071_JCMT-850_i+var_Jy-arcs2.fits', 
                    folder + 'N2071_JCMT-850_q+var_Jy-arcs2.fits', 
                    folder + 'N2071_JCMT-850_u+var_Jy-arcs2.fits']
i_d_raw, q_d_raw, u_d_raw, di_d_raw, dq_d_raw, du_d_raw = readfilesfinal(files_d_raw, units_in = 'Jy/arcsec2')  # Reading + unit conversion
i_e_raw, q_e_raw, u_e_raw, di_e_raw, dq_e_raw, du_e_raw = readfilesfinal(files_e_raw, units_in = 'Jy/arcsec2')
i_850_raw, q_850_raw, u_850_raw, di_850_raw, dq_850_raw, du_850_raw = readfilesfinal(files_850_raw, units_in = 'Jy/arcsec2')
# Smoothed + regridded version
# Rebinned version ('rebin = yes')
files_d_smooth = [folder + 'N2071_HAWC+D_i+var_Jy-arcs2_beam18.2_pix4.0-sincsinc-rb.fits',
                  folder + 'N2071_HAWC+D_q+var_Jy-arcs2_beam18.2_pix4.0-sincsinc-rb.fits',
                  folder + 'N2071_HAWC+D_u+var_Jy-arcs2_beam18.2_pix4.0-sincsinc-rb.fits'] 
files_e_smooth = [folder + 'N2071_HAWC+E_i+var_Jy-arcs2_beam18.2_pix4.0-sincsinc-rb.fits', 
                  folder + 'N2071_HAWC+E_q+var_Jy-arcs2_beam18.2_pix4.0-sincsinc-rb.fits', 
                  folder + 'N2071_HAWC+E_u+var_Jy-arcs2_beam18.2_pix4.0-sincsinc-rb.fits']
'''
# Resampled version ('rebin = no')
files_d_smooth = [folder + 'N2071_HAWC+D_i+var_Jy-arcs2_beam18.2_pix4.0-sincsinc-rb.fits',
                  folder + 'N2071_HAWC+D_q+var_Jy-arcs2_beam18.2_pix4.0-sincsinc-rb.fits',
                  folder + 'N2071_HAWC+D_u+var_Jy-arcs2_beam18.2_pix4.0-sincsinc-rb.fits'] 
files_e_smooth = [folder + 'N2071_HAWC+E_i+var_Jy-arcs2_beam18.2_pix4.0-sincsinc-rb.fits', 
                  folder + 'N2071_HAWC+E_q+var_Jy-arcs2_beam18.2_pix4.0-sincsinc-rb.fits', 
                  folder + 'N2071_HAWC+E_u+var_Jy-arcs2_beam18.2_pix4.0-sincsinc-rb.fits']
'''
files_850_smooth = [folder + 'N2071_JCMT-850_i+var_Jy-arcs2_beam18.2_pix4.0.fits',      # JCMT files have no regridding anyway
                    folder + 'N2071_JCMT-850_q+var_Jy-arcs2_beam18.2_pix4.0.fits', 
                    folder + 'N2071_JCMT-850_u+var_Jy-arcs2_beam18.2_pix4.0.fits']
i_d, q_d, u_d, di_d, dq_d, du_d = readfilesfinal(files_d_smooth, units_in = 'Jy/arcsec2')
i_e, q_e, u_e, di_e, dq_e, du_e = readfilesfinal(files_e_smooth, units_in = 'Jy/arcsec2')
i_850, q_850, u_850, di_850, dq_850, du_850 = readfilesfinal(files_850_smooth, units_in = 'Jy/arcsec2')
'''
# Previous version
i_d, q_d, u_d, di_d, dq_d, du_d = readfilesfinal(files_d_smooth, beam = 13.8)
i_e, q_e, u_e, di_e, dq_e, du_e = readfilesfinal(files_e_smooth, beam = 18.2)
#i_450, q_450, u_450, di_450, dq_450, du_450 = readfilesfinal(files_450_smooth, beam = 9.8, save='N2071_JCMT-450')
i_850, q_850, u_850, di_850, dq_850, du_850 = readfilesfinal(files_850_smooth, beam = 14.1)
'''

# Create polarized quantity maps, save as files if requested
# Unsmoothed
Ip_d_raw, p_d_raw, theta_d_raw, dIp_d_raw, dp_d_raw, dtheta_d_raw = make_p_theta(i_d_raw, q_d_raw, u_d_raw, di_d_raw, dq_d_raw, du_d_raw)
Ip_e_raw, p_e_raw, theta_e_raw, dIp_e_raw, dp_e_raw, dtheta_e_raw = make_p_theta(i_e_raw, q_e_raw, u_e_raw, di_e_raw, dq_e_raw, du_e_raw)
Ip_850_raw, p_850_raw, theta_850_raw, dIp_850_raw, dp_850_raw, dtheta_850_raw = make_p_theta(i_850_raw, q_850_raw, u_850_raw, 
                                                                                             di_850_raw, dq_850_raw, du_850_raw)
# Smoothed
Ip_d, p_d, theta_d, dIp_d, dp_d, dtheta_d = make_p_theta(i_d, q_d, u_d, di_d, dq_d, du_d)#, save='N2071_HAWC+D')
Ip_e, p_e, theta_e, dIp_e, dp_e, dtheta_e = make_p_theta(i_e, q_e, u_e, di_e, dq_e, du_e)#, save='N2071_HAWC+E')
Ip_850, p_850, theta_850, dIp_850, dp_850, dtheta_850 = make_p_theta(i_850, q_850, u_850, di_850, dq_850, du_850)#, save='N2071_JCMT-850')


# PLOTS:

# A1. Plotting S/N before/after smoothing, HAWC+ D band
fig_SNR_d_1 = plot_levels(i_d_raw, di_d_raw, title = 'N2071 I with SNR contours, pre-smoothing (HAWC+ D, 154 um)', 
                          contours = 'SNR', crop = -1., showbeam = False)
fig_SNR_d_2 = plot_levels(i_d, di_d, title = 'N2071 I with SNR contours,  after smoothing (HAWC+ D, 154 um)', 
                          contours = 'SNR', crop = -1., showbeam = False)
#np.nanmax(i_d_raw.data), np.nanmax(i_d.data)
#np.nanmean(i_d_raw.data), np.nanmean(i_d.data)
#fig_SNR_d_1.save('N2071_HAWC+154um_before-smoothing_contours.png')
#fig_SNR_d_2.save('N2071_HAWC+154um_after-smoothing_contours.png')

# A2. Plotting S/N before/after smoothing, HAWC+ E band
fig_SNR_e_1 = plot_levels(i_e_raw, di_e_raw, title = 'N2071 I with SNR contours, pre-smoothing (HAWC+ E, 214 um)', 
                          contours = 'SNR', crop = -1., showbeam = False)
fig_SNR_e_2 = plot_levels(i_e, di_e, title = 'N2071 I with SNR contours,  after smoothing (HAWC+ E, 214 um)', 
                          contours = 'SNR', crop = -1., showbeam = False)
#np.nanmax(i_e_raw.data), np.nanmax(i_e.data)
#np.nanmean(i_e_raw.data), np.nanmean(i_e.data)
#fig_SNR_e_1.save('N2071_HAWC+214um_before-smoothing_contours.png')
#fig_SNR_e_2.save('N2071_HAWC+214um_after-smoothing_contours.png')

# A3. Plotting S/N before/after smoothing, SCUBA2 850 um band
fig_SNR_850_1 = plot_levels(i_850_raw, di_850_raw, title = 'N2071 I with SNR contours, pre-smoothing (HAWC+ E, 214 um)', 
                          contours = 'SNR', crop = 7.5, showbeam = False)
fig_SNR_850_2 = plot_levels(i_850, di_850, title = 'N2071 I with SNR contours,  after smoothing (HAWC+ E, 214 um)', 
                          contours = 'SNR', crop = 7.5, showbeam = False)
#np.nanmax(i_850_raw.data), np.nanmax(i_850.data)
#np.nanmean(i_850_raw.data), np.nanmean(i_850.data)
#fig_SNR_850_1.save('N2071_850um_before-smoothing_contours.png')
#fig_SNR_850_2.save('N2071_850um_after-smoothing_contours.png')


# B. Contour levels of Herschel (PACS) map
imap_h160 = fits.open(file_h160) # Ancillary Herschel data
imap_h250 = fits.open(file_h250)
fig_contours_h160 = plot_levels(imap_h160, title = 'N2071 contours for I (PACS 160 um)', contours = 'I', 
                                contour_values = [1e-3, 3e-3, 1e-2], colors = ['white', 'darkgrey', 'black'], 
                                crop = -1., showbeam = False)
#fig_contours_h160.save('N2071_PACS_contours.png')
fig_contours_h250 = plot_levels(imap_h250, title = 'N2071 contours for I (SPIRE 250 um)', contours = 'I',
                                contour_values = [1e-3, 3e-3, 1e-2], colors = ['white', 'darkgrey', 'black'], 
                                crop = -1., showbeam = False)
#fig_contours_h250.save('N2071_SPIRE_contours.png')


# C.A Plotting I map + pol. vectors them with PLOT_I_POLVEC, save figures if needed.
fig_850 = plot_I_polvec(i_850, di_850, p_850, dp_850, theta_850, dtheta_850,
                                                      title = 'POL2 850 um', i_SNR_cut = 20., p_SNR_cut = 3., # 
                                                      crop = 7.5, scalevec = .6, stepvec = 2)
#fig_850.save('N2071_jcmt850um_I+polvec_SNR-I20-p3.png')
fig_d = plot_I_polvec(i_d, di_d, p_d, dp_d, theta_d, dtheta_d, title = 'HAWC+ band D (154 um)', i_SNR_cut = 20.,
                                                   thresh = 10e-3, p_SNR_cut = 3., showbeam = False, scalevec = .6, stepvec = 3) # 
#fig_d.save('N2071_HAWC+154um_I+polvec_SNR-I20-p3_thr10.png')
fig_e = plot_I_polvec(i_e, di_e, p_e, dp_e, theta_e, dtheta_e, title = 'HAWC+ band E (214 um)', i_SNR_cut = 20.,
                                                   thresh = 10e-3, p_SNR_cut = 3., showbeam = False, scalevec = .6, stepvec = 3) # 
#fig_e.save('N2071_HAWC+214um_I+polvec_SNR-I20-p3_thr10.png')

