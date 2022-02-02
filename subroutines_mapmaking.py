
# For debugging:
# pdb.set_trace()

import numpy as np
from astropy.io import fits
from aplpy import FITSFigure
import matplotlib.pyplot as plt
import copy
import pdb



### READ ORIGINAL FILES AND CONVERT INTO PROPER UNITS ###

def convertfiles(fname, instru = 'undefined', units = 'Jy/arcsec2', beam = -1., savefile = ''):
    # HAWC+ & JCMT data: Reads the original files and converts them in a Starlink-usable format
    #  (HDUList object for each Stokes parameter, with variance map)
    # Herschel data (ancillary, no need for smoothing): Extract and convert to Jy/arcsec2

    instru = instru.upper()
    if instru == 'HAWC+' or instru == 'SOFIA':
        # Need to: convert error --> variance, units from Jy/pixel --> whatever is asked
        fcontent = fits.open(fname)
        # Compute factor for unit conversion. HAWC+ native unit: Jy/pixel
        pxscale = fcontent['STOKES I'].header['CDELT2']*3600  # map scale in arcsec/pix. CDELT2 always in deg
        if units.upper() == 'JY/ARCSEC2':
            convfact = 1./pxscale**2
            unit_string = 'Jy-arcs2'
        elif units.upper() == 'JY/BEAM':
            pxscale = fcontent['STOKES I'].header['CDELT2']*3600  # map scale in arcsec/pix. CDELT2 always in deg
            beam = float(fcontent['STOKES I'].header['BEAM'])     # Beam FWHM in arcsec
            beam_sig = beam / 2*np.sqrt(2*np.log(2))              # FWHM --> sigma (assuming Gaussian beam)
            beam_area = 2 * np.pi * beam_sig**2                   # Beam area in arcsec^2
            convfact = beam_area / pxscale**2                     # Beam area in pixel^2
            unit_string = 'Jy-beam'
        else:
            print()
            print('ERROR: requested output unit {} not available. Please choose "Jy/arcsec2" or "Jy/beam".'.format(units))
            print()
        
        # Reading and converting IQU maps
        # Nota: Only 'STOKES I' is a PrimaryHDU. Converting 'STOKES Q/U' to PrimaryHDU as well
        ## Data
        i = fcontent['STOKES I']
        i.data *= convfact
        i.header['BUNIT'] = units
        i.header['EXTNAME'] = 'PRIMARY'
        q = fits.PrimaryHDU()
        q.data = fcontent['STOKES Q'].data * convfact
        q.header = copy.deepcopy(i.header)
        u = fits.PrimaryHDU()
        u.data = fcontent['STOKES U'].data * convfact
        u.header = copy.deepcopy(i.header)
        ## Variances
        di = fits.ImageHDU()
        di.data = (fcontent['ERROR I'].data**2) * (convfact**2)
        di.header = fcontent['ERROR I'].header
        di.header['BUNIT'] = units + '**2'
        di.header['EXTNAME'] = 'VARIANCE'
        dq = fits.ImageHDU()
        dq.data = (fcontent['ERROR Q'].data**2) * (convfact**2)
        dq.header = copy.deepcopy(di.header)
        du = fits.ImageHDU()
        du.data = (fcontent['ERROR U'].data**2) * (convfact**2)
        du.header = copy.deepcopy(di.header)
        
        # Save maps if requested
        imap = fits.HDUList([i, di])
        qmap = fits.HDUList([q, dq])
        umap = fits.HDUList([u, du])
        if savefile != '':
            imap.writeto(savefile + '_i+var_' + unit_string + '.fits', overwrite=True)
            qmap.writeto(savefile + '_q+var_' + unit_string + '.fits', overwrite=True)
            umap.writeto(savefile + '_u+var_' + unit_string + '.fits', overwrite=True)
        else:
            print('Files not saved (specify filename as SAVEFILE keyword if you wanted to)')
        
        # Output
        return imap, qmap, umap

    elif instru == 'JCMT' or (instru == 'SCUBA2' or instru == 'POL2'):
        # Read data
        i_content = fits.open(fname[0])
        q_content = fits.open(fname[1])
        u_content = fits.open(fname[2])

        #pdb.set_trace()

        # Native units: Jy/beam

        # If array is 3D (we assume if one is, all are):
        # Reformat (3D --> 2D) Stokes parameters and variance maps
        if i_content['PRIMARY'].data.ndim == 3:
            i_out = fits.PrimaryHDU()
            i_out.data = i_content['PRIMARY'].data[0,:,:]   # This automatically creates a header
            i_out.header = i_content['PRIMARY'].header      # Copied header is atomatically edited to be 2D
            di_out = fits.ImageHDU()
            di_out.data = i_content['VARIANCE'].data[0,:,:]
            di_out.header = i_content['VARIANCE'].header
            q_out = fits.PrimaryHDU()
            q_out.data = q_content['PRIMARY'].data[0,:,:]
            q_out.header = i_out.header
            dq_out = fits.ImageHDU()
            dq_out.data = q_content['VARIANCE'].data[0,:,:]
            dq_out.header = di_out.header
            u_out = fits.PrimaryHDU()
            u_out.data = u_content['PRIMARY'].data[0,:,:]
            u_out.header = i_out.header
            du_out = fits.ImageHDU()
            du_out.data = u_content['VARIANCE'].data[0,:,:]
            du_out.header = di_out.header
            #unit_string = ''  # To be updated
        else:
            i_out = fits.PrimaryHDU()
            i_out.data = i_content['PRIMARY'].data
            i_out.header = i_content['PRIMARY'].header
            di_out = fits.ImageHDU()
            di_out.data = i_content['VARIANCE'].data
            di_out.header = i_content['VARIANCE'].header
            q_out = fits.PrimaryHDU()
            q_out.data = q_content['PRIMARY'].data
            q_out.header = i_out.header
            dq_out = fits.ImageHDU()
            dq_out.data = q_content['VARIANCE'].data
            dq_out.header = di_out.header
            u_out = fits.PrimaryHDU()
            u_out.data = u_content['PRIMARY'].data
            u_out.header = i_out.header
            du_out = fits.ImageHDU()
            du_out.data = u_content['VARIANCE'].data
            du_out.header = di_out.header
            #unit_string = ''  # To be updated

        if units.upper() == 'JY/ARCSEC2':
            if beam <= 0.:
                print('ERROR: Variable "beam" is needed to convert JCMT maps to units of Jy/arcsec^2.')
            beam_sig = beam / 2*np.sqrt(2*np.log(2))   # FWHM --> sigma (assuming Gaussian beam)
            beam_area = 2 * np.pi * beam_sig**2        # Beam area in arcsec^2
            convfact = 1./beam_area
            i_out.data *= convfact
            di_out.data *= convfact**2
            i_out.header['BUNIT'] = units
            di_out.header['BUNIT'] = '(' + units + ')' + '**2'
            q_out.data *= convfact
            dq_out.data *= convfact**2
            q_out.header['BUNIT'] = units
            dq_out.header['BUNIT'] = '(' + units + ')' + '**2'
            u_out.data *= convfact
            du_out.data *= convfact**2
            u_out.header['BUNIT'] = units
            du_out.header['BUNIT'] = '(' + units + ')' + '**2'
            unit_string = 'Jy-arcs2'
        elif units.upper() == 'JY/BEAM':
            unit_string = 'Jy-beam'
        else:
            print()
            print('ERROR: requested output unit {} not available. Please choose "Jy/arcsec2" or "Jy/beam".'.format(units))
            print()
        
        # Make & (if requested) save HDUs; each Stokes parameter with its variance
        imap = fits.HDUList([i_out, di_out])
        qmap = fits.HDUList([q_out, dq_out])
        umap = fits.HDUList([u_out, du_out])
        if savefile != '':
            imap.writeto(savefile + '_i+var_' + unit_string + '.fits', overwrite=True)
            qmap.writeto(savefile + '_q+var_' + unit_string + '.fits', overwrite=True)
            umap.writeto(savefile + '_u+var_' + unit_string + '.fits', overwrite=True)
        else:
            print('Files not saved (specify filename as SAVEFILE keyword if you wanted to)')
            
        # Output
        return imap, qmap, umap

    elif instru == 'HERSCHEL' or (instru == 'PACS' or instru == 'SPIRE'):
        imap_raw = fits.open(fname)
        imap = fits.PrimaryHDU()
        imap.data = imap_raw['PRIMARY'].data     # This automatically creates a header
        imap.header = imap_raw['PRIMARY'].header
        imap.data *= 1e6                         # Conversion MJy/sr --> Jy/arcsec**2
        imap.data /= 3600**2*(180./np.pi)**2
        imap.header['BUNIT'] = 'Jy/arcsec2'
        # Output
        return imap
    
    else:
        print()
        print('ERROR: Instrument {} not recognized. Please use HAWC+/SOFIA, JCMT/SCUBA2/POL2 or HERSCHEL/PACS/SPIRE.'.format(instru))
        print()



### READ FILES OF STARLINK-SMOOTHED MAPS ###

def readfilesfinal(fnames, units_in = 'Jy/beam', units_out = 'Jy/arcsec2', beam = 0.):
    # Read Starlink-smoothed files, output I, Q, U maps in HDU format
    # Converts map units to unit of choice: Jy/arcsec^2 (default), Jy/beam
    # Beam (FWHM in arcsec) needed for JCMT conversion; see Dempsey et al. 2013

    units_in = units_in.upper()
    units_out = units_out.upper()

    i_content = fits.open(fnames[0])
    q_content = fits.open(fnames[1])
    u_content = fits.open(fnames[2])
    
    if units_out == units_in:
        convfact = 1.
    elif (units_in == 'JY/BEAM' and units_out == 'JY/ARCSEC2'):
        if beam <= 0:                                          # Beam FWHM in arcsec
            beam = float(i_content['PRIMARY'].header['BEAM'])
        beam_sig = beam / 2*np.sqrt(2*np.log(2))               # FWHM --> sigma (assuming Gaussian beam)
        convfact = 2 * np.pi * beam_sig**2                     # Beam area in arcsec^2
    elif (units_in == 'JY/ARCSEC2' and units_out == 'JY/BEAM'):
        if beam <= 0:                                          # Beam FWHM in arcsec
            beam = float(i_content['PRIMARY'].header['BEAM'])
        beam_sig = beam / 2*np.sqrt(2*np.log(2))               # FWHM --> sigma (assuming Gaussian beam)
        convfact = 1./(2 * np.pi * beam_sig**2)                # Inverse of beam area in arcsec^2
    else:
        #print('Either the UNITS_IN value "%s" or the UNITS_OUT variable "%s" is unrecognized.' (units_in, units_out))
        print('Either the UNITS_IN or the UNITS_OUT variable is unrecognized.')
        print('Please choose them from the following selection: "Jy/arcsec2", "Jy/beam".')

    #pdb.set_trace()

    ## I, dI
    i = fits.PrimaryHDU()
    i.data = i_content['PRIMARY'].data / convfact              # This automatically creates a minimal header (BITPIX, NAXIS, NAXIS1/2, EXTEND)
    i.header['CDELT1'] = i_content['PRIMARY'].header['CDELT1'] # Completing the header
    i.header['CDELT2'] = i_content['PRIMARY'].header['CDELT2']
    i.header['CRVAL1'] = i_content['PRIMARY'].header['CRVAL1']
    i.header['CRVAL2'] = i_content['PRIMARY'].header['CRVAL2']
    i.header['CRPIX1'] = i_content['PRIMARY'].header['CRPIX1']
    i.header['CRPIX2'] = i_content['PRIMARY'].header['CRPIX2']
    i.header['CTYPE1'] = i_content['PRIMARY'].header['CTYPE1']
    i.header['CTYPE2'] = i_content['PRIMARY'].header['CTYPE2']
    i.header['EQUINOX'] = i_content['PRIMARY'].header['EQUINOX']
    i.header['EXTNAME'] = 'DATA    '
    i.header['BUNIT'] = 'Jy/arcsec2'
    #i.header['PCOUNT'] = i_content['PRIMARY'].header['PCOUNT'] # These are present in HAWC+ but not JCMT. What are they? Are they important?
    #i.header['GCOUNT'] = i_content['PRIMARY'].header['GCOUNT']
    #i.header['NWHP'] = i_content['PRIMARY'].header['NWHP']
    di = fits.PrimaryHDU()
    di.data = np.sqrt(i_content['VARIANCE'].data) / convfact
    di.header = copy.deepcopy(i.header)
    di.header['EXTNAME'] = 'ERROR   '
    di.header['BUNIT'] = 'Jy/arcsec2'
    ## Q, dQ
    q = fits.PrimaryHDU()
    q.data = q_content['PRIMARY'].data / convfact
    q.header = copy.deepcopy(i.header)
    dq = fits.PrimaryHDU()
    dq.data = np.sqrt(q_content['VARIANCE'].data) / convfact
    dq.header = copy.deepcopy(i.header)
    dq.header['EXTNAME'] = 'ERROR   '
    dq.header['BUNIT'] = 'Jy/arcsec2'
    ## U, dU
    u = fits.PrimaryHDU()
    u.data = u_content['PRIMARY'].data / convfact
    u.header = copy.deepcopy(i.header)
    du = fits.PrimaryHDU()
    du.data = np.sqrt(u_content['VARIANCE'].data) / convfact
    du.header = copy.deepcopy(i.header)
    du.header['EXTNAME'] = 'ERROR   '
    du.header['BUNIT'] = 'Jy/arcsec2'
    
    return i, q, u, di, dq, du
    


### MAKE PI, P,  THETA MAPS FROM I, Q, U ###

def make_p_theta(i, q, u, di, dq, du, save = '', magfield = 'True'):
    # Obtain polarized quantities (I_p, p, theta) using formulae from 30 Dor tutorial
    #  We assume for now that error cross-correlations are 0

    # POLARIZED INTENSITY
    Ip_raw = fits.PrimaryHDU()
    Ip = fits.PrimaryHDU()
    dIp = fits.PrimaryHDU()
    Ip_raw.data = np.sqrt(q.data**2 + u.data**2)
    dIp.data = np.sqrt((q.data * dq.data)**2 + (u.data * du.data)**2) / Ip_raw.data
    Ip.data = np.sqrt(Ip_raw.data**2 - dIp.data**2)
    Ip_raw.header = i.header
    dIp.header = i.header
    Ip.header = i.header
    
    ## POL. FRACTION (%)
    pfrac_raw = fits.PrimaryHDU()
    pfrac = fits.PrimaryHDU()
    dpfrac = fits.PrimaryHDU()
    #pdb.set_trace()
    pfrac_raw.data = 100 * np.sqrt((q.data/i.data)**2 + (u.data/i.data)**2)
    dpfrac.data = 100 * np.sqrt(
        ((q.data * dq.data)**2 + (u.data * du.data)**2) / ((q.data)**2 + (u.data)**2) +
    ((q.data / i.data)**2 + (u.data * i.data)**2) * di.data**2
    ) / i.data
    pfrac.data = np.sqrt(pfrac_raw.data**2 - dpfrac.data**2)
    pfrac_raw.header = i.header
    dpfrac.header = i.header
    pfrac.header = i.header
    
    ## POL. ANGLES
    theta = fits.PrimaryHDU()
    dtheta = fits.PrimaryHDU()
    theta.data = (90./np.pi) * np.arctan2(u.data, q.data) #+ 90.
    if magfield == 'True':
        theta.data = theta.data + 90.
    dtheta.data = (90./np.pi) * np.sqrt((q.data * dq.data)**2 + (u.data * du.data)**2) / (q.data**2 + u.data**2)
    theta.header = i.header
    dtheta.header = i.header

    if save != '':
        filename = save + '_allvars_smoothed+regridded.fits'
        i_map = fits.PrimaryHDU()
        i_map.data = i.data
        i_map.header = i.header
        di_map = fits.ImageHDU()
        di_map.data = di.data
        #di_map.header = di.header
        q_map = fits.ImageHDU()
        q_map.data = q.data
        #q_map.header = q.header
        dq_map = fits.ImageHDU()
        dq_map.data = dq.data
        #dq_map.header = dq.header
        u_map = fits.ImageHDU()
        u_map.data = u.data
        #u_map.header = u.header
        du_map = fits.ImageHDU()
        du_map.data = du.data
        #du_map.header = du.header
        Ip_raw_map = fits.ImageHDU()
        Ip_raw_map.data = Ip_raw.data
        #Ip_raw_map.header = Ip_raw.header
        Ip_map = fits.ImageHDU()
        Ip_map.data = Ip.data
        #Ip_map.header = Ip.header
        dIp_map = fits.ImageHDU()
        dIp_map.data = dIp.data
        #dIp_map.header = dIp.header
        pfrac_raw_map = fits.ImageHDU()
        pfrac_raw_map.data = pfrac_raw.data
        #pfrac_raw_map.header = pfrac_raw.header
        pfrac_map = fits.ImageHDU()
        pfrac_map.data = pfrac.data
        #pfrac_map.header = pfrac.header
        dpfrac_map = fits.ImageHDU()
        dpfrac_map.data = dpfrac.data
        #dpfrac_map.header = dpfrac.header
        theta_map = fits.ImageHDU()
        theta_map.data = theta.data
        #theta_map.header = theta.header
        dtheta_map = fits.ImageHDU()
        dtheta_map.data = dtheta.data
        #dtheta_map.header = dtheta.header
        allvars = fits.HDUList([i_map, di_map, q_map, dq_map, u_map, du_map, Ip_raw_map, Ip_map, dIp_map, pfrac_raw_map, pfrac_map, dpfrac_map,
                                theta_map, dtheta_map])
        allvars.writeto(filename, overwrite=True)
    
    return Ip, pfrac, theta, dIp, dpfrac, dtheta



### PLOT CONTOUR LEVEL PICTURES (FOR S/N ASSESSMENT) ###

def plot_levels(array, darray = None, title = None, contours = 'SNR', contour_values = [3., 10., 20.],
                colors = ['darkgrey', 'dimgrey', 'black'], crop = -1., showbeam = False, cmap = 'rainbow'):

    if darray == None and contours == 'SNR':
        print('ERROR (plot_levels): You need a noise array to do anything SNR-related with this function')
        
    if contours == 'SNR':
        SNR = fits.ImageHDU()
        SNR.data = array.data/darray.data
        SNR.header = array.header
        map4levels = SNR
    else:
        map4levels = array

    #if plotSNR:
    #    arr2plot = SNR #copy.deepcopy(SNR)
    #else:
    #    arr2plot = array #copy.deepcopy(array)

    # Crop figure if required
    if crop > 0.:
        xlen, ylen = array.data.shape
        xarr, yarr = np.ogrid[0:xlen, 0:ylen]
        cropmask = ((xarr - xlen / 2) * array.header['CDELT1']*60) ** 2 + ((yarr - ylen / 2) * array.header['CDELT2']*60) ** 2 > crop**2
        array.data[cropmask] = np.nan

    # Plot(s)
    fig = FITSFigure(array, subplot=(1,1,1))
    fig.show_colorscale(cmap=cmap)
    fig.set_title(title)
    fig.add_colorbar()
    fig.colorbar.set_axis_label_text('Flux (Jy/arcsec$^2$)')
    fig.show_contour(map4levels, colors=colors, levels=contour_values, smooth=1, kernel='box', linewidths=1.5)
    #if contours == 'SNR':
    #    fig.show_contour(SNR, colors=colors, levels=contour_values, smooth=1, kernel='box', linewidths=1.5)
    #else:
    #    fig.show_contour(array, colors=colors, levels=contour_values, smooth=1, kernel='box', linewidths=1.5)
    '''
    # Older version
    fig = FITSFigure(array, subplot=(1,1,1))
    # Show intensity
    fig.show_colorscale(cmap=cmap)
    fig.set_title(title)
    # Add colorbar
    fig.add_colorbar()
    if plotSNR:
        fig.colorbar.set_axis_label_text('S/N')
    else:
        fig.colorbar.set_axis_label_text('Flux (Jy/arcsec$^2$)')
    # Levels contours
    fig.show_contour(SNR, colors='white', levels=contour_values, smooth=1, kernel='box', linewidths=1.6)
    fig.show_contour(SNR, colors=colors, levels=contour_values, smooth=1, kernel='box', linewidths=1.)
    '''
    
    return fig



### PLOT INTENSITY MAP + POLARIZATION VECTORS ###

def plot_I_polvec(i_in, di_in, p_in, dp_in, theta_in, dtheta_in, mask = None, title = None, subplot = (1, 1, 1), i_SNR_cut = -1., p_SNR_cut = -1.,
                     crop = -1., thresh = -1., scalevec = .6, stepvec = 1, showbeam = False, cmap = 'rainbow'):

    #SNRmask = 0.
    i = copy.deepcopy(i_in)
    di = copy.deepcopy(di_in)
    p = copy.deepcopy(p_in)
    dp = copy.deepcopy(dp_in)
    theta = copy.deepcopy(theta_in)
    dtheta = copy.deepcopy(dtheta_in)
    #pdb.set_trace()

    #if 'mask' in locals():
    if hasattr(mask, "__len__"):
        mask_final = mask
        if crop > 0.:
            print('Variable ''crop'' not needed when ''mask'' is defined. Ignoring.')
        if i_SNR_cut > 0.:
            print('Variable ''i_SNR_cut'' not needed when ''mask'' is defined. Ignoring.')
        if p_SNR_cut > 0.:
            print('Variable ''p_SNR_cut'' not needed when ''mask'' is defined. Ignoring.')
        if thresh > -1.:
            print('Variable ''thresh'' not needed when ''mask'' is defined. Ignoring.')
    else:
        # Crop if required
        if crop >= 0.:
            xlen, ylen = i.data.shape
            xarr, yarr = np.ogrid[0:xlen, 0:ylen]
            cropmask = ((xarr - xlen / 2) * i.header['CDELT1']*60) ** 2 + ((yarr - ylen / 2) * i.header['CDELT2']*60) ** 2 > crop**2
        else:
            xlen, ylen = i.data.shape
            xarr, yarr = np.ogrid[0:xlen, 0:ylen]
            cropmask = ((xarr - xlen / 2) * i.header['CDELT1']*60) ** 2 + ((yarr - ylen / 2) * i.header['CDELT2']*60) ** 2 < 0.
        # Select by S/N if required:
        if i_SNR_cut >= 0. and p_SNR_cut >= 0.:
            #pdb.set_trace()
            SNRmask = (i.data/di.data < i_SNR_cut) | (p.data/dp.data < p_SNR_cut)
        elif i_SNR_cut >= 0.:
            #pdb.set_trace()
            SNRmask = i.data/di.data < i_SNR_cut
        elif p_SNR_cut >= 0.:
            #pdb.set_trace()
            SNRmask = p.data/dp.data < p_SNR_cut
        else:
            SNRmask = i.data/di.data < -1.
        # Select by I threshold if required:
        thresh_mask = i.data < thresh
        # Putting all together
        mask_final = np.any([cropmask, SNRmask, thresh_mask], axis = 0)
        
    #i.data[mask_final] = np.nan
    p.data[mask_final] = np.nan
    theta.data[mask_final] = np.nan

    # Show image
    fig = FITSFigure(i, subplot = subplot)#(p.data/dp.data, subplot = subplot)#
    fig.show_colorscale(cmap = cmap)
    fig.set_title(title)
    # Add colorbar
    fig.add_colorbar()
    fig.colorbar.set_axis_label_text('Flux (Jy/arcsec$^2$)')
    # Add vectors
    fig.show_vectors(p, theta, scale = scalevec, step = stepvec)
    pxscale = i.header['CDELT2']*3600
    vectscale = scalevec * pxscale / 3600
    fig.add_scalebar(20 * vectscale, "p = 20%",corner = 'top right', frame = True)
    # Add beam indicator if required
    if showbeam:
        fig.add_beam(facecolor = 'red', edgecolor = 'black', linewidth = 2, pad = 1, corner = 'bottom left')
        fig.add_label(0.02, 0.02, 'Beam FWHM', horizontalalignment = 'left', weight = 'bold', relative = True, size = 'small')
    # Format
    fig.tick_labels.set_font(size='large')
    fig.axis_labels.set_font(size='large')

    '''
    if 'cropmask' in locals() and 'SNRmask' in locals():
        return cropmask, SNRmask, fig
    elif 'cropmask' in locals():
        return cropmask, fig
    elif 'SNRmask' in locals():
        return SNRmask, fig
    else:
        return fig
    '''
    
    return fig



### LINEAR FUNCTION FIT WITH 2D ERRORS ###
import corner

def lnprior(theta):
    m,b=theta
    if 0. < b < 1. and -1. < m < 0.:
        return 0.0
    return -np.inf

def lnlike(theta,x,y,xerr,yerr):
    m,b=theta
    model=m*x + b   #no longer necessary, but I've left it in so you can still see what it is.
    angle=np.arctan(m)
    delta=-1.*np.sin(angle)*x + np.cos(angle)*y - b*np.cos(angle)
    sigmasq=xerr*np.sin(angle)**2 + yerr*np.cos(angle)**2
    return -np.sum(0.5* delta**2 / sigmasq)#0.5 * (y - model)**2/yerr**2 + np.log(2*np.pi*yerr**2)

def lnprob(theta,x,y,xerr,yerr):
    lp=lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, xerr, yerr)

def run_emcee(sampler,pos,ndim,labels,steps=500,prefix=""):
    print("Running MCMC...")
    sampler.run_mcmc(pos,steps, rstate0=np.random.get_state())
    print("Done.")

    plt.clf()    
    fig, axes = plt.subplots(ndim, 1, sharex=True, figsize=(8, 9))
    
    for i in range(ndim):
        axes[i].plot(sampler.chain[:, :, i].T, color="k", alpha=0.4)
        axes[i].set_ylabel(labels[i])

    fig.tight_layout(h_pad=0.0)
    fig.savefig(prefix+"line-time.png")
    return sampler

def mcmc_results(sampler,ndim,percentiles=[16, 50, 84],burnin=200,labels="",prefix=""):

    samples = sampler.chain[:, burnin:, :].reshape((-1, ndim))
    print(samples.shape)
    
    fig = corner.corner(samples, color='blue',labels=labels[0:ndim],quantiles=[0.16, 0.5, 0.84],show_titles=True,cmap='blues')
    fig.savefig(prefix+"line-triangle.png")
    credible_interval=[]
    for i in range(ndim):
        credible_interval.append(np.percentile(samples[:,i], percentiles))
        credible_interval[i][2] -= credible_interval[i][1]
        credible_interval[i][0] = credible_interval[i][1] - credible_interval[i][0]
    print("MCMC results:")
    for i in range(ndim):
        print("{0}  = {1[1]} + {1[2]} - {1[0]}".format(labels[i],credible_interval[i]))

    #now produce output plots of the distribution of lines
    fig=plt.figure()
    ax=fig.add_subplot(111)
    xplot=np.arange(-1,7,0.3)
    try:
        for m,b in samples[np.random.randint(len(samples),size=1000),0:2]:
            ax.plot(xplot,m*xplot+b,color="k",alpha=0.02)
    except:
        for m,b in samples[np.random.randint(len(samples),size=1000)]:
            ax.plot(xplot,m*xplot+b,color="k",alpha=0.02)
    ax.errorbar(x,y,xerr=xerr,yerr=yerr,fmt="ob")
    ax.set_xlim([-1,2])
    ax.set_ylim([-1,1])
    fig.savefig(prefix+"line-mcmc.png")


### LIMIT ANGLE TO 0-MAX RANGE ###
def wrapangle(angles, limit = 180.):
    """ Wraps an angle, or a 2D array of angles, to between two limiting values."""
    """ 'limit' can be a scalar, in which case the range chosen is 0 to limit, or a 1D array, in which case  the range is first to last 'limit' value."""
    if np.isscalar(limit):
        limit1 = 0
        limit2 = limit
    else:
        limit1 = limit[0]
        limit2 = limit[len(limit) - 1]
    span = limit2 - limit1
    
    if np.isscalar(angles):
        while(angles >= limit2):
            angles += -span
        while(angles < limit1):
            angles += span
        angles_out = angles
    else:
        if angles.ndim==1:
            for i, angle in enumerate(angles):
                while(angle >= limit2):
                    angle += -span
                while(angle < limit1):
                    angle += span
                angles[i] = angle
        elif angles.ndim==2:
            for i, row in enumerate(angles):
                for j, angle in enumerate(row):
                    while(angle >= limit2):
                        angle += -span
                        while(angle < limit1):
                            angle += span
                        angles[i,j] = angle
    return angles


### ANGLE DIFFERENCES ###
### Adapted from angle wrapping function
def polangle_diff(angles1, angles2, liminf = 0., limsup = 180.):
    """ Calculates the difference between *polarization* angles, or 2D arrays of angles (i.e. wrapping 0-180)."""
    angles = angles1 - angles2
    span = limsup - liminf
    if np.isscalar(angles):
        while(angles >= limsup):
            angles += -span
        while(angles < liminf):
            angles += span
        angles_out = angles
    else:
        for i, row in enumerate(angles):
            for j, angle in enumerate(row):
                while(angle >= limsup):
                    angle += -span
                while(angle < liminf):
                    angle += span
                angles[i,j] = angle
    return angles


### LINEAR FUNCTION FOR ODR FITTING ###

def linear_func(p, x):
   m, c = p
   return m*x + c

# Alternative:
def line(x, a, b):
    return a * x + b
