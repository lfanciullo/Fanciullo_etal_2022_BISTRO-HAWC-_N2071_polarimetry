
import numpy as np
from astropy.io import fits
from aplpy import FITSFigure
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
import copy
#import pdb

from varnames import *


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

def make_p_theta(i, q, u, di, dq, du, magfield = 'True'):
    # Obtain polarized quantities (I_p, p, theta) using formulae from 30 Dor tutorial
    #  Scripts assumes that error cross-correlations are 0

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
    
    return Ip, pfrac, theta, dIp, dpfrac, dtheta



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


### HERSCHEL MAPS CUTOUTS ###
def herschel_cutout(fname_in, fname_out = '', fname_ref = folder_reprocessed_data + files_850[0], size = 15.):
    """Script to make a cutout of the N2071 region from the full Orion B Herschel maps."""
    """Parameters used:"""
    """ fname_in: """
    """ fname_out: """
    """ fname_ref: """
    """ size: """
    refmap = fits.open(fname_ref)     # Reference map for central cutout coordinate
    size_deg = size/60.
    orion_b_data = fits.open(fname_in) # Orion B full file
    orion_b_map = orion_b_data[0]      # Orion B full file
    w = WCS(orion_b_map)               # WCS object
    position = SkyCoord(refmap[0].header['CRVAL1'], refmap[0].header['CRVAL2'], unit = 'deg') # Center coordinates for cutout 
    xsize = int(round(np.abs(size_deg / orion_b_map.header['CDELT1'])))                       # Cutout size in pixels
    ysize = int(round(np.abs(size_deg / orion_b_map.header['CDELT2'])))
    cutsize = (xsize, ysize)
    cutout = Cutout2D(orion_b_map.data, position, cutsize, wcs = w, mode = 'trim')            # Map cutout creation
    cuthdu = fits.PrimaryHDU()                                                                # HDU for storing cutout
    cuthdu.data = cutout.data
    cuthdu.header.update(cutout.wcs.to_header())                                              # Header creation for HDU
    if fname_out == '':
        fname_out = 'herschel_cutout.fits'
        print('Cutout map saved as {}. To save with a different filename use the FNAME_OUT parameter.'.format(fname_out))
    cuthdu.writeto(folder_reprocessed_data + fname_out, overwrite=True)
