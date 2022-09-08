# Fanciullo_etal_2022_BISTRO-HAWCplus_N2071_polarimetry

**Nota: The page is still being completed. If you do not find something you were looking for, check again in a few days.**

Software used
-------------

* Python. The scripts have been tested up to version 3.9.7, using IPython 7.25.0. Packages used:
  * aplpy
  * astropy
  * copy
  * matplotlib
  * numpy
  * scipy
  * statsmodels
* Starlink, a data reduction and analysis suite for JCMT data (links: [download page](http://starlink.eao.hawaii.edu/starlink), [introduction](https://www.eaobservatory.org/jcmt/observing/getting-started/#Starlink_analysis_and_reduction_software), [cheat sheet](https://www.eaobservatory.org//jcmt/wp-content/uploads/sites/2/2016/04/StarlinkBeginner.pdf)). The scripts on this page have been developed using version 2018A. Packages used:
  * convert
  * kappa
  * smurf


Available scripts and files
---------------------------

### Folders ###
The scripts, as they are provided, use the following folders to store data:
* The N2071 maps downloaded in their original form (see ['Online data'](#Online-data)) go into *data_downloaded*.
* The maps after conversion into the final format (both with and without smoothing) go into *data_reprocessed*.
* The output figures are saved in the main folder.

The user can modify the folders and/or paths for the data; if so, however, the pathnames in [varnames.py](varnames.py) should be updated.


### Files and scripts ###
* [varnames.py](varnames.py) defines all the variables used in Python scripts. Any modifications to the standard variables -- e.g. to modify the file paths -- needs to be done here.
* [Starlink_commands_HAWC+JCMT_sample.txt](Starlink_commands_HAWC+JCMT_sample.txt) and [Starlink_commands_Herschel_sample.txt](Starlink_commands_Herschel_sample.txt) contain examples of the *Starlink* commands used to resampla nad smooth the maps (see also ['Data formatting'](#Data-formatting)).
  * [tab](tab) and [tab2](tab2): format tables needed to smooth the maps with *Starlink*.
* [all_functions.py](all_functions.py) contains all the functions for reformatting the data used in other Python scripts.  
* [data_analysis.py](data_analysis.py): Script to perform the analysis of the N2071  data and create the plots shown in the article.
  * Uses functions from [all_functions.py](all_functions.py).


Online data
-----------

* JCMT observations of NGC 2071 can be found via the [JCMT Science Archive](https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/jcmt/) hosted at CADC (See also [the JCMT data access page](https://www.eaobservatory.org/jcmt/data-access/)). The following instrutions assume you have dowmloaded the raw data from project code M16AL004, object name: NGC2071IR, instrument: SCUBA-2 POL-2.
  * The *Starlink* commands to produce maps binned to a pixel size of 8'' (as opposed to the default 4'') are shown below. See also Section 2 of [Lyo et al. 2021](https://iopscience.iop.org/article/10.3847/1538-4357/ac0ce9) for a description of the data reduction procedure.
  
    > pol2map in=^infiles iout=iauto qout=! uout=! mapdir=maps qudir=qudata mapvar=no skyloop=no obsweight=no normalise=no pixsize=4
    >
    > pol2map in=qudata/\* iout=iext qout=qext uout=uext mapdir=maps mask=iauto maskout1=astmask maskout2=pcamask ipref=iext cat=mycat mapvar=true skyloop=true obsweight=no normalise=no pixsize=4
    >
    > pol2map in=^outfiles iout=! qout=! uout=! cat=mycat_8bin ipcor=false binsize=8
    
* The SOFIA (HAWC+) observations can be found at the [the SOFIA IRSA database](https://irsa.ipac.caltech.edu/applications/sofia/). The relevant data is in AORs 07_130_07 and 07_130_08 (processing level 4). Assuming the two files are downloaded together, the path to the .fits files is
  > (download folder name)/data/SOFIA/HAWC_PLUS/L4/p10196/data/g5/F0615_HA_POL_0701308_HAWDHWPD_PMP_071-114.fits
  
  for band D, and

  > (download folder name)/data/SOFIA/HAWC_PLUS/L4/p10197/data/g1/F0616_HA_POL_0701309_HAWEHWPE_PMP_067-102.fits'

  for band E, respectively. 
* The Herschel data can be downloaded from [the Herschel Gould Belt Survey Archive](http://gouldbelt-herschel.cea.fr/archives). This work used the temperature map and and the high-resolution column density map for Orion B:
   * HGBS_orionB_dust_temperature_map.fits
   * HGBS_orionB_hires_column_density_map.fits

   Since these maps cover a much larger area than is needed for the present work, we suggest to create a fits file of the N2071 area alone (see ['Data formatting'](#Data-formatting) below).


Data formatting
---------------

The data files recovered in [the previous section](#Online-data) all have different units and formats. This section contains the instructions to fix their format and make two sets of data:
1. Maps with a common format and units, retaining their original pixel and beam size ('original pixel' data);
2. Maps resampled and smoothed to a common pixel and beam size (8'' and 18''.9, respectively), for direct inter-band comparison.
By default all these are saved in the *data_reprocessed/* folder, (the value of *folder_reprocessed_data* in [varnames.py](varnames.py)).

### HAWC+ and JCMT data ###
The 'original pixel' maps can be created in Python running the following commands in the main repository:
~~~python
from all_functions import convertfiles
from varnames import *

filed = folder_download + file_d_download  # Read HAWC+ D band map
imap_d, qmap_d, umap_d = convertfiles(filed, instru = 'HAWC+', units = 'Jy/arcsec2', savefile = folder_reprocessed_data + 'N2071_HAWC+D')  # Convert and save

filee = folder_download + file_e_download  # Read HAWC+ E band map
imap_e, qmap_e, umap_e = convertfiles(filee, instru = 'HAWC+', units = 'Jy/arcsec2', savefile = folder_reprocessed_data + 'N2071_HAWC+E')

files_850 = [folder_download + s for s in files_850_download]  # Read JCMT 850 um map (3 separate files for I, Q and U)
imap_850, qmap_850, umap_850 = convertfiles(files_850, instru = 'SCUBA2', units = 'Jy/arcsec2', beam = 14.1, savefile = folder_reprocessed_data + 'N2071_JCMT-850-8as')
~~~
The 'savefile' strings used in this example are chosen so that the output filenames will be the same as the default *files_(band)_originalpixel* variables in [varnames.py](varnames.py).

To resample and smooth the data, run the *Starlink* commands in [Starlink_commands_HAWC+JCMT_sample.txt](Starlink_commands_HAWC+JCMT_sample.txt) from inside the *data_reprocessed/* folder. Running the *wcsalign* commands will return a comment similar to the following, which will require to press Enter to proceed:
  > LBND - Lower pixel index bounds for output images /[-64,-64]/ > 
  > UBND - Upper pixel index bounds for output images /[69,68]/ >
  >
  > The 1st NDF (filename) does not have the selected QUALITY component.
  > The output FITS file will be produced for the remaining array components that were specified.

The values used in [Starlink_commands_HAWC+JCMT_sample.txt](Starlink_commands_HAWC+JCMT_sample.txt) for the *gausmooth* 'fwhm' parameter are the ones needed to smooth all bands to the same resolution as the HAWC+E map (18''.9). The values depend on the original resolution and pixel size of the maps; if you want to change the maps' final resolution, modify the fwhm value according to section 2.3.1 of Fanciullo et al. 2022.


### *Herschel* data ###

To create a coutout of the N2071 region out of the whole Orion B map, one can run the following Python commands from the main repository:
~~~python
from all_functions import herschel_cutout
from varnames import *

herschel_cutout(folder_download + file_Herschel_h2_large, fname_out = file_Herschel_h2_cutout, size = 20.)
herschel_cutout(folder_download + file_Herschel_T_large, fname_out = file_Herschel_T_cutout, size = 20.)
~~~

The two cutouts thus created -- one for the column density, one for dust temperature -- are saved by default in the *data_reprocessed/* folder. We recommend to use a value of 20 or higher for the 'size' parameter (cutout edge length, in arcmin). 

To make a *P* vs. &tau; plot, as in section 3.4 of Fanciullo et al. 2022, the cutout maps need to be resampled to the pixel size of the HAWC+ and JCMT bands. This can be done by running the *Starlink* commands in [Starlink_commands_Herschel_sample.txt](Starlink_commands_Herschel_sample.txt) from inside the *data_reprocessed/* (or equivalent) folder. Note that these commands assume the default filenames from [varnames.py](varnames.py) are used; change the commands as appropriate if you are not using the default names.


Data analysis
-------------

Running lines 1-442 of [data_analysis.py](data_analysis.py) reads the reformatted IQU files (for HAWC+ and JCMT), uses them calculate the polarized quantities of interest (polarization fraction and angle in all bands, difference in polarization agle across bads) and creates masks to select the data by S/N and by position.

This provides the user with the quantities needed to recreate the analysis in Fanciullo et al. 2022.

<!--- Works in progress:
1. Herschel T and NH plots
...
--->


Plots
-----

Running lines 454-1704 of [data_analysis.py](data_analysis.py) recreates the figures used in Fanciullo et al. 2022.


