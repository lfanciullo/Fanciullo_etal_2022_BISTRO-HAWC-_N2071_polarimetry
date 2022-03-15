# Fanciullo_etal_2022_BISTRO-HAWCplus_N2071_polarimetry

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


Available scripts
-----------------

* [mapmaker.py](mapmaker.py): Converts the data files from their original format to one that can be used by the other scripts.
  * Uses subroutines from [subroutines_mapmaking.py](subroutines_mapmaking.py).


Online data
-----------

* JCMT observations of NGC 2071 can be found via the [JCMT Science Archive](https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/jcmt/) hosted at CADC. See also [the JCMT data access page](https://www.eaobservatory.org/jcmt/data-access/).
  * The data reduction procedure is described in [Lyo et al. 2021](https://iopscience.iop.org/article/10.3847/1538-4357/ac0ce9), Section 2.
  * In the present work, the maps have been further rebinned from a pixel size of 4'' to 8''. The Starlink commands used to do so are as follows:
  
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

The data files recovered in [the previous section](#Online-data) all have different units and formats. This section contains the instructions to fix their format to make sets of data:
1. Maps with a common format and units, but each retaining their original ixel and beam size ('orignal pixel' data);
2. Maps resampled and smoothed to a common pixel and beam size (8'' and 18''.9, respectively), for direct inter-band comparison. 


Data analysis
-------------


Plots
-----


