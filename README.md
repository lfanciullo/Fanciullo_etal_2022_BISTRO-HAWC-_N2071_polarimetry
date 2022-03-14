# Fanciullo_etal_2022_BISTRO-HAWCplus_N2071_polarimetry

Software used:
--------------

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


Available scripts:
------------------

* [mapmaker.py](mapmaker.py): Converts the data files from their original format to one that can be used by the other scripts.
  * Uses subroutines from [subroutines_mapmaking.py](subroutines_mapmaking.py).


Data recovery:
--------------

* JCMT observations of NGC 2071 can be found via the [JCMT Science Archive](https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/jcmt/) hosted at CADC. See also [the JCMT data access page](https://www.eaobservatory.org/jcmt/data-access/).
  * The data reduction procedure is described in [Lyo et al. 2021](https://iopscience.iop.org/article/10.3847/1538-4357/ac0ce9), Section 2.
  * In the present work, the maps have been further rebinned from a pixel size of 4'' to 8''. The Starlink commands used to do so are as follows:
  
    pol2map in=^infiles iout=iauto qout=! uout=! mapdir=maps qudir=qudata mapvar=no skyloop=no obsweight=no normalise=no pixsize=4
    pol2map in=qudata/\* iout=iext qout=qext uout=uext mapdir=maps mask=iauto maskout1=astmask maskout2=pcamask ipref=iext cat=mycat mapvar=true skyloop=true obsweight=no normalise=no pixsize=4
    pol2map in=^outfiles iout=! qout=! uout=! cat=mycat_8bin ipcor=false binsize=8
    
* SOFIA HAWC+ data
* Herschel data


Data formatting:
----------------


Data anlysis:
-------------


Plots:
------


