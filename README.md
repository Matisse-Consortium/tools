# [MATISSE python tools](https://github.com/Matisse-Consortium/tools/wiki)

## You are on the MATISSE tools distribution page.
The MATISSE tools are an addition to the [standard data reduction software](http://www.eso.org/sci/software/pipelines/matisse/matisse-pipe-recipes.html) that is distributed separately by ESO. These scripts were developped for the MATISSE commissionning and may be of some use for general users to reduce and visualize data from [MATISSE](http://www.eso.org/sci/facilities/paranal/instruments/matisse.html). However, please note that they should be used at your own risks!

## Python tools

### Quick install:
`pip install git+https://github.com/Matisse-Consortium/tools.git`

### Quick use (with list of main Python routines)
1. `mat_autoPipeline.py <raw data directory>`: runs the MATISSE pipeline to obtain reduced oifits files.  
2. `mat_tidyupOiFits.py <reduced data directoy>`: renames the reduced oifits files in a more explicit way and put them in a new ***_OIFITS directory.
3. `mat_showTransFunc.py <oifits directory>`: shows the raw visibilities of the science targets and the instrumental visibilities of the calibrators over the considered MATISSE observing sequences.
3. `mat_autoCalib.py <oifits directoy>`: calibrates the interferometric observables and store them in a new ***_CALIBRATED.
4. `mat_fluxCal.py <oifits directory>`: performs (if needed) flux (correlated or total) flux calibration.
5. `mat_mergeAllOiFits.py <calibrated oifits directory>`: perfoms the final BCD calibration, i.e. merges the calibrated oifits files over the MATISSE observing cycle.

### More info:
- A step by step description of the reduction and and visualization of [MATISSE data](https://www.eso.org/public/news/eso1808/) with the Python routines, together with the MATISSE DRS, is given [here](https://github.com/Matisse-Consortium/tools/wiki/Getting-Started).
- Installation instructions are given [there](https://github.com/Matisse-Consortium/tools/wiki/Installation).
- The development status of the pipeline (including the ESO version of the code) is provided [there](https://github.com/Matisse-Consortium/tools/wiki/Known-bugs-and-development-plan).

## IRBIS interface

We also distribute an interface to the IRBIS image reconstruction software (part of MATISSE DRS). Please refer to [this page](https://github.com/Matisse-Consortium/tools/tree/master/imarec) for more information on installation and use.
