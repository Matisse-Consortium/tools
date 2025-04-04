# [MATISSE python tools](https://github.com/Matisse-Consortium/tools/wiki)

## You are on the MATISSE python tools distribution page.
The [MATISSE](http://www.eso.org/sci/facilities/paranal/instruments/matisse.html) python tools provide to the general MATISSE user an environment for running in a more friendly way the [standard data reduction software](http://www.eso.org/sci/software/pipelines/matisse/matisse-pipe-recipes.html) that is distributed separately by ESO, as well as visualizing intermediate and final products of the pipeline. These scripts were originally developed for the MATISSE commissionning. 

## Python tools

### Possible installation:
[`pip install git+https://github.com/Matisse-Consortium/tools.git`]:#
- (if necessary) install gtk+-3.0 on your system (sudo apt-get install gtk+-3.0)
- (if necessary) install the following Python packages: attrdict tqdm astroquery openpyxl statistics ObjectListView numpy matplotlib shapely wxpython skycalc_cli
- download the tools ZIP file by clicking on the green button named 'Code'
- Unzip tools-master.zip on your system in the directory of your choice (for instance $HOME/your_dir)
- add the following path $HOME/your_dir/tools-master/mat_tools in your PATH environment variable.

### Important note
The Python scripts of the mat_tools package were designed to be run in a terminal (batch mode). No execution in an interactive console mode is offered for the moment.

### Step-by-step description of a MATISSE data reduction with the Python tools: [here](https://github.com/Matisse-Consortium/tools/wiki/Getting-Started) 

### Quick use (with list of main Python routines)
1. `mat_autoPipeline.py <raw data directory>`: runs the MATISSE pipeline to obtain reduced oifits files.  
2. `mat_tidyupOiFits.py <reduced data directoy>`: renames the reduced oifits files in a more explicit way and put them in a new ***_OIFITS directory.
3. `mat_showTransFunc.py <oifits directory>`: shows the raw visibilities of the science targets and the instrumental visibilities of the calibrators over the considered MATISSE observing sequences.
3. `mat_autoCalib.py <oifits directoy>`: calibrates the interferometric observables and store them in a new ***_CALIBRATED.
4. `mat_fluxCal.py <oifits directory>`: performs (if needed) (correlated or total) flux calibration.
5. `mat_mergeAllOiFits.py <calibrated oifits directory>`: perfoms the final BCD calibration, i.e. merges the calibrated oifits files over the MATISSE observing cycle.

### More info:
- Installation instructions are given [there](https://github.com/Matisse-Consortium/tools/wiki/Installation).
- The development status of the pipeline (including the ESO version of the code) is provided [there](https://github.com/Matisse-Consortium/tools/wiki/Known-bugs-and-development-plan).

## IRBIS interface

We also distribute an interface to the IRBIS image reconstruction software (part of MATISSE DRS). Please refer to [this page](https://github.com/Matisse-Consortium/tools/tree/master/imarec) for more information on installation and use.
