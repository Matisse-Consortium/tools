# [MATISSE python tools](https://github.com/Matisse-Consortium/tools/wiki)

## You are on the MATISSE tools distribution page.
The MATISSE tools are an addition to the [standard data reduction software](http://www.eso.org/sci/software/pipelines/matisse/matisse-pipe-recipes.html) that is distributed separately by ESO. These scripts were developped for the MATISSE commissionning and may be of some use for general users to reduce and visualize data from [MATISSE](http://www.eso.org/sci/facilities/paranal/instruments/matisse.html). However, please note that they should be used at your own risks!

## Python tools

### Quick install:
`pip install git+https://github.com/Matisse-Consortium/tools.git`

### Quick use:
1. `mat_autoPipeline.py <raw data directory>`
2. `mat_tidyupOiFits.py <product directoy>`
3. `mat_reflagData.py <oifits directory>`
3. `mat_autoCalib.py <reflagged directoy>`
4. `mat_showOiData.py <a MATISSE oifits file>`

### More info:
- A description of how to reduce and visualize [MATISSE data](https://www.eso.org/public/news/eso1808/) with these scripts together with the MATISSE DRS is given [here](https://github.com/Matisse-Consortium/tools/wiki/Getting-Started).
- Installation instructions are given [there](https://github.com/Matisse-Consortium/tools/wiki/Installation).
- The development status of the pipeline (including the ESO version of the code) is provided [there](https://github.com/Matisse-Consortium/tools/wiki/Known-bugs-and-development-plan).

## IRBIS interface

We also distribute an interface to the IRBIS image reconstruction software (part of MATISSE DRS). Please refer to [this page](https://github.com/Matisse-Consortium/tools/tree/master/imarec) for more information on installation and use.
