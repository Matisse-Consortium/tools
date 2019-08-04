# [MATISSE python tools](https://gitlab.oca.eu/MATISSE/tools/wikis/home)

## You are on the MATISSE python tools distribution page.
The MATISSE python tools are an addition to the [standard data reduction software](http://www.eso.org/sci/software/pipelines/matisse/matisse-pipe-recipes.html) that is distributed separately by ESO. These scripts were developped for the MATISSE commissionning and may be of some use for general users to reduce and visualize data from [MATISSE](http://www.eso.org/sci/facilities/paranal/instruments/matisse.html). However, please note that they should be used at your own risks!

### Quick install:
`pip install --user https://gitlab.oca.eu/MATISSE/tools/-/archive/master/tools-<version>.tar`

`<version>` can be `master` or `0.1`

### Quick use:
1. `mat_autoPipeline.py <raw data directory>`
2. `mat_tidyupOiFits.py <product directoy>`
3. `mat_autoCalib.py <oifits directoy>`
4. `mat_showOiData.py <a MATISSE oifits file>`

### More info:
- A description of how to reduce and visualize [MATISSE data](https://www.eso.org/public/news/eso1808/) with these scripts together with the MATISSE DRS is given [here](https://gitlab.oca.eu/MATISSE/tools/wikis/Using%20the%20pipeline).
- Installation instructions are given [there](https://gitlab.oca.eu/MATISSE/tools/wikis/Installation).
- The development status of the pipeline (including the ESO version of the code) is provided [there](https://gitlab.oca.eu/MATISSE/tools/wikis/Known%20bugs%20and%20development%20plan)
