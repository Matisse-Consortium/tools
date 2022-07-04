C-Shell script for easy usage of the MATISSE image reconstruction software
==========================================================================
(Karl-Heinz Hofmann, Dieter Schertl, Max-Planck-Institut for Radio Astronomy, Bonn, Germany, February 2018)
modified F. Millour, J. Drevon July 2022


1. Both C-shell scripts `mat_cal_imarec_all.csh` and `mat_cal_imarec_all.2.csh` are easy-to-use 
   callers of the MATISSE image reconstruction recipe mat_cal_imarec.

1. You need to have csh or tcsh running, i.e. if you are using bash (check that typing `echo $SHELL`).

1. On Ubuntu, try installing most of the necessary software by tying the following command: `sudo apt install tcsh ftools-fv ftools-pow wcstools gnuplot original-awk imagemagick gv libreadline-dev libncurses5-dev ncurses-dev curl libcurl4 libcurl4-gnutls-dev xorg-dev make gcc g++ gfortran perl-modules python3-dev python3-astropy python3-numpy python3-scipy psutils ps2eps`

2. All necessary scripts and files are stored in one folder, named `imarec`, anywhere on your computer.
The enviroment variable `$SCRIPTS` must point to this folder, i.e.
    `setenv SCRIPTS <your path>/imarec`
or
    `set SCRIPTS <your path>/imarec`
    `export SCRIPTS`
Put this command in the file `echo 'setenv SCRIPTS <your path>/imarec' >> .cshrc`

3. The following software is required:

   - The ESO pipeline for MATISSE must be installed. Follow the instructions there: https://www.eso.org/sci/software/pipelines/matisse/

   - The HEASARC ftools package
		- In ubuntu, install it with `sudo apt install ftools-fv ftools-pow`
   		- It can be downloaded from https://heasarc.gsfc.nasa.gov/ftools/ftools_menu.html 
     (for installation instructions see https://heasarc.gsfc.nasa.gov/lheasoft/install.html )
     As noted there, you must initialize the package by sourcing `$HEADAS/headas-init.sh`, i.e.
        `source $HEADAS/headas-init.csh`
     or 
        `$HEADAS/headas-init.sh`

   - The Image Utility Programs from WCSTools
		- In ubuntu, install it with `sudo apt install wcstools`
   		-  (see http://tdc-www.harvard.edu/software/wcstools/wcsprogsi.html),
     especially imrot (to rotate and reflect fits images).

   - Our scripts need the following LINUX programs, too:
     `gnuplot`, `awk`, `latex`, `convert`, `sort`, `psmerge`, `ps2eps`, `ftcopy, `gv`.
     You can check if all this software is available with:
        `$SCRIPTS/swtest.csh`
     
4. Overall reconstruction procedure.
   During the run, this software automatically uses different image reconstruction parameters to find the "best" reconstruction.
   These parameters are:
     - a number of regularisation functions $regFuncs
     - different regularisation hyperparameters $mu, which control the strength of the regularisation
     - a number of radii of a circular binary object mask, which resctricts the image space, where the reconstruction can have non-zero values

   You can use one of the two scripts `mat_cal_imarec_all.csh` or `mat_cal_imarec_all.2.csh`, their difference is explained below.
   Both scripts use a parameter file (`mat_cal_imarec_all.par` or `mat_cal_imarec_all.2.par`) where all parameters of the reconstruction run are described.

   Their overall reconstruction scheme is illustrated by the the following pseudocode: (`$name` denotes the name of the corresponding parameter in the parameter file)

	>     FOR EACH "regularization function" $regFuncs DO`
	>     	FOR EACH "start value of the hyperparameter" $muStarts DO
	>          # Call IRBis with the actual value of the regularization function 
	>          #and the hyperparameter
	>          # Internally IRBis then performes two loops:
	>            FOR EACH oradius = $oradiusStart+i*$stepSize (i = 0..$oradiusNumber-1) DO
	>               FOR EACH mu = $muStart*$muFactor^j (j = 0..$muNumber-1) DO
	>                  Reconstruct one Image with the current parameters oradius, 
	>                  #mu and regularization function
	>                  Calculate its qrec ("quality" parameter of the reconstruction,
	>										 explained below)
	>                  Manage an internal list of all reconstructions sorted by qrec 
	>               ENDFOR (mu)
	>            ENDFOR (oradius)
	>            RETURN reconstruction with the best (=smallest) qrec (explained below)
	>          # end of IRBis call
	>          Maintain an external list of all best reconstructions returned from IRBis,
	>          and continue to loop over $regFuncs and $muStarts
	>        ENDFOR ($muStarts)
	>      ENDFOR ($regFuncs)

	Look for the reconstruction with the smallest qrec in the external list and return this reconstruction as the best one.

   As a criterion to select the best reconstruction, the minimum of a value named qrec is used. It is calculated from the measured Chi-squares
   and the residual ratio values (for more details see [Hofmann, Weigelt, & Schertl, D. 2014, A&A, 565, A48](https://www.aanda.org/articles/aa/abs/2014/05/aa23234-13/aa23234-13.html)).
   There are two versions how qrec is calculated:
   * `$qrecmode=1`: it is calculated using the $\chi^2$ and the residual ratio values of both the visibilities and the phases (Closure or Fourier phase)
    * `$qrecmode=2`: it is calculated using the $\chi^2$ and the residual ratio values of only phases (Closure or Fourier phase)
--> experience shows, usually the qrec from $qrecmode=2 yields better results than qrec from $qrecmode=1! (but you should test $qrecmode=1, too)


5. The Difference between `mat_cal_imarec_all.csh` and `mat_cal_imarec_all.2.csh` is the different way, the $muStarts are defined and used.

   - `mat_cal_imarec_all.csh` loop over all $regFuncs using a fixed list of different start values $muStart.
     It will call IRBis with each of this $regFuncs/$muStart pairs.

   - `mat_cal_imarec_all.2.csh` loop over all $regFuncs, too, but using start values $muStart calculated as:
        $muStart(next) = $muStart(actual) * $muFactor0, starting with $muStart0 
     and calls IRBis with this $regFuncs/$muStart pairs.
     This loop stops, if the returned qrec starts to increase, i.e. the actual qrec is larger than the one of the run before.

   NOTE: the parameter $qrecmode does change the output of the script "mat_cal_imarec_all.2.csh", because the termination of the hyper parameter loop depends on
         the calculated qrec values which are different for $qrecmode=1 or 2;
	 the parameter $qrecmode does not change the output of the script "mat_cal_imarec_all.csh".

6. An image reconstruction run is started with: 
    $SCRIPTS/mat_cal_imarec_all.csh    parfile # (a template of the parameter file parfile is $SCRIPTS/mat_cal_imarec_all.par)
or
    $SCRIPTS/mat_cal_imarec_all.2.csh  parfile # (a template of the parameter file parfile is $SCRIPTS/mat_cal_imarec_all.2.par)


7. Tutorial -- Usage of the scripts `mat_cal_imarec_all.csh` and `mat_cal_imarec_all.2.csh`:

    * Step 0: Define the enviroment variable $SCRIPTS with the absolute path of the script folder, i.e.
    
        % setenv SCRIPTS /MYPATH/mat_cal_imarec.scripts00
    or
        $ set SCRIPTS /MYPATH/mat_cal_imarec.scripts00
        $ export SCRIPTS

	* Step 1: Create a work folder for the image reconstruction session (usually named as the target to be processed)
    and cd to this folder.

	* Step 2: Copy the parameter file `$SCRIPTS/mat_cal_imarec_all.par` (`$SCRIPTS/mat_cal_imarec_all.2.par`) of the C-shell script
	`mat_cal_imarec_all.csh` (`mat_cal_imarec_all.2.csh`) into your work folder
	    cp $SCRIPTS/mat_cal_imarec_all.par .

	* Step 3: Estimation of the size of the target and FOV to be used for image reconstruction:
		- The size of the target is estimated by fitting a circular Gaussian, uniform disk, fully darkened disk or a modified Lorentz
                  intensity distribution into the observed squared visibilities of the target and their corresponding reduced $\chi^2$.
        - Calculation of the proposed size for the FOV of the reconstruction (in pixel and milliarcseconds).

		1) Edit the parameter file `mat_cal_imarec_all.par`:

            set data       = (myfile1.oifits myfile2.oifits ...)
                                          # list all oifits data file you want to use for image reconstruction
            set lambdaList = (0.1 20)     # The wavelength intervall to be processed: use (0.1 20) to cover all data first, regardless of its wavelength;
			                              # in a next step you can select one (or more) wavelength intervall you want to use for image reconstrcution
            set guess      =  1           # insert the number 1 to switch on the calculations of the first action (estimating the quality of the data
		                                  # and image reconstruction parameters)

		2) Run the script:
		    $SCRIPTS/mat_cal_imarec_all.csh  mat_cal_imarec_all.par

		3) Have a look to the uv coverage of the data, to the wavelengths in the data, and to the fits of the geometrical models to the measured visibilities:
            		% gv Parameter.Estimation/uv.ps
            		% gv Parameter.Estimation/wavelengths.ps
		        % gv Parameter.Estimation/gaussudfdda.ps

		4) Select one or more wavelength intervalls you want to use for image reconstruction, 
        and change the parfile, i.e.
            set lambdaList = (3.5 3.8)

		5) Run the script again:
            $SCRIPTS/mat_cal_imarec_all.csh  mat_cal_imarec_all.par

		6) Check the selected wavelengths (Parameter.Estimation/wavelengths.ps) and the resulting uv coverage (Parameter.Estimation/uv.ps) of the data.
		The text file `Parameter.Estimation/data.parameter` contains information about the size of the target
        and the recommended FOV. The FOV for image reconstruction should be roughly ~2 to 4 x the target size!

        Be careful with the proposed sizes of the target. 
        If the Targets is composed of many distinct small objects far apart (i.e. double stars),
        then the estimation of the size is too small usually. 
        In this case, you should set the FOV to about 2-4 times the expected separation.

	* Step 4: First reconstruction run

        - Carefully read the hints for FOV and number of pixels in the Parameter.Estimation/data.parameter
		- Start with the smallest pixel number possible for the target, i.e. 64x64 pixels ($npix)
          	- The Image mask start radius should be about 0.5 to 2 times larger than the estimated target size ($oradiusStart)

          	1) edit the image reconstruction parameters into the parameter file `mat_cal_imarec_all.par` or `mat_cal_imarec_all.2.par`, respectively:
            (some values below are just examples)

             		set guess         = 0    Set to 0 to switch to the image reconstruction run
					set engine        = 2    Specify the optimization engine used; 1: ASA-CG, 2: L-BFGS-B. [1]
		        	set algoMode      = 1    Specify if bispectrum or complex visibilities will be used for reconstruction. 
                                              1 = use bispectrum only, 2 = use complex visibilities only, 3 = use bispectrum and complex visibilities. [1]
			        set costFunc      = 1   # Specify which chi square function Q[ok(x)] of the measured data should be minimized (ok(x) is the actual iterated image):
						# 1: Q[ok(x)] = Sum{ |ibis(u,v) - mbis(u,v)|^2/var(u,v) }, with ibis(u,v) & mbis(u,v) the iterated and measured bispectrum, repectively, and
						#    var(u,v) the variance of mbis(u,v)
						# 2: Q[ok(x)] = Sum{ |exp(i iph(u,v)) - exp(i mph(u,v))|^2/varph(u,v) + costWeight*|imod(u,v) - mmod(u,v)|^2/varmod(u,v) }, with 
						#    exp(i iph(u,v)) & exp(i mph(u,v)) the phasors of the iterated and measured bispectrum, repectively, and varph(u,v) the variance of exp(i mph(u,v)),
						#    and with imod(u,v) & mmod(u,v) the moduli of the iterated and measured bispectrum, repectively, and varmod(u,v) the variance of mmod(u,v),
						#    and costWeight is a weight for the modulus part of Q[ok(x)], [1.0]
					        # 3: Q[ok(x)] = Sum{ |exp(i iph(u,v)) - exp(i mph(u,v))|^2/varph(u,v) } + costWeight*Sum{ |ipow(u) - mpow(u)|^2/varpow(u) }, with
						#    ipow(u) & mpow(u) the iterated and measured power spectrum, repectively, and varpow(u) the variance of mpow(u), and
					    #    with costWeight is a weight for the power spectrum part of Q[ok(x)], [1.0]
             		set fov           = 50  # FOV (in mas) of the reconstruction 
                                                #  with respect to the recommendations in "Parameter.Estimation/data.parameter":
                                     	        #  use that FOV (or a smaller one) listed for different 2^n array sizes which best fits to the
                                     		#  expected size of the target (FOV ~ 2-4x target size, be careful with binaries as mentioned above!)
             		set npix          = 64  # Number of pixels corresponding to the FOV chosen above.
                                                #  For speed purpose, only powers of 2 should be used.
             		set oradiusStart  = 20  # First radius (in mas) of the binary image mask, which should be larger than the target.
                                                   (It could even be larger than the chosen FOV)
             		set stepSize      = 2   # Step size in mas (usually >0 to increase the image mask radius #oradiusNumber times)
					set oradiusNumber = 6   # Number of object mask radii to be testet (6 is mostly ok)

		- for "mat_cal_imarec_all.csh" only --------------------------------------------------------------------------
             		set muStarts      = (1.0 0.1 0.01) # Define the list of start regularization parameters muStart you want to use.
		- for "mat_cal_imarec_all.2.csh" only ------------------------------------------------------------------------
			set muStart0      = 1.0 # First value of muStart to be used.
						#  The later muStart values are calculated as: muStart(next) = muStart(actual) * muFactor0
                        set muFactor0     = 0.5 # Factor to calculate the next muStart
		
			    	set muFactor      = 0.5     # Multiplication factor to calculate the next regularization parameter mu in IRBis
			    	set muNumber      = 12      # Defines how often the multiplication factor is applied in IRBis (12 is normally ok)

			    	set qrecmode      = 2       # Selection criterion to be used for the estimation of the quality of the reconstructions
			            		            #  qrecmode = 1 : qrec is calculated using the $\chi^2$/residual ratios of both the measured phases (Closure or Fourier phase) and the squared visibilities.
					                   	    #  qrecmode = 2 : qrec is calculated using the $\chi^2$/residual ratios of the measured phases (Closure or Fourier phase) only
           			set regFuncs      = (-3 -4) # Define the numbers of the regularization functions you want to use (here you can test several reg.functions)
           			set weightPower   = 0.0     # Power for the uv density weight 
                                            #  (best experience with 0.0 for bispectrum and 0.5 for 0.0 complex visibilities)
           			set startmode     = 2       # Define then type of the start object: 0=read from file, 1=point source, 2=gaussian disk, 3=uniform disk, 
                                            #  4=fully darkened disk, 5=modified Lorentz function
           			set startparam    = 3.7     # Define the size of the start object chosen above.
                                   		    #  Use the value of the fit listed at the end of `Parameter.Estimation/data.parameter`
                                   		    #  (startmode= 2 -> FWHM [mas], 3 -> diameter [mas], 4 -> diameter [mas], 5 -> FWHM [mas])

	            --> some additional inputs:
				 	set model        = MODEL.fits  # For simulated interferometric data, you can insert the model fits file. [no]
                                               #  Define it's pixel scale in mas in "set modelPixelScale =".
             		set startima     = START.fits  # If you want to use another, better fitting, start image, for example, a model image generated by
		                				       #  a radiative transfer code, you can insert the fits file here, and "set startmode = 0" above. [no]
                                 		       #  Define it's pixel scale in mas in "set startimaPixelScale =".
                	set priorima     = PRIOR.fits  # If you want to use another prior image, you can insert a fits file here, and "set priormode = 0 above. [no]
                                               #  Define it's pixel scale in mas in "set priorimaPixelScale =".
                                               # If you don't want to use an inputfile, set it to "no"

		2) run the script:
	        `$SCRIPTS/mat_cal_imarec_all.csh  mat_cal_imarec_all.par`
	    or
		    `$SCRIPTS/mat_cal_imarec_all.2.csh  mat_cal_imarec_all.2.par`


		3) The outputs of the script `mat_cal_imarec_all.csh` and `mat_cal_imarec_all.2.csh` are described below


8. Short summary of the results of the scripts:

	1) The results are stored in result folders named 
        * `BIS.*.Script*.E.1/`   (`$algoMode=1` : only the bispectrum was used), 
        * `FT.*.Script*.E.1/`    (`$algoMode=2` : only the complex visibilities were used),
        * `FTBIS.*.Script*.E.1/` (`$algoMode=3` : both bispectrum and complex visibilities were used).
    If later the script is started again, the extention of the folder will automatically increase, i.e. -> .2 .3 .. (if the folder name
	could be identical to another one)
	    * `Script* == Script1` means produced by script `mat_cal_imarec_all.csh`,
	    * `Script* == Script2` means produced by script `mat_cal_imarec_all.2.csh`.

    2) In the result folder, each subfolder `E.1`, `E.2`, .. contains the best reconstruction of one run of IRBis
    Contents of the subfolders `E.1`, `E.2` ...:
	    - `rec_*.fits` : contains the "best" reconstruction of each run. This is the direct outcome of the recipe mat_cal_imarec (IRBis).
        It not only contains the unconvolved and convolved reconstruction, but also other data. See the manual of mat_cal_imarec for more.
		Man page call: `esorex --man-page mat_cal_imarec`
	   	- `bestrec.fits, bestrecconv.fits` : the "best" reconstruction without this other data, unconvolved and convolved, respectively.
	   	- `model.fits, modelconv.fits` : the model image (unconvolved/convolved) if specified.

	3) In the result folder, the textfile E.liste contains the $\chi^2$/ResidualRatio values and the image reconstruction parameters
    (regularisation function (Reg-Fct.), hyperparameter mu, etc..., for more details see Appendix below) for each IRBis run, i.e. in E.1, E.2, ...
	Each run is listed in one line; the runs are sorted with
	    1. increasing image quality parameter qrec (qrecmode=1: phase-visibiliy-qrec, derived from $\chi^2$ and residual ratios of the
	    squared visibilities and CP's or Fourier phases),
	    2. increasing image quality parameter qrec (qrecmode=2: phase-qrec, derived from $\chi^2$ and residual ratios of the phases (CPs or
		Fourier phases only)

	4) In the result folder the best reconstruction out of all reconstructions stored in E.1, E.2, etc.. is named:
	   	- `bestrec.qrec1.fits` (best according to the measure phase+visibility-qrec, i.e. qrecmode=1)
		- `bestrec.qrec2.fits` (best according to the measure phase-qrec, i.e. qrecmode=2)
	  (both are fits files and the direct outcome of the recipe mat_cal_imarec (i.e. the fits file rec_*.fits in the subdirectories E.*/), and
	   they contain, for example, the unconvolved, convolved reconstruction)

	5) In the result folder there are postscript files which contain all reconstructions (from E.1, E.2, ...) sorted according to the increasing
	   qrec measure as listed in E.liste:
		a) in folder qrecmode=1/ the reconstructions sorted according to the increasing measure phase+visibility-qrec
	 	b) in folder qrecmode=2/ the reconstructions sorted according to the increasing measure phase-qrec
		   the postscript files in these two folders are 
			- *.lin.ps (linear display)
			- *.sqrt.ps (sqrt display)
			- *.log.ps (log display)

	6) Furthermore, in the result folder there are postscript files for the case of image reconstruction simulations (if the known theoretical target
	is set in the *.par files), which display the correlation between the qrec parameters and the direct distance between the reconstructed and
	theoretical object (for the definition of the distance see dist1.0 in the Appendix below):
		- `*.distqrec.1.plot.ps` (for the image quality parameter qrec including phases and visibilties, i.e. qrecmode=1)
		- `*.distqrec.2.plot.ps` (for the image quality parameter qrec including phases only, i.e. qrecmode=2)

## Appendix: E.liste parameter  
	- qrec                   : image quality parameter derived from the Chi^2 and residual ratio values of V^2 and CP or absolut phase
	- cost                   : value of the cos function
	- chi2bis, rresbis       : Chi^2 and residual ratio of bispectrum
	- chi2vis2, rresvis2     : Chi^2 and residual ratio of squared visibility
	- chi2cp, rrescp         : Chi^2 and residual ratio of closure phase
	- dist1.0		 : distance between theoretical (m(x)) and reconstructed object (ok(x)); dist = Sqrt(Sum{ (m(x,y) - ok(x,y))^2 }/Sum{ m(x,y)^2 })
	- FOV (mas)		 : FOV of the reconstrcution
	- Reg-Fct.		 : number of the regularization function used
	- oradius step number	 : list of radii of the binary object masks; oradius(next) = oradius(actual) + step*i, with i = 0, number
	- mu factor number       : list of regularization parameter mu; mu(next) = mu(actual) * factor^j, with j = 0, number
	- calcVis2f0		 : =1: an artificial squared visibility and error for spatial frequency 0 should be calculated; =0: is not calculated
	- weightPower            : power for the uv density weight
	- npix			 : size of the reconstructed image in pixels
	- startmode	 	 : start image: 0 = read from file, 1 = point source, 2 = gaussian disc, 3 = uniform disc, 4 = fully darkened disc,
				   5 = modified Lorentz function
	- startparam             : startmode=0 -> scale [mas/px], mode=2 -> FWHM [mas], mode=3 -> diameter [mas], mode=4 -> diameter [mas], mode=5 -> FWHM [mas]
	- directory		 : reconstruction folder E.1, E.2, ...
	- cpqrec		 : image quality parameter derived from the Chi^2 and residual ratio values of the CP or absolut phase only

	
	man page of mat_cal_imarec: esorex --man-page mat_cal_imarec


j) Interferometric test data plus example parameter files:

   The directory TestData/ contains real and simulated interferometric data (CPs, V^2 in oifits files) of different targets plus the corresponding parameter files (mat_cal_imarec_all.par/mat_cal_imarec_all.2.par)
   to run the above described image reconstrcution scripts. The target directories contain the interferometric data in  data/ and some of the postscript files with sorted reconstrcutions in plots/.


