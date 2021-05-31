#!/bin/tcsh
#
#
# ----------- estimation of image reconstruction parameter --------------------------------------- A -------------------------
  set lambdaintervalls = $1
  set lambdaList0 = (); set lambdaList = ()
  set calcT3amp = 0
  if( $?CALCT3AMP ) then
    set calcT3amp = $CALCT3AMP
  endif
  set bltol = 0.05
  if( $?BLTOL ) then
    set bltol = $BLTOL
  endif
  set mjdtol = 0.0001
  if( $?MJDTOL ) then
    set mjdtol = $MJDTOL
  endif
  set fitfwhm = 2.35482004503095
  if( $?FITFWHM ) then
    set fitfwhm = $FITFWHM
  endif
  set i = 1
  while( $i <= $lambdaintervalls )
    set j = `echo $i | awk '{print 2*$1;}'`
    set jj = `echo $j | awk '{print $1+1;}'`
    set lambdaA = `echo $argv[1-]:q | awk -v j=$j '{print $j}'`
    set lambdaE = `echo $argv[1-]:q | awk -v jj=$jj '{print $jj}'`
    if($i <  $lambdaintervalls) set lambdaList0 = ($lambdaList0$lambdaA,$lambdaE,)
    if($i == $lambdaintervalls) set lambdaList0 = ($lambdaList0$lambdaA,$lambdaE)
    set lambdaList = ($lambdaList $lambdaA $lambdaE)
  @ i++
  end
  set nr = `echo $lambdaintervalls | awk '{print 2*$1+2;}'`
  set data = ($argv[$nr-])

  set info          = "istd,param,result,prepare"
  ls -l $data
  set est = Parameter.Estimation
  if( -d $est ) rm -rf $est
  mkdir $est
  set estft = Parameter.Estimation.FT
  if( -d $estft ) rm -rf $estft
  mkdir $estft
  set DARGSbis = (--nbresult=3 --lambda_list=$lambdaList0 \
	     --algo_mode=1 --calc_t3amp=$calcT3amp --calc_vis2f0=0 --model_scale=0.1 \
	     --om_start=10 --om_step=1 --om_count=1 \
	     --mu_start=1.0 --mu_factor=1.0 --mu_count=1 --mjd_tol=$mjdtol --bl_tol=$bltol --fit_fwhm=$fitfwhm \
	     --cost_func=3 --cost_weight=1.0 --reg_func=0 --reg_eps=1.0 \
	     --grad_tol=0.00000001 --weight_power=0.0 --conv_scale=1.0 --precision=0 --info_flags=$info)
  set DARGSft = (--nbresult=3 --lambda_list=$lambdaList0 \
             --algo_mode=2 --calc_t3amp=$calcT3amp --calc_vis2f0=0 --model_scale=0.1 \
             --om_start=10 --om_step=1 --om_count=1 \
             --mu_start=1.0 --mu_factor=1.0 --mu_count=1 --mjd_tol=$mjdtol --bl_tol=$bltol --fit_fwhm=$fitfwhm \
             --cost_func=3 --cost_weight=1.0 --reg_func=0 --reg_eps=1.0 \
             --grad_tol=0.00000001 --weight_power=0.0 --conv_scale=1.0 --precision=0 --info_flags=$info)
  set NSTART = (--start_mode=2 --start_param=2.0)
  set NPRIOR = (--prior_mode=2 --prior_param=2.0)
  set NFOV   = (--fov=50) 
  set NDIM   = (--npix=32)
  set OUTPUT = (--vis2_name=${est}/A.vis2.dat --cp_name=${est}/A.cp.dat --vis_name=${estft}/A.ft.dat)
  set ARGS = ($NSTART $NPRIOR $NFOV $NDIM $OUTPUT)

  ls -l $data
  set anzahl = $#data
  set oits = ()
  rm -f sof
  set i = 1
  while($i <= $anzahl)
    echo "$data[$i] TARGET_CAL_INT" >> sof
#   estimation of the data types in the oifits file and therefore which image reconstructions are possilble
    rm -f txt0; $SCRIPTS/DataType.csh $data[$i] >> txt0
    set oit = `awk '{if ($0 ~ /# oit = /) {print $4;}}' txt0`
    set oits = ($oits $oit)
  @ i++
  end
  set posbiss = ()
  set posfts  = ()
  foreach oit ($oits)
    set posbis = 0
    set posft  = 0
    echo "oit = $oit"
    if(($oit == 1)||($oit == 2)) set posbis = 1
    if(($oit == 1)||($oit == 3)||($oit == 4)) set posft  = 1
    set posbiss = ($posbiss $posbis)
    set posfts  = ($posfts  $posft)
  end
# estimation, if OI_T3 and OI_VIS2 are available in each OIFITS file  --> bispectrum reconstruction is tested below
  set posbis = 1
  foreach posbis0 ($posbiss)
    set posbis1 = `echo $posbis $posbis0 | awk '{print $1*$2;}'`; set posbis = $posbis1
  end
# estimation, if OI_VIS and OI_VIS2 are available in each OIFITS file  --> complex visibility reconstruction is tested below
  set posft  = 1
  foreach posft0 ($posfts)
    set posft1 = `echo $posft $posft0 | awk '{print $1*$2;}'`; set posft = $posft1
  end
#  echo "posbis, posft = " $posbis $posft
#  echo "weiter [cr]"; set ja = $<

  if($posbis == 1) esorex --log-dir=${est} --output-dir=${est} mat_cal_imarec $DARGSbis $ARGS sof
  if($posft  == 1) esorex --log-dir=${estft} --output-dir=${estft} mat_cal_imarec $DARGSft $ARGS sof
  if( $status != 0 ) exit 1
  if($posbis == 0) then
    echo "no OI_VIS2 and OI_T3 data available ---> EXIT"
    exit 1
  endif

  if($posft == 1) then
    cd $estft
    cp A.ft.dat ../$est
    cd ..
  endif 
  rm -rf $estft
  cd $est
  sed -i 's/\[tid=000\]//' esorex.log
#    estimation of the optimal FOV in pixel and mas
                    set tgauss = `awk '{if ($0 ~ /mat_fit_start_image_vis2: fit gaussian /) {print $12;}}' esorex.log`
  if($#tgauss == 0) set tgauss = `awk '{if ($0 ~ /mat_fit_start_image_vis2:  fit gaussian /) {print $12;}}' esorex.log`
                        set tchi2gauss = `awk '{if ($0 ~ /mat_fit_start_image_vis2: fit gaussian /) {print $19;}}' esorex.log`
  if($#tchi2gauss == 0) set tchi2gauss = `awk '{if ($0 ~ /mat_fit_start_image_vis2:  fit gaussian /) {print $19;}}' esorex.log`
                 set tud    = `awk '{if ($0 ~ /mat_fit_start_image_vis2: fit uniform disc /) {print $12;}}' esorex.log`
  if($#tud == 0) set tud    = `awk '{if ($0 ~ /mat_fit_start_image_vis2:  fit uniform disc /) {print $12;}}' esorex.log`
                     set tchi2ud    = `awk '{if ($0 ~ /mat_fit_start_image_vis2: fit uniform disc /) {print $19;}}' esorex.log`
  if($#tchi2ud == 0) set tchi2ud    = `awk '{if ($0 ~ /mat_fit_start_image_vis2:  fit uniform disc /) {print $19;}}' esorex.log`
                  set tfdd   = `awk '{if ($0 ~ /mat_fit_start_image_vis2: fit fully darkened /) {print $13;}}' esorex.log`
  if($#tfdd == 0) set tfdd   = `awk '{if ($0 ~ /mat_fit_start_image_vis2:  fit fully darkened /) {print $13;}}' esorex.log`
                      set tchi2fdd   = `awk '{if ($0 ~ /mat_fit_start_image_vis2: fit fully darkened /) {print $20;}}' esorex.log`
  if($#tchi2fdd == 0) set tchi2fdd   = `awk '{if ($0 ~ /mat_fit_start_image_vis2:  fit fully darkened /) {print $20;}}' esorex.log`
                 set tld    = `awk '{if ($0 ~ /mat_fit_start_image_vis2: fit Lorentz /) {print $12;}}' esorex.log`
  if($#tld == 0) set tld    = `awk '{if ($0 ~ /mat_fit_start_image_vis2:  fit Lorentz /) {print $12;}}' esorex.log`
                     set tchi2ld    = `awk '{if ($0 ~ /mat_fit_start_image_vis2: fit Lorentz /) {print $19;}}' esorex.log`
  if($#tchi2ld == 0) set tchi2ld    = `awk '{if ($0 ~ /mat_fit_start_image_vis2:  fit Lorentz /) {print $19;}}' esorex.log`
                  set tfov   = `awk '{if ($0 ~ /mat_guess_parameters: fov   = /) {print $0;}}' esorex.log`
  if($#tfov == 0) set tfov   = `awk '{if ($0 ~ /mat_guess_parameters:  fov   = /) {print $0;}}' esorex.log`
                   set tnpix  = `awk '{if ($0 ~ /mat_guess_parameters: npix  = /) {print $0;}}' esorex.log | tr '*' x`
  if($#tnpix == 0) set tnpix  = `awk '{if ($0 ~ /mat_guess_parameters:  npix  = /) {print $0;}}' esorex.log | tr '*' x`
                   set ttfov  = `awk '{if ($0 ~ /mat_guess_parameters: --fov=/) {print $0;}}' esorex.log`
  if($#ttfov == 0) set ttfov  = `awk '{if ($0 ~ /mat_guess_parameters:  --fov=/) {print $0;}}' esorex.log`
                    set ttnpix = `awk '{if ($0 ~ /mat_guess_parameters: --npix=/) {print $0;}}' esorex.log`
  if($#ttnpix == 0) set ttnpix = `awk '{if ($0 ~ /mat_guess_parameters:  --npix=/) {print $0;}}' esorex.log`

# ----------- Plott der gemessenen Visibility und der gefittetn Modelle ------------------------------------- A -----------------
set vis20 = `echo A.vis2.dat`
set vis2 = $vis20.new; rm -f $vis2; awk '{ if($1!="#") {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$12,$13,$14,$11;}; }' $vis20 > $vis2
ls -l $vis2
rm -f $vis2.0
awk 'BEGIN{z=0;} { if(($1=="#")&&(z==0)) {print "# Nr. |", "U |", "V |", "Vis |", "VisErr |", "Gaussian |", "Uniform d |", "Fully dark d |", "Vis (1.rec) |", "Vis (2.rec) |", "Vis (3.rec) |", "f (1/arcsec) |", "Lorentz ";}; \
                               if(($1=="#")&&(z==1)) {print "# 1   |", "2 |", "3 |", "4   |", "5      |", "6        |", "7        |", "8            |", "9           |", "10          |", "11          |", "12          |", "13   ";}; \
                               if($1!="#") {rad=sqrt($3^2+$4^2); radc=rad*1000.; vis=sqrt(sqrt($6^2)); if($6==0.0) {viserr=0.0;}; if($6!=0.0) {viserr=0.5*$7/vis;}; \
                                            visgauss=sqrt(sqrt($8^2)); visud=sqrt(sqrt($9^2)); visfdd=sqrt(sqrt($10^2)); visld=sqrt(sqrt($14^2)); visa=sqrt(sqrt($11^2)); visb=sqrt(sqrt($12^2)); visc=sqrt(sqrt($13^2)); \
                                            print $1,$3,$4,vis,viserr,visgauss,visud,visfdd, visa,visb,visc,radc,visld;}; z=z+1;}' $vis2 > $vis2.0

rm -f $vis2.0.sorted ttt; awk '{if($1!="#") {print $0;};}' $vis2.0 > ttt; sort -n -k 12 ttt > $vis2.0.sorted

$SCRIPTS/gnupl2.csh $vis2.0.sorted "measured visibilities" 12 4 5 "spatial frequency (1/arcsec)" "Visibility" gaussudfdda.ps \
                            "fitted Uniform disk model (diameter={$tud}mas)" 7 "fitted Gaussian model (FWHM={$tgauss}mas)" 6 "fitted Fully darkened disk model (diameter={$tfdd}mas)" 8 \
                            "fitted Lorentz function (FWHM={$tld}mas)" 13

# ----------- Plott der gemessenen Visibility und der gefittetn Modelle ------------------------------------- E -----------------

  set vis2 = A.vis2.dat
  set cp   = A.cp.dat
  set ft   = A.ft.dat
#  set fits = `echo *.fits` 
  set fits = ()
  foreach ff ($data)
     set fits = ($fits $ff:t)
  end
# cp:
# measured closure phases vs reconstructed closure phases
# nr lambda u1 v1 u2 v2 amp cp cperr biserr rec1_amp rec1_cp rec2_amp rec2_cp rec3_amp rec3_cp
# 1  2      3  4  5  6  7   8  9     10     11       12      13       14      15       16
# (lambda in m)
# cp: neu seit 13-11-15:
# nr lambda u1 v1 u2 v2 amp amperr cp cperr biserr rec1_amp rec1_cp rec2_amp rec2_cp rec3_amp rec3_cp
# 1  2      3  4  5  6  7   8      9  10    11     12       13       14      15      16       17
#
# vis2:
# measured squared visibilities vs reconstructed squared visibilities
# nr lambda u v bl vis2 err gd ud fdd rec1 rec2 rec3
# 1  2      3 4 5  6    7   8  9  10  11   12   13
# ft:
# measured complex visibilities vs reconstructed complex visibilities
# nr lambda u v bl amp amperr phi phierr gd ud fdd ld rec1_amp rec1_phi rec2_amp rec2_phi rec3_amp rec3_phi
# 1  2      3 4 5  6   7      8   9      10 11 12  13 14       15       16       17       18       19
#
# (lambda in m)
     rm -f lambdas; awk '{ if($1!="#") {print $2*1000000.,1.;}; }' $vis2 > lambdas
     rm -f wavelengths.ps; $SCRIPTS/gplot -psfl -pscolor -o wavelengths.ps -grid -ds impulses -headline "Plot of all wavelengths in the data" -xlabel "wavelengths (in mu)" -x 1 -y 2 lambdas
     set minbaseline = `awk 'BEGIN{min=1e36;} { if(($1!="#")&&($5>0.)) {frad=sqrt($3*$3+$4*$4); fac=$5/frad; ux=$3*fac; uy=$4*fac; rad=sqrt(ux^2+uy^2); if(rad<min) {min=rad;};}; } END{print min;}' $vis2`
     set maxbaseline = `awk 'BEGIN{max=-10.;} { if(($1!="#")&&($5>0.)) {frad=sqrt($3*$3+$4*$4); fac=$5/frad; ux=$3*fac; uy=$4*fac; rad=sqrt(ux^2+uy^2); if(rad>max) {max=rad;};}; } END{print max;}' $vis2`
     set minlambda   = `awk 'BEGIN{min=1e36;} { if(($1!="#")&&($5>0.)) {lam=$2; if(lam<min) {min=lam;};}; } END{print min;}' $vis2`
     set maxlambda   = `awk 'BEGIN{max=-10.;} { if(($1!="#")&&($5>0.)) {lam=$2; if(lam>max) {max=lam;};}; } END{print max;}' $vis2`
     set vis2anzahl  = `awk 'BEGIN{z=0;}      { if(($1!="#")&&($5>0.)) {z=z+1;}; } END{print z;}' $vis2`
     set cpanzahl    = `awk 'BEGIN{z=0;}      { if($1!="#") {z=z+1;}; } END{print z;}' $cp`
     set ftanzahl = 0; if($posft == 1) set ftanzahl    = `awk 'BEGIN{z=0;}      { if(($1!="#")&&($5>0.)) {z=z+1;}; } END{print z;}' $ft`
     set vis2snr     = `awk 'BEGIN{sum1=0.; sum2=0.;} { if(($1!="#")&&($5>0.)) {vis2=sqrt($6^2); err=sqrt($7^2); if(err>0.) {snr=vis2/err; sum1=sum1+snr; sum2=sum2+1.;}; }; } END{snr=sum1/sum2; print snr;}' $vis2`
     set cperr       = `awk 'BEGIN{sum1=0.; sum2=0.;} { if($1!="#") {pi=atan2(1.,1.)*4.; err=sqrt($10^2)*180./pi; sum1=sum1+err; sum2=sum2+1.;}; } END{err=sum1/sum2; print err;}' $cp`
     set ftpherr = 0; if($posft == 1) set ftpherr     = `awk 'BEGIN{sum1=0.; sum2=0.;} { if($1!="#") {pi=atan2(1.,1.)*4.; err=sqrt($9^2)*180./pi; sum1=sum1+err; sum2=sum2+1.;}; } END{err=sum1/sum2; print err;}' $ft`

     set pi = `echo 1. | awk '{ pi=atan2($1,$1)*4.; print pi; }'`
     set alpha = 0.25
     set maxobjsize = `echo $minbaseline $maxlambda $pi | awk '{ size = (1.22/2.)*($2/$1)*180.*60.*60.*1000./$3; print size; }'`
#    set maxfov = `echo $minbaseline $maxlambda $pi | awk '{ fov = ($2/$1)*180.*60.*60.*1000./$3; print fov; }'`
     set maxfov = `echo $minbaseline $maxlambda $pi | awk '{ fov = 2.44*($2/$1)*180.*60.*60.*1000./$3; print fov; }'`
#    set npix0  = `echo $minbaseline $maxbaseline $minlambda $maxlambda $alpha | awk '{ npix = (1./$5)*($2/$1)*($4/$3); print npix; }'`
#    set npix0  = `echo $minbaseline $maxbaseline $minlambda $maxlambda $alpha | awk '{ npix = (2.44/$5)*($2/$1)*($4/$3); print npix; }'`
     set npix0  = `echo $minbaseline $maxbaseline $minlambda $maxlambda $alpha | awk '{ npix = (2./$5)*($2/$1)*($4/$3); print npix; }'`
     set maxfov2 = `echo $maxfov | awk '{print $1/2;}'`

#    max. FOV, wenn max. spatiale Frequenz bei npix/4 liegen soll: 
#      --> max. FOV[mas] = 0.25 * npix * minlambda[m]/maxbaseline[m] * (180*60*60*1000)/pi
#    set maxfov16 =  `echo $maxbaseline $minlambda $pi 16 $alpha | awk '{ fov=$5*$4*($2/$1)*180.*60.*60.*1000./$3; print fov; }'`
     set maxfov16 =  `echo $maxbaseline $minlambda $pi 16 $alpha | awk '{ fov=$5*$4*($2/$1)*180.*60.*60.*1000./$3; print fov; }'`  
     set maxfov32 =  `echo $maxbaseline $minlambda $pi 32 $alpha | awk '{ fov=$5*$4*($2/$1)*180.*60.*60.*1000./$3; print fov; }'`  
     set maxfov64 =  `echo $maxbaseline $minlambda $pi 64 $alpha | awk '{ fov=$5*$4*($2/$1)*180.*60.*60.*1000./$3; print fov; }'`  
     set maxfov128 =  `echo $maxbaseline $minlambda $pi 128 $alpha | awk '{ fov=$5*$4*($2/$1)*180.*60.*60.*1000./$3; print fov; }'`  
     set maxfov256 =  `echo $maxbaseline $minlambda $pi 256 $alpha | awk '{ fov=$5*$4*($2/$1)*180.*60.*60.*1000./$3; print fov; }'`  
     set maxfov512 =  `echo $maxbaseline $minlambda $pi 512 $alpha | awk '{ fov=$5*$4*($2/$1)*180.*60.*60.*1000./$3; print fov; }'`  

#     set maxfov16 =  `echo $maxfov $npix0  16 | awk '{ fov=($3/$2)*$1; print fov; }'`
#     set maxfov32 =  `echo $maxfov $npix0  32 | awk '{ fov=($3/$2)*$1; print fov; }'`
#     set maxfov64 =  `echo $maxfov $npix0  64 | awk '{ fov=($3/$2)*$1; print fov; }'`
#     set maxfov128 = `echo $maxfov $npix0 128 | awk '{ fov=($3/$2)*$1; print fov; }'`
#     set maxfov256 = `echo $maxfov $npix0 256 | awk '{ fov=($3/$2)*$1; print fov; }'`
#     set maxfov512 = `echo $maxfov $npix0 512 | awk '{ fov=($3/$2)*$1; print fov; }'`

     set fovminbl  = $tfov[8]     # FOV[mas] calculated from the smallest baseline and max. wavelength
     set npixminbl = $tnpix[8]    # calculated from $fovminbl 

	set res = `echo $maxbaseline $minlambda $maxlambda $pi | awk '{ lam=1e6*($2+$3)/2.; res = (lam/$1) * (18.*6.*6.)/$4; print res; }'`

     rm -f data.parameter
     echo "# parameter of the interferometric data in the oifits file(s):" >> data.parameter
     echo "# $fits " >> data.parameter
     echo "# ====================================================================================================" >> data.parameter
     echo " " >> data.parameter
     echo "A) Information about the data: " >> data.parameter
     echo "   - minimum and maximum wavelength (in m)     : $minlambda $maxlambda " >> data.parameter
     echo "   - minimum and maximum baseline length (in m): $minbaseline $maxbaseline " >> data.parameter
     echo "   - resolution lambda / B_max (in mas)        :  $res " >> data.parameter
     echo "   - number of measured squared visibilities   : $vis2anzahl " >> data.parameter
     echo "   - number of measured closure phases         : $cpanzahl " >> data.parameter
     echo "   - number of measured Fourier phases         : $ftanzahl " >> data.parameter
     echo "   - average SNR of measured squared visibilities: $vis2snr " >> data.parameter
     echo "   - average error of the closure phase (in deg) : $cperr " >> data.parameter
     echo "   - average error of the Fourier phase (in deg) : $ftpherr " >> data.parameter
     if( ($posbis == 1)&&($posft == 1) ) echo " --> closure phase imaging and absolut phase imaging are possible!" >> data.parameter
     if( ($posbis == 0)&&($posft == 1) ) echo " --> absolut phase imaging is possible only!" >> data.parameter
     if( ($posbis == 1)&&($posft == 0) ) echo " --> closure phase imaging is possible only!" >> data.parameter
     if( ($posbis == 0)&&($posft == 0) ) echo " --> no imaging is possible!" >> data.parameter
     echo " " >> data.parameter
     echo "B) Estimated target size by fitting the V^2 data: " >> data.parameter
     echo "   * (2) Gaussian            --> FWHM = $tgauss mas (red. Chi^2 = $tchi2gauss) " >> data.parameter
     echo "   * (3) Uniform disk        --> diameter = $tud mas (red. Chi^2 = $tchi2ud) " >> data.parameter
     echo "   * (4) Fully darkened disk --> diameter = $tfdd mas  (red. Chi^2 = $tchi2fdd) " >> data.parameter
     echo "   * (5) Lorentzian function --> FWHM = $tld mas (red. Chi^2 = $tchi2ld) " >> data.parameter
     echo " " >> data.parameter
     echo "C) Recommended size of the angular FOV and the size of the NxN pixel grid for the image reconstruction run: " >> data.parameter
     echo "    - the optimal FOV[mas] should be about ~ 2 to 4 times the size of the target (see estimations in B)) " >> data.parameter
     echo "    - the optimal FOV[mas] is covered by a NxN pixel grid where in the Fourier plane all uv points " >> data.parameter
     echo "      are lying within the cut-off frequency f_pixel=N/4, i.e. where all uv points have distances to " >> data.parameter
     echo "      the Fourier center which are smaller or equal to N/4 (to avoid aliasing) " >> data.parameter
     echo "    - the maximum target size derived from the present uv coverage for reliable image reconstruction " >> data.parameter
     echo "      should be not larger than Theta_max = (1.22/2)*lambda_max/B_min = $maxobjsize mas " >> data.parameter
     echo "      (in order to have a few data points within the first loop of the Fourier transform of the target, e.g. a Uniform Disk.) " >> data.parameter
     echo "      --> therefore the max. FOV[mas] for image reconstruction should be ~ 2 to 4 x Theta_max, i.e. between $maxfov2 and $maxfov mas " >> data.parameter
     echo " " >> data.parameter
     echo "==> a collection of NxN pixel grids and their FOVs[mas] providing that the uv point with the largest distance to the Fourier center" >> data.parameter
     echo "    has a distance of N/4 pixels to it: " >> data.parameter
     echo "      * 16x16 pixels --> FOV = $maxfov16 mas " >> data.parameter
     echo "      * 32x32 pixels --> FOV = $maxfov32 mas " >> data.parameter
     echo "      * 64x64 pixels --> FOV = $maxfov64 mas " >> data.parameter
     echo "      * 128x128 pixels --> FOV = $maxfov128 mas " >> data.parameter
     echo "      * 256x256 pixels --> FOV = $maxfov256 mas " >> data.parameter
     echo "      * 512x512 pixels --> FOV = $maxfov512 mas " >> data.parameter
     echo " Note: a) Choose that NxN pixel grid, where the corresponding FOV is nearly equal to about 2 to 4 times the target size " >> data.parameter
     echo "       b) For the next runs you may want to use larger NxN pixel grids (this may improve the fits to the data): if you take " >> data.parameter
     echo "           a smaller FOV for the chosen NxN pixel grid as stated in the list above, this means that the length of uv point " >> data.parameter
     echo "           of the longest baseline is smaller than N/4, which is OK" >> data.parameter
     echo " " >> data.parameter
     echo "==>    settings in the image reconstruction script: " >> data.parameter
     echo "        - set the variable 'npix' to the size of the chosen NXN pixel grid, " >> data.parameter
     echo "        - set the variable 'fov' to ~ 2 to 4 times the target size (= FOV[mas] of the image reconstruction run) " >> data.parameter
     echo " " >> data.parameter
     echo " " >> data.parameter
     
     $SCRIPTS/uvplot.csh

     # Erzeugung des ASCII-Files *.vis2 (altes Format fuer Fortran-Code von IRBis):
     # Spaltenbelegung von *.vis2:
     #  ID     Vis^2   Vis^2Err   u[m]  v[m]
     set lam0 = `awk 'BEGIN{sum1=0.; sum2=0.;} { if($1!="#") {sum1=sum1+$2; sum2=sum2+1.;}; } END{print (sum1/sum2)*1000000.;}' $vis2`
     rm -f $vis2:r.vis2.0; awk '{ if($1!="#") {frad=sqrt($3*$3+$4*$4); fac=$5/frad; ux=$3*fac; uy=$4*fac; print $1, $6, $7, ux,uy;}; }' $vis2 > $vis2:r.vis2.0
     rm -f head; echo "# uv-Vektoren zur Bezugswellenlaenge $lam0 mu" >> head
     rm -f $vis2:r.vis2; cat head $vis2:r.vis2.0 > $vis2:r.vis2

     if( $?JPLOT ) then
      if( $JPLOT == 1 ) then
       # Darstellung der Daten
##       $SCRIPTS/plotvis2.csh     # --> $vis2.ps  =  A.vis2.dat.ps
##       $SCRIPTS/plotvis.csh      # -->              A.vis.dat.ps
##       $SCRIPTS/plotft.csh       # --> $ft.ps    =  A.ft.dat.ps
##       $SCRIPTS/plotcp.csh       # --> $cp.ps    =  A.cp.dat.ps
      endif
     endif

  cd ..

  echo "- some proposed image reconstruction parameter are listed in data.parameter "
  echo "- uv coverage as postscript plot uv.ps "
  ls -l $est/data.parameter $est/uv.ps $est/wavelengths.ps $est/gaussudfdda.ps
  if( $?JPLOT ) then
##    if( $JPLOT == 1 ) ls -l $est/A.vis.dat.ps $est/A.ft.dat.ps $est/A.cp.dat.ps
  endif

# ----------- estimation of image reconstruction parameter --------------------------------------- E -------------------------
