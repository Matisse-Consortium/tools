#!/bin/tcsh
#
#
# ------------------------------------------------------------------------------------
# ----------- START of INPUT ------------------------------------------------------------
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# ----------- END of INPUT ------------------------------------------------------------
# ------------------------------------------------------------------------------------

set parfile = $1

source $parfile

if( $model != "no" ) then
  rm -f output2 ; fdump $model \!output2 NAXIS1 - ; set cdelt1 = `awk '{if ($0 ~ /CDELT1 /) {print $3;};}' output2`
  rm -f output2 ; fdump $model \!output2 NAXIS1 - ; set cdelt2 = `awk '{if ($0 ~ /CDELT2 /) {print $3;};}' output2`
  rm -f output2 ; fdump $model \!output2 NAXIS1 - ; set naxis1 = `awk '{if ($0 ~ /NAXIS1 /) {print $3;};}' output2`
  rm -f output2 ; fdump $model \!output2 NAXIS1 - ; set naxis2 = `awk '{if ($0 ~ /NAXIS2 /) {print $3;};}' output2`
  set qq1 = `echo $cdelt1 | awk '{ qq = 1; if($1 < 0.) {qq = -1;}; print qq; }'`
  set qq2 = `echo $cdelt2 | awk '{ qq = 1; if($1 < 0.) {qq = -1;}; print qq; }'`
  if( ($qq1 < 0) && ($qq2 < 0) ) then
    rm -f $model.99; ftcopy ''$model'['$naxis1':1,'$naxis2':1]' $model.99
  endif
  if( ($qq1 < 0) && ($qq2 > 0) ) then
    rm -f $model.99; ftcopy ''$model'['$naxis1':1,1:'$naxis2']' $model.99
  endif
  if( ($qq1 > 0) && ($qq2 < 0) ) then
    rm -f $model.99; ftcopy ''$model'[1:'$naxis1':1,'$naxis2':1]' $model.99
  endif
  if( ($qq1 > 0) && ($qq2 > 0) ) then
    rm -f $model.99; cp $model $model.99
  endif
  set model0 = $model.99
  ls -l $model0
else
  set model0 = $model
endif

# set data = ()
# foreach f ($data0)
#    set data = ($data $f)
# end
ls -l $data
if( $status != 0 ) exit 1

rm -rf E.*

unset lambdaList0
set lambdaList0 = ()
set lambdaintervalls = `echo $#lambdaList | awk '{print int($1/2);}'`
set lambdas = `echo $lambdaintervalls | awk '{print 2*$1;}'`
if( $#lambdaList != $lambdas ) then
  echo "Anzahl der Intervallgrenzen ist ungerade: $#lambdaList"
  goto Label999
endif
set lambdaFrom = $lambdaList[1]
set lambdaTo   = $lambdaList[$#lambdaList]
set i = 1
while( $i <= $lambdaintervalls )
  set j = `echo $i | awk '{print 2*$1-1;}'`
  set jj = `echo $j | awk '{print $1+1;}'`
  if($i <  $lambdaintervalls) set lambdaList0 = ($lambdaList0$lambdaList[$j],$lambdaList[$jj],)
  if($i == $lambdaintervalls) set lambdaList0 = ($lambdaList0$lambdaList[$j],$lambdaList[$jj])
@ i++
end
echo $lambdaList0

# ----------- estimation of image reconstruction parameter --------------------------------------- A -------------------------
if($guess != 0 ) then
  set jplot = 0
  setenv JPLOT $jplot
  setenv CALCT3AMP $calcT3amp2
  setenv MJDTOL $mjdtol2
  setenv BLTOL  $bltol2
  setenv FITFWHM $fitfwhm2
  $SCRIPTS/parameter.est.csh $lambdaintervalls $lambdaList $data
  if( $status != 0 ) exit 1
#  cp $0 Parameter.Estimation/
  cp $parfile Parameter.Estimation/
  goto Label999
endif
# ----------- estimation of image reconstruction parameter --------------------------------------- E -------------------------

echo "------------1-------------"
if($algoMode == 1) then
  set resultdir = BIS.$objname.F$fov.pix$npix.Cost$costFunc.weightPow$weightPower.qrec$qrecmode.mu$muStart0.reg$regFuncs[1].wiener$wienerfilter.engine$engine.Script2.E.1
  set name = BIS.$objname.F$fov.pix$npix.weightPow$weightPower.qrec$qrecmode.mu$muStart0.reg$regFuncs[1].wiener$wienerfilter.engine$engine
endif
if($algoMode == 2) then
  set resultdir = FT.$objname.F$fov.pix$npix.Cost$costFunc.weightPow$weightPower.qrec$qrecmode.mu$muStart0.reg$regFuncs[1].wiener$wienerfilter.engine$engine.Script2.E.1
  set name = FT.$objname.F$fov.pix$npix.weightPow$weightPower.qrec$qrecmode.mu$muStart0.reg$regFuncs[1].wiener$wienerfilter.engine$engine
endif
if($algoMode == 3) then
  set resultdir = FTBIS.$objname.F$fov.pix$npix.Cost$costFunc.weightPow$weightPower.qrec$qrecmode.mu$muStart0.reg$regFuncs[1].wiener$wienerfilter.engine$engine.Script2.E.1
  set name = FTBIS.$objname.F$fov.pix$npix.weightPow$weightPower.qrec$qrecmode.mu$muStart0.reg$regFuncs[1].wiener$wienerfilter.engine$engine
endif
set i = 1
label1:
if(-d $resultdir) then
  echo "$resultdir existiert bereits."
  @ i++
  set resultdir = $resultdir:r.$i
  goto label1
endif
echo "$resultdir existiert noch nicht."
mkdir -p $resultdir
echo $resultdir

echo "------------2-------------"

set starttime = `date`
Label111:

foreach regFunc ( $regFuncs )
# ------------ loop over chosen regularization functions ------------------------------- A --------------
set qrecvorher = 1e36
set muStart = $muStart0
while( 1 > 0 )
# foreach muStart ( $muStarts )
# ------------ loop over chosen start values of mu ----------------------------------- A ------------


ls -l $data $model0
set anzahl = $#data
rm -f sof
set i = 1
while($i <= $anzahl)
  echo "$data[$i] TARGET_CAL_INT" >> sof
@ i++
end 
if( $model0 != "no" ) then
  echo "$model0 MODEL_IMAGE" >> sof
endif
if( $startima != "no" ) then
  echo "$startima START_IMAGE" >> sof
  set startmode = 0
  set startparam = $startimaPixelScale 
endif
if( $priorima != "no" ) then
  echo "$priorima PRIOR_IMAGE" >> sof
  set priormode = 0
  set priorparam = $priorimaPixelScale
endif

echo "---1---"

#  --nbresult            : Number of reconstructions written to the result file.
#                          The best result is always stored. If 0 is given, all
#                          created reconstructions are stored in the result
#                          file. [0]
#  --vis2_name           : ASCII file for measured and reconstructed squared
#                          visibilities. []
#  --cp_name             : ASCII file for measured and reconstructed closure
#                          phases. []
#  --vis_name            : ASCII file for measured and reconstructed complex visibilities. []
# ------------------------------------------------------------------------------------
if( $engine > 1 ) then
  set DARGS = (--nbresult=3 --lambda_list=$lambdaList0 \
	       --algo_mode=$algoMode \
	       --calc_visamp=$calcVisamp --calc_visf0=$calcVisf0 \
	       --calc_t3amp=$calcT3amp --calc_vis2f0=$calcVis2f0 --model_scale=$modelPixelScale \
	       --om_start=$oradiusStart --om_step=$stepSize --om_count=$oradiusNumber --engine=$engine \
	       --mu_start=$muStart --mu_factor=$muFactor --mu_count=$muNumber --wiener_filter=$wienerfilter --fit_fwhm=$fitfwhm \
	       --cost_func=$costFunc --cost_weight=$costWeight --reg_func=$regFunc --reg_eps=$regEps --start_select=$startselect --prior_select=$priorselect \
	       --grad_tol=$gradTol --pg_tol=$pgTol --factr=$factr --ncorr=$ncorr --weight_power=$weightPower --conv_scale=$convScale --precision=$precision --info_flags=$info \
	       --mjd_tol=$mjdtol --bl_tol=$bltol --filter_fwhm=$filterfwhm --filter_factor=$filterfactor --guess=$guess --noise_factor=$noisefactor)
else
  set DARGS = (--nbresult=3 --lambda_list=$lambdaList0 \
               --algo_mode=$algoMode \
               --calc_visamp=$calcVisamp --calc_visf0=$calcVisf0 \
               --calc_t3amp=$calcT3amp --calc_vis2f0=$calcVis2f0 --model_scale=$modelPixelScale \
               --om_start=$oradiusStart --om_step=$stepSize --om_count=$oradiusNumber \
               --mu_start=$muStart --mu_factor=$muFactor --mu_count=$muNumber --wiener_filter=$wienerfilter --fit_fwhm=$fitfwhm \
               --cost_func=$costFunc --cost_weight=$costWeight --reg_func=$regFunc --reg_eps=$regEps --start_select=$startselect --prior_select=$priorselect \
               --grad_tol=$gradTol --pg_tol=$pgTol --factr=$factr --ncorr=$ncorr --weight_power=$weightPower --conv_scale=$convScale --precision=$precision --info_flags=$info \
               --mjd_tol=$mjdtol --bl_tol=$bltol --filter_fwhm=$filterfwhm --filter_factor=$filterfactor --guess=$guess --noise_factor=$noisefactor)
endif
# ------------------------------------------------------------------------------------
echo $DARGS

echo "---2---"
# ------------------------------------------------------------------------------------
set NSTART = (--start_mode=$startmode --start_param=$startparam)
set NPRIOR = (--prior_mode=$priormode --prior_param=$priorparam)
set NFOV   = (--fov=$fov)
set NDIM   = (--npix=$npix)
set results = E.1
set i = 1
label2:
if(-d $results) then
  echo "$results existiert bereits."
  @ i++
  set results = $results:r.$i
  goto label2
endif
echo "$results existiert noch nicht."
mkdir -p $results
echo $results

set OUTPUT = (--vis2_name=${results}/A.vis2.dat --cp_name=${results}/A.cp.dat --vis_name=${results}/A.ft.dat)
# ------------------------------------------------------------------------------------

echo "---3---"

set ARGS = ($NSTART $NPRIOR $NFOV $NDIM $OUTPUT)
echo $ARGS

echo "---4---"

esorex --log-dir=${results} --output-dir=${results} mat_cal_imarec $DARGS $ARGS sof
if( $status != 0 ) exit 1

# cp $0 ${results}/
cp $parfile ${results}/

cd ${results}/
#  if($algoMode == 1) then
  if(($algoMode == 1)||($costFunc == 3)) then
   echo "-------------- mat_cal_imarec_all.2.csh : algoMode == 1 || costFunc == 3 -------------------- "
    set tgauss = `awk '{if ($0 ~ /mat_fit_start_image_vis2: fit gaussian /) {print $12;}}' esorex.log`
    set tchi2gauss = `awk '{if ($0 ~ /mat_fit_start_image_vis2: fit gaussian /) {print $19;}}' esorex.log`
    set tud    = `awk '{if ($0 ~ /mat_fit_start_image_vis2: fit uniform /) {print $12;}}' esorex.log`
    set tchi2ud    = `awk '{if ($0 ~ /mat_fit_start_image_vis2: fit uniform /) {print $19;}}' esorex.log`
    set tfdd   = `awk '{if ($0 ~ /mat_fit_start_image_vis2: fit fully darkened /) {print $13;}}' esorex.log`
    set tchi2fdd   = `awk '{if ($0 ~ /mat_fit_start_image_vis2: fit fully darkened /) {print $20;}}' esorex.log`
    set tld    = `awk '{if ($0 ~ /mat_fit_start_image_vis2: fit Lorentz disk /) {print $12;}}' esorex.log`
    set tchi2ld    = `awk '{if ($0 ~ /mat_fit_start_image_vis2: fit Lorentz disk /) {print $19;}}' esorex.log`
  endif
#  if($algoMode != 1) then
  if(($algoMode != 1)&&($costFunc != 3)) then
   echo "-------------- mat_cal_imarec_all.2.csh : algoMode != 1 && costFunc != 3 -------------------- "
    set tgauss = `awk '{if ($0 ~ /mat_fit_start_image_vis: fit gaussian /) {print $12;}}' esorex.log`
    set tchi2gauss = `awk '{if ($0 ~ /mat_fit_start_image_vis: fit gaussian /) {print $19;}}' esorex.log`
    set tud    = `awk '{if ($0 ~ /mat_fit_start_image_vis: fit uniform /) {print $12;}}' esorex.log`
    set tchi2ud    = `awk '{if ($0 ~ /mat_fit_start_image_vis: fit uniform /) {print $19;}}' esorex.log`
    set tfdd   = `awk '{if ($0 ~ /mat_fit_start_image_vis: fit fully darkened /) {print $13;}}' esorex.log`
    set tchi2fdd   = `awk '{if ($0 ~ /mat_fit_start_image_vis: fit fully darkened /) {print $20;}}' esorex.log`
    set tld    = `awk '{if ($0 ~ /mat_fit_start_image_vis: fit Lorentz disk /) {print $12;}}' esorex.log`
    set tchi2ld    = `awk '{if ($0 ~ /mat_fit_start_image_vis: fit Lorentz disk /) {print $19;}}' esorex.log`
  endif
  set tfov   = `awk '{if ($0 ~ /mat_guess_parameters: fov   = /) {print $0;}}' esorex.log`
  set tnpix  = `awk '{if ($0 ~ /mat_guess_parameters: npix  = /) {print $0;}}' esorex.log | tr '*' x`
  set ttfov  = `awk '{if ($0 ~ /mat_guess_parameters: --fov=/) {print $0;}}' esorex.log`
  set ttnpix = `awk '{if ($0 ~ /mat_guess_parameters: --npix=/) {print $0;}}' esorex.log`
  set ttpx   = `awk '{if ($0 ~ /mat_guess_parameters: -> psize = /) {print $0;}}' esorex.log`
  echo "------------------------------------------------------------------------------------------------------"
  echo "tgauss = $tgauss"; echo "tud = $tud"; echo "tfdd = $tfdd"; echo "tld = $tld"; echo "tfov = $tfov"; echo "tnpix = $tnpix"; echo "ttfov = $ttfov"; echo "ttnpix = $ttnpix"; echo "ttpx = $ttpx"
  if( X$tgauss == "X" ) set tgauss = -1.0
  if( X$tud == "X" ) set tud = -1.0
  if( X$tfdd == "X" ) set tfdd = -1.0
  if( X$tld == "X" ) set tld = -1.0
  echo "------------------------------------------------------------------------------------------------------"
  set fits = `echo rec_*.fits`; ls -l $fits
  if($algoMode == 1) then  
    echo $fov $fits $convScale $costFunc $regFunc $oradiusStart $stepSize $oradiusNumber $muStart $muFactor $muNumber $npix $startmode $startparam $objname $model0 $tgauss $tud $tfdd $weightPower $calcVis2f0 $tld
    $SCRIPTS/IRBis.display.Mac.nt.csh \
    $fov $fits $convScale $costFunc $regFunc $oradiusStart $stepSize $oradiusNumber $muStart $muFactor $muNumber $npix $startmode $startparam $objname $model0 \
    $tgauss $tud $tfdd $weightPower $calcVis2f0 $tld
  endif
  if($algoMode == 2) then
    echo $fov $fits $convScale $costFunc $regFunc $oradiusStart $stepSize $oradiusNumber $muStart $muFactor $muNumber $npix $startmode $startparam $objname $model0 $tgauss $tud $tfdd $weightPower $calcVisf0 $tld
    $SCRIPTS/IRBis.display.Mac.nt.ft.csh \
    $fov $fits $convScale $costFunc $regFunc $oradiusStart $stepSize $oradiusNumber $muStart $muFactor $muNumber $npix $startmode $startparam $objname $model0 $tgauss $tud $tfdd $weightPower $calcVisf0 $tld
  endif
  if($algoMode == 3) then
    echo $fov $fits $convScale $costFunc $regFunc $oradiusStart $stepSize $oradiusNumber $muStart $muFactor $muNumber $npix $startmode $startparam $objname $model0 $tgauss $tud $tfdd $weightPower $calcVisf0 $tld
    $SCRIPTS/IRBis.display.Mac.nt.bisft.csh \
    $fov $fits $convScale $costFunc $regFunc $oradiusStart $stepSize $oradiusNumber $muStart $muFactor $muNumber $npix $startmode $startparam $objname $model0 $tgauss $tud $tfdd $weightPower $calcVisf0 $tld
  endif
#  echo "weiter mit CR"
#  set ja = $<
cd ..
mv ${results}/RR.konv.iterated.rek1.lin.ps ${results}/E.lin.ps; mv ${results}/RR.konv.iterated.rek1.sqrt.ps ${results}/E.sqrt.ps; mv ${results}/RR.konv.iterated.rek1.log.ps ${results}/E.log.ps
ls -l  ${results}/E.lin.ps ${results}/E.sqrt.ps ${results}/E.log.ps ${results}/liste
# removing all *.ps, but not E.lin.ps, E.sqrt.ps, E.log.ps
cd $results; $SCRIPTS/delps.csh; cd ..
# -----------
if($algoMode == 1) then
  set chi2cpc = `awk '{ if($1!="#") {print $7;}; }' ${results}/liste`
  set rrescp  = `awk '{ if($1!="#") {print $8;}; }' ${results}/liste`
  set cpqrec  = `echo $chi2cpc $rrescp | awk '{ q=(($1-1.)+($2-1.))/2.; if($1<1.) {q=(((1./$1)-1.)+($2-1.))/2.;}; print q;}'`
  set phiqrec  = -1.
  set cpphiqrec = -1.
endif
if($algoMode == 2) then
  set chi2phi  = `awk '{ if($1!="#") {print $7;}; }' ${results}/liste`
  set rresphi  = `awk '{ if($1!="#") {print $8;}; }' ${results}/liste`
  set phiqrec  = `echo $chi2phi $rresphi | awk '{ q=(($1-1.)+($2-1.))/2.; if($1<1.) {q=(((1./$1)-1.)+($2-1.))/2.;}; print q;}'`
  set cpqrec   = -1.
  set cpphiqrec = -1.
endif
if($algoMode == 3) then
  set chi2phi  = `awk '{ if($1!="#") {print $7;}; }' ${results}/liste`
  set rresphi  = `awk '{ if($1!="#") {print $8;}; }' ${results}/liste`
  set chi2cpc = `awk '{ if($1!="#") {print $3;}; }' ${results}/liste`
  set rrescp  = `awk '{ if($1!="#") {print $4;}; }' ${results}/liste`
  set cpphiqrec = `echo $chi2cpc $rrescp $chi2phi $rresphi | awk '{ cb=$1; if($1<1.) {cb=1./$1;}; cf=$3; if($3<1.) {cf=1./$3;}; q=((cb-1.)+($2-1.)+(cf-1.)+($4-1.))/4.; print q;}'`
  set phiqrec  = -1.
  set cpqrec   = -1.
endif
if($qrecmode == 1) set qrec = `awk '{ if($1!="#") {print $1;}; }' ${results}/liste`
if(($qrecmode == 2)&&($algoMode == 1)) set qrec = $cpqrec
if(($qrecmode == 2)&&($algoMode == 2)) set qrec = $phiqrec
if(($qrecmode == 2)&&($algoMode == 3)) set qrec = $cpphiqrec
set weiter = `echo $qrec $qrecvorher | awk '{ weiter=0; if($1<$2) {weiter=1;}; print weiter; }'`
set qrecvorher = $qrec
if( $weiter == 0 ) goto Label222
set muStart1 = `echo $muStart $muFactor0 | awk '{print $1*$2;}'`
set muStart  = $muStart1
# -----------
echo "Ergebnisse stehen in $results in $resultdir"
cat ${results}/liste
# gv ${results}/RR.konv.iterated.rek1.lin.ps
# ------------ loop over chosen start values of mu ----------------------------------- E ------------
end
Label222:
# ------------ loop over chosen regularization functions ------------------------------- E --------------
end

set endtime = `date`
set host = `hostname`
rm -f runtime.txt; echo "run from day $starttime[3] $starttime[4] until day $endtime[3] $endtime[4] on $host" >> runtime.txt

set dir0s = `echo E.*`
foreach f ($dir0s)
   mv $f $resultdir
end
mv runtime.txt $resultdir
# cp $0 $resultdir
cp $parfile $resultdir
cd $resultdir
   set listen = `echo E.*/liste`
   rm -f E.liste.0; cat $listen >> E.liste.0
   awk '{ if($1!="#") {print $0;}; }' E.liste.0 > E.liste.00

   rm -f head; awk 'BEGIN{z=1;} { if((z<4)&&($1=="#")) {print $0;}; z=z+1; }' E.liste.0 > head
   if($algoMode == 1) rm -f head00; awk 'BEGIN{z=1;} { if(z==3) {print $0,"| cpqrec";}; z=z+1; }' head > head00
   if($algoMode == 2) rm -f head00; awk 'BEGIN{z=1;} { if(z==3) {print $0,"| phqrec";}; z=z+1; }' head > head00
   if($algoMode == 3) rm -f head00; awk 'BEGIN{z=1;} { if(z==3) {print $0,"| cpphqrec";}; z=z+1; }' head > head00
   rm -f E.liste.00.sorted; sort -k 1 -n E.liste.00 > E.liste.00.sorted
   if(($algoMode == 1)||($algoMode == 2)) then
     rm -f E.liste.03; awk '{ cb=$7; if($7<1.) {cb=1./$7;}; q=( (cb-1.)+($8-1.) )/2.; print $0,q;}' E.liste.00 > E.liste.03
   endif
   if($algoMode == 3) then
     rm -f E.liste.03; awk '{ cb=$3; if($3<1.) {cb=1./$3;}; cf=$7; if($7<1.) {cf=1./$7;}; q=( (cb-1.)+($4-1.)+(cf-1.)+($8-1.) )/4.; print $0,q;}' E.liste.00 > E.liste.03
   endif
   rm -f E.liste.03.sorted; sort -k 24 -n E.liste.03 > E.liste.03.sorted
   if( $model0 != "no" ) then
#       qrec built from all measuruments ($qrecmode == 1)
        rm -f E.liste.00.sorted.plot; awk '{ if($9!=0.) {print $0;}; }' E.liste.00.sorted > E.liste.00.sorted.plot
	set distqrecplot = $resultdir.distqrec.1.plot.ps
	rm -f $distqrecplot; $SCRIPTS/gplot -psfl -o $distqrecplot -pscolor -headline "$name" \
                               -xlabel "qrec" -ylabel "distance to theoretical object (khh)" -grid -x 1 -y 9 -t "Distance on qrec: qrec from all measurements" E.liste.00.sorted.plot
#       qrec built from the phase measurements only ($qrecmode == 2)
        rm -f E.liste.03.sorted.plot; awk '{ if($9!=0.) {print $0;}; }' E.liste.03.sorted > E.liste.03.sorted.plot
	set distqrecplot = $resultdir.distqrec.2.plot.ps
	rm -f $distqrecplot; $SCRIPTS/gplot -psfl -o $distqrecplot -pscolor -headline "$name" \
			       -xlabel "qrec" -ylabel "distance to theoretical object (khh)" -grid -x 24 -y 9 -t "Distance on qrec: qrec from phase measurements only" E.liste.03.sorted.plot
	rm -f E.liste.03.sorted.plot E.liste.00.sorted.plot
   endif
   if($algoMode == 1) then
     rm -f head03; echo "#  " >> head03; echo "# 2) reconstructions sorted with increasing CP-qrec (CP-qrec in the last column):" >> head03
   endif
   if($algoMode == 2) then
     rm -f head03; echo "#  " >> head03; echo "# 2) reconstructions sorted with increasing PH-qrec (PH-qrec in the last column):" >> head03
   endif
   if($algoMode == 3) then
     rm -f head03; echo "#  " >> head03; echo "# 2) reconstructions sorted with increasing CP-PH-qrec (CP-PH-qrec in the last column):" >> head03
   endif
   rm -f head04; echo "#  " >> head04; echo "# 1) reconstructions sorted with increasing qrec:" >> head04
   rm -f head05; echo "#  " >> head05
   rm -f E.liste; cat head head04 E.liste.00.sorted head05 head00 head03 E.liste.03.sorted > E.liste
   set dir0qrec   = `awk 'BEGIN{z=0;} { z=z+1; if(z==1) {print $23;}; }' E.liste.00.sorted`
   set dir0cpqrec = `awk 'BEGIN{z=0;} { z=z+1; if(z==1) {print $23;}; }' E.liste.03.sorted`
#  set fits = `echo rec_*.fits`; ls -l $fits
   cp $dir0qrec/rec_*.fits   bestrec.qrec1.fits
   cp $dir0cpqrec/rec_*.fits bestrec.qrec2.fits
   cp $SCRIPTS/displ.csh .
#
   rm -rf qrecmode=1 qrecmode=2; mkdir qrecmode=1 qrecmode=2
   ./displ.csh E.liste.00.sorted
   mv $cwd:t.lin.ps $cwd:t.sqrt.ps $cwd:t.log.ps qrecmode=1/
   ./displ.csh E.liste.03.sorted
   mv $cwd:t.lin.ps $cwd:t.sqrt.ps $cwd:t.log.ps qrecmode=2/
#   if( $qrecmode == 1 ) ./displ.csh E.liste.00.sorted
#   if( $qrecmode == 2 ) ./displ.csh E.liste.03.sorted
   rm -f head03 head04 E.liste.00 E.liste.0
cd ..
# echo "- some proposed image reconstruction parameter are listed in data.parameter "
# ls -l $resultdir/data.parameter
echo "=============================================================================================================================================================================="
echo "---> chi^2 values and more of all runs are listed in E.liste "
ls -l $resultdir/E.liste
echo "------------------------------------------------------------------------------------------------------------------------------------------------"
# if($qrecmode == 2) echo "- the reconstructions of all runs are displayed sorted with according decreasing quality (using  phase qrec)"
# if($qrecmode == 1) echo "- the reconstructions of all runs are displayed sorted with according decreasing quality (using phase+visibility qrec)"
echo "   a) with a linear lookuptable in $resultdir.lin.ps "
echo "   b) with a sqrt   lookuptable in $resultdir.sqrt.ps "
echo "   c) with a log    lookuptable in $resultdir.log.ps "
echo "------------------------------------------------------------------------------------------------------------------------------------------------"
echo "---> best reconstruction according to qrec from qrecmode=1: bestrec.qrec1.fits "
ls -l $resultdir/bestrec.qrec1.fits
echo "---> best reconstruction according to qrec from qrecmode=2: bestrec.qrec2.fits "
ls -l $resultdir/bestrec.qrec2.fits
echo "------------------------------------------------------------------------------------------------------------------------------------------------"
if($model0 != "no") echo "---> Distance reconstructed & theoretical object on qrec in $resultdir.distqrec.1.plot.ps and $resultdir.distqrec.2.plot.ps"
if($model0 != "no") ls -l $resultdir/$resultdir.distqrec.1.plot.ps $resultdir/$resultdir.distqrec.2.plot.ps
echo "- the reconstructions of all runs are displayed sorted according to decreasing quality (using  phase qrec):"
ls -l $resultdir/qrecmode=2/$resultdir.lin.ps $resultdir/qrecmode=2/$resultdir.sqrt.ps $resultdir/qrecmode=2/$resultdir.log.ps
echo "- the reconstructions of all runs are displayed sorted according to decreasing quality (using phase+visibility qrec):"
ls -l $resultdir/qrecmode=1/$resultdir.lin.ps $resultdir/qrecmode=1/$resultdir.sqrt.ps $resultdir/qrecmode=1/$resultdir.log.ps
echo "=============================================================================================================================================================================="

Label999:
