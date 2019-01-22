#!/bin/tcsh
#
set marke = $#argv
if($marke == 0) then
  echo "Aufgabe: Bildrekonstruktion mit ESO-IRBis: a) Erzeugung der Plotts, welche die gemessenen"
  echo "         Visbilities/FT phases vergleichen mit den aus den 3 besten Rekonstruktionen abgeleiteten;"
  echo "         b) Abbildung der besten Rek. (ungefaltet/gefaltet) und das Original, falls vorhanden, wenn das Original"
  echo "            nicht vorhanden, dann werden das Startbild & der Start-Prior dargestellt."
  echo
  echo "Eingaben: [FOV of reconstruction in mas] [result fits file from ESO-IRBis] [scaling factor for the convolution of the reconstruction]"
  echo "          [cost funcion, regularization function, oradius,step,count, mu,mufac,count, npix, startmode,startparam] [OBJ-Name] [model image/no]"
  echo "          [FWHM of fitted Gaussian (mas)] [diameter of fitted UD (mas)] [diameter of fitted FDD (mas)] [power for the UV density weight]"
  echo "          [1: for frequency 0 a visibility is calculated, 0: this is not done] [FWHM of the fitted Lorentz function (mas)]"
  echo "Ausgaben: Visibility-Plotts: visa.ps visb.ps visc.ps; FT-Phases-Plotts: cpa.ps cpb.ps cpc.ps"

else

set fov  = $1  # fov in mas
set fits = $2  # Ergebniss-Fits-File aus ESO-IRBis
set convscale = $3  # Scale factor for the convolution (1.0 means telescope diameter = max basline; >1.0 means telescope diameter > max basline)
set costFunc = $4
set regFunc  = $5
set oradiusStart = $6
set stepSize     = $7
set oradiusNumber = $8
set muStart  = $9
set muFactor = $10
set muNumber = $11
set npix = $12
set startmode = $13
set startparam = $14
set objname = $15
set model = $16
set tgauss = $17
set tud    = $18
set tfdd   = $19
set weightPower = $20
set calcVisf0 = $21
set tld    = $22

echo "Inputs: "
echo "fov fits convscale costFunc regFunc oradiusStart stepSize oradiusNumber muStart muFactor muNumber npix startmode startparam objname model tgauss tud tfdd weightPower calcVisf0: "
echo $fov $fits $convscale $costFunc $regFunc $oradiusStart $stepSize $oradiusNumber $muStart $muFactor $muNumber $npix $startmode $startparam $objname $model $tgauss $tud $tfdd $weightPower $calcVisf0
echo "--- IRBis.display.Mac.nt.ft.csh --------------------------------------------------------------"

echo "a) Erzeugung der Plotts: Visbilities/ft phases"
set ft0 = `echo *.ft.dat`
ls -l $ft0
# Umwandeln in altes Format:
set ft = $ft0.new; rm -f $ft; awk '{ if($1!="#") {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$14,$15,$16,$17,$18,$19,$13;}; }' $ft0 > $ft
echo "FOV = $fov mas"
#
# ft - new (11-12-15):
# measured complex visibilities vs reconstructed complex visibilities
# nr lambda u v bl amp amperr phi phierr gd ud fdd ld rec1_amp rec1_phi rec2_amp rec2_phi rec3_amp rec3_phi
# 1  2      3 4 5  6   7      8   9      10 11 12  13 14       15       16       17       18       19
#
# ft - old:
# measured complex visibilities vs reconstructed complex visibilities
# nr lambda u v bl amp amperr phi phierr gd ud fdd rec1_amp rec1_phi rec2_amp rec2_phi rec3_amp rec3_phi
# 1  2      3 4 5  6   7      8   9      10 11 12  13       14       15       16       17       18
#
# $ft0.new  == $ft
# nr lambda u v bl amp amperr phi phierr gd ud fdd rec1_amp rec1_phi rec2_amp rec2_phi rec3_amp rec3_phi ld
# 1  2      3 4 5  6   7      8   9      10 11 12  13       14       15       16       17       18       19
#
# vis2:
# measured squared visibilities vs reconstructed squared visibilities
# nr lambda u v bl vis2 err gd ud fdd rec1 rec2 rec3
# 1  2      3 4 5  6    7   8  9  10  11   12   13
# cp: neu seit 13-11-15:
# nr lambda u1 v1 u2 v2 amp amperr cp cperr biserr rec1_amp rec1_cp rec2_amp rec2_cp rec3_amp rec3_cp
# 1  2      3  4  5  6  7   8      9  10    11     12       13       14      15      16       17
rm -f $ft.vis.0 $ft.ph.0

set chivisa = `awk 'BEGIN{sum=0.;z=0.} { if($1!="#") {res=($6-$13)/$7; sum=sum+res^2; z=z+1.;}; } END{printf "%8.3f\n", sum/z;}' $ft`
set chivisb = `awk 'BEGIN{sum=0.;z=0.} { if($1!="#") {res=($6-$15)/$7; sum=sum+res^2; z=z+1.;}; } END{printf "%8.3f\n", sum/z;}' $ft`
set chivisc = `awk 'BEGIN{sum=0.;z=0.} { if($1!="#") {res=($6-$17)/$7; sum=sum+res^2; z=z+1.;}; } END{printf "%8.3f\n", sum/z;}' $ft`
set resratioa = `awk 'BEGIN{sump=0.;summ=0.;} { if($1!="#") {res=($6-$13)/$7; if(res>0.) {sump=sump+res;}; if(res<=0.) {summ=summ-res;}; }; } END{rr=10000.; if((sump>0.)&&(summ>0.)) {rr=sump/summ; if(rr<1.) {rr=summ/sump;};}; printf "%8.3f\n",rr;}' $ft`
set resratiob = `awk 'BEGIN{sump=0.;summ=0.;} { if($1!="#") {res=($6-$15)/$7; if(res>0.) {sump=sump+res;}; if(res<=0.) {summ=summ-res;}; }; } END{rr=10000.; if((sump>0.)&&(summ>0.)) {rr=sump/summ; if(rr<1.) {rr=summ/sump;};}; printf "%8.3f\n",rr;}' $ft`
set resratioc = `awk 'BEGIN{sump=0.;summ=0.;} { if($1!="#") {res=($6-$17)/$7; if(res>0.) {sump=sump+res;}; if(res<=0.) {summ=summ-res;}; }; } END{rr=10000.; if((sump>0.)&&(summ>0.)) {rr=sump/summ; if(rr<1.) {rr=summ/sump;};}; printf "%8.3f\n",rr;}' $ft`

awk -v fov=$fov 'BEGIN{z=0;} { if(($1=="#")&&(z==0)) {print "# Nr. |", "U |", "V |", "Vis |", "VisErr |", "Gaussian |", "Uniform d |", "Fully dark d |", "Vis (1.rec) |", "Vis (2.rec) |", "Vis (3.rec) |", "f (1/arcsec) |", "Lorentz ";}; \
                               if(($1=="#")&&(z==1)) {print "# 1   |", "2 |", "3 |", "4   |", "5      |", "6        |", "7        |", "8            |", "9           |", "10          |", "11          |", "12           |", "13  ";}; \
                               if($1!="#") {rad=sqrt($3^2+$4^2); radc=rad*1000.; vis=sqrt($6^2); if($6==0.0) {viserr=0.0;}; if($6!=0.0) {viserr=$7;}; \
                                            visgauss=sqrt($10^2); visud=sqrt($11^2); visfdd=sqrt($12^2); visld=sqrt($19^2); visa=sqrt($13^2); visb=sqrt($15^2); visc=sqrt($17^2); \
                                            print $1,$3,$4,vis,viserr,visgauss,visud,visfdd, visa,visb,visc,radc,visld;}; z=z+1;}' $ft > $ft.vis.0


set pi = `echo 1. | awk '{print atan2($1,$1)*4.;}'`
set chi2fta   = `awk -v p=$pi 'BEGIN{sum=0.;z=0.} { if($1!="#") {d0=$8-$14;d1=d0+2.*p;d2=d0-2.*p;d00=d0^2;d11=d1^2;d22=d2^2;d=d0;dd=d00;if(dd>d11) {d=d1;dd=d11;};if(dd>d22) {d=d2;dd=d22};res=d/$9; sum=sum+res^2; z=z+1.;}; } END{printf "%8.3f\n", sum/z;}' $ft`
set chi2ftb   = `awk -v p=$pi 'BEGIN{sum=0.;z=0.} { if($1!="#") {d0=$8-$16;d1=d0+2.*p;d2=d0-2.*p;d00=d0^2;d11=d1^2;d22=d2^2;d=d0;dd=d00;if(dd>d11) {d=d1;dd=d11;};if(dd>d22) {d=d2;dd=d22};res=d/$9; sum=sum+res^2; z=z+1.;}; } END{printf "%8.3f\n", sum/z;}' $ft`
set chi2ftc   = `awk -v p=$pi 'BEGIN{sum=0.;z=0.} { if($1!="#") {d0=$8-$18;d1=d0+2.*p;d2=d0-2.*p;d00=d0^2;d11=d1^2;d22=d2^2;d=d0;dd=d00;if(dd>d11) {d=d1;dd=d11;};if(dd>d22) {d=d2;dd=d22};res=d/$9; sum=sum+res^2; z=z+1.;}; } END{printf "%8.3f\n", sum/z;}' $ft`
set ftresa   = `awk -v p=$pi 'BEGIN{sump=0.;summ=0.;} { if($1!="#") {d0=$8-$14;d1=d0+2.*p;d2=d0-2.*p;d00=d0^2;d11=d1^2;d22=d2^2;d=d0;dd=d00;if(dd>d11) {d=d1;dd=d11;};if(dd>d22) {d=d2;dd=d22};res=d/$9; if(res>0.) {sump=sump+res;}; if(res<=0.) {summ=summ-res;}; }; } END{rr=10000.; if((sump>0.)&&(summ>0.)) {rr=sump/summ; if(rr<1.) {rr=summ/sump;};}; printf "%8.3f\n",rr;}' $ft`
set ftresb   = `awk -v p=$pi 'BEGIN{sump=0.;summ=0.;} { if($1!="#") {d0=$8-$16;d1=d0+2.*p;d2=d0-2.*p;d00=d0^2;d11=d1^2;d22=d2^2;d=d0;dd=d00;if(dd>d11) {d=d1;dd=d11;};if(dd>d22) {d=d2;dd=d22};res=d/$9; if(res>0.) {sump=sump+res;}; if(res<=0.) {summ=summ-res;}; }; } END{rr=10000.; if((sump>0.)&&(summ>0.)) {rr=sump/summ; if(rr<1.) {rr=summ/sump;};}; printf "%8.3f\n",rr;}' $ft`
set ftresc   = `awk -v p=$pi 'BEGIN{sump=0.;summ=0.;} { if($1!="#") {d0=$8-$18;d1=d0+2.*p;d2=d0-2.*p;d00=d0^2;d11=d1^2;d22=d2^2;d=d0;dd=d00;if(dd>d11) {d=d1;dd=d11;};if(dd>d22) {d=d2;dd=d22};res=d/$9; if(res>0.) {sump=sump+res;}; if(res<=0.) {summ=summ-res;}; }; } END{rr=10000.; if((sump>0.)&&(summ>0.)) {rr=sump/summ; if(rr<1.) {rr=summ/sump;};}; printf "%8.3f\n",rr;}' $ft`

# qrec-Berechnung aus den obigen direkt berechneten Chi^2-/RR-Werten
# mit qrec = ( |Cv-1|+|RR{V^2}-1| + |Cb-1|+|RR{CP}-1| )/4; mit Cv=Chi^2{V^2} oder Cv=1/Chi^2{V^2} bei Chi^2{V^2}<1, Cb=Chi^2{CP} oder Cb=1/Chi^2{CP} bei Chi^2{CP}<1
set qreca = `echo $chivisa $resratioa $chi2fta $ftresa | awk '{Cv=$1; if(Cv<1.) {Cv=1./$1;}; Cb=$3; if(Cb<1.) {Cb=1./$3;}; qrec = ( sqrt((Cv-1.)^2)+sqrt(($2-1.)^2)+sqrt((Cb-1.)^2)+sqrt(($4-1.)^2) )/4.; printf "%12.3f", qrec;}'`
echo $qreca
set phqreca = `echo $chivisa $resratioa $chi2fta $ftresa | awk '{Cv=$1; if(Cv<1.) {Cv=1./$1;}; Cb=$3; if(Cb<1.) {Cb=1./$3;}; qrec = ( sqrt((Cb-1.)^2)+sqrt(($4-1.)^2) )/2.; printf "%12.3f", qrec;}'`
echo $phqreca


# ft:
# measured complex visibilities vs reconstructed complex visibilities
# nr lambda u v bl amp amperr phi phierr gd ud fdd rec1_amp rec1_phi rec2_amp rec2_phi rec3_amp rec3_phi
# 1  2      3 4 5  6   7      8   9      10 11 12  13       14       15       16       17       18
# fuer Plott aufbereiten: die rek. PH durch Addition von +-360Grad anpassen a die gemessenen PHs (moeglichst kleiner Abstand)
set pi = `echo 1. | awk '{print atan2($1,$1)*4.;}'`
awk -v p=$pi 'BEGIN{z=0;} { if(($1=="#")&&(z==0)) {print "# Nr. |", "U |", "V |", "PH [deg]|", "PH err |", "PH (1.rec) |", "PH (2.rec) |", "PH (3.rec) |"}; \
                             if(($1=="#")&&(z==1)) {print "# 1   |", "2  |", "3  |", "4       |", "5      |", "6          |", "7          |", "8         |";}; \
                             if($1!="#") {d0=$8-$14;d1=d0+2.*p;d2=d0-2.*p;d00=d0^2;d11=d1^2;d22=d2^2;d8=$14;dd=d00;if(dd>d11) {d8=$14-2.*p;dd=d11;};if(dd>d22) {d8=$14+2.*p;dd=d22}; \
                                          d0=$8-$16;d1=d0+2.*p;d2=d0-2.*p;d00=d0^2;d11=d1^2;d22=d2^2;d9=$16;dd=d00;if(dd>d11) {d9=$16-2.*p;dd=d11;};if(dd>d22) {d9=$16+2.*p;dd=d22}; \
                                          d0=$8-$18;d1=d0+2.*p;d2=d0-2.*p;d00=d0^2;d11=d1^2;d22=d2^2;d10=$18;dd=d00;if(dd>d11) {d10=$18-2.*p;dd=d11;};if(dd>d22) {d10=$18+2.*p;dd=d22}; \
                                          ph=$8*180./p; pherr=$9*180./p; pha=d8*180./p; phb=d9*180./p; phc=d10*180./p; print $1,$3,$4,ph,pherr,pha,phb,phc;}; z=z+1;}' $ft > $ft.ph.0

# A) plott of the visibilities as function of the length of the spatial frequency
rm -f visa.ps visb.ps visc.ps
rm -f $ft.vis.0.sorted ttt; awk '{if($1!="#") {print $0;};}' $ft.vis.0 > ttt; sort -n -k 12 ttt > $ft.vis.0.sorted

ls -l {$ft}*

set psplotva = visa.ps
rm -f $psplotva
# ./gnupl.csh
# Input: [name of the ascii file of the data] [text for the measured data] [text for the fitted data] 
#        [column of the X axis in the file] [column of the measured data] [column of the errors of the measured data]
#        [column of the fitted data] [text on X axis] [text on Y axis] [name of the ps plot]
$SCRIPTS/gnupl2.csh $ft.vis.0.sorted "measured visibilities" 12 4 5 "spatial frequency (1/arcsec)" "Visibility" $psplotva "Chi\^2/RR = $chivisa/$resratioa" 9
ls -l $psplotva

set psplotvb = visb.ps
rm -f $psplotvb
$SCRIPTS/gnupl2.csh $ft.vis.0.sorted "measured visibilities" 12 4 5 "spatial frequency (1/arcsec)" "Visibility" $psplotvb "Chi\^2/RR = $chivisb/$resratiob" 10
ls -l $psplotva $psplotvb

set psplotvc = visc.ps
rm -f $psplotvc
$SCRIPTS/gnupl2.csh $ft.vis.0.sorted "measured visibilities" 12 4 5 "spatial frequency (1/arcsec)" "Visibility" $psplotvc "Chi\^2/RR = $chivisc/$resratioc" 11

ls -l $psplotva $psplotvb $psplotvc
rm -f t1.eps t2.eps; cp $psplotva t1.eps; cp $psplotvb t2.eps

# ./gnupl2.csh
# Input: [name of the ascii file of the data] [text for the measured data]
#        [column of the X axis in the file] [column of the measured data] [column of the errors of the measured data]
#        [text on X axis] [text on Y axis] [name of the ps plot]
#        [text for the fitted data] [column of the fitted data] ...

# A2) plot of the fitted model visibilities and measured visibilities as a function of the length of the spatial frequency
rm -f gaussa.ps uda.ps fdda.ps lda.ps gaussudfdda.ps
$SCRIPTS/gnupl2.csh $ft.vis.0.sorted "measured visibilities" 12 4 5 "spatial frequency (1/arcsec)" "Visibility" gaussa.ps "fitted Gaussian model (FWHM={$tgauss}mas)" 6

$SCRIPTS/gnupl2.csh $ft.vis.0.sorted "measured visibilities" 12 4 5 "spatial frequency (1/arcsec)" "Visibility" uda.ps "fitted Uniform disk model (diameter={$tud}mas)" 7

$SCRIPTS/gnupl2.csh $ft.vis.0.sorted "measured visibilities" 12 4 5 "spatial frequency (1/arcsec)" "Visibility" fdda.ps "fitted Fully darkened disk model (diameter={$tfdd}mas)" 8

$SCRIPTS/gnupl2.csh $ft.vis.0.sorted "measured visibilities" 12 4 5 "spatial frequency (1/arcsec)" "Visibility" lda.ps "fitted Lorentz function (FWHM={$tld}mas)" 13


$SCRIPTS/gnupl2.csh $ft.vis.0.sorted "measured visibilities" 12 4 5 "spatial frequency (1/arcsec)" "Visibility" gaussudfdda.ps \
			    "fitted Uniform disk model (diameter={$tud}mas)" 7 "fitted Gaussian model (FWHM={$tgauss}mas)" 6 "fitted Fully darkened disk model (diameter={$tfdd}mas)" 8 \
			    "fitted Lorentz function (FWHM={$tld}mas)" 13

ls -l gaussa.ps uda.ps fdda.ps lda.ps gaussudfdda.ps

# $ft.ph.0 :
# Nr. |", "U |", "V |", "PH [deg]|", "PH err |", "PH (1.rec) |", "PH (2.rec) |", "PH (3.rec) |
# 1   |", "2 |", "3 |", "4       |", "5      |", "6          |", "7          |", "8          |
# C) plott of the closure phases ordered according to increasing values of the measured closure phases
# $cp.0 :
#  Nr. |", "U1 |", "V1 |", "U2 |", "V2 |", "CP [deg]|", "CP err |", "CP (1.rec) |", "CP (2.rec) |", "CP (3.rec) |
#  1   |", "2  |", "3  |", "4  |", "5  |", "6       |", "7      |", "8          |", "9          |", "10         |
# a) sorting the PHs according to increasing values of the measuered PHs --> $ft.ph.0.sorted.0
rm -f $ft.ph.0.sorted ttt; awk '{if($1!="#") {print $0;};}' $ft.ph.0 > ttt; sort -n -k 4 ttt > $ft.ph.0.sorted
rm -f $ft.ph.0.sorted.0; awk 'BEGIN{z=0;} { z=z+1; print $0, z; }' $ft.ph.0.sorted > $ft.ph.0.sorted.0

# b) plotting the CPs as function of the longest baseline vector:
# $cp.0.lbl :
# length of longest baseline | measured CP | CP error | CP (1.rec) | CP (2.rec) | CP (3.rec) |
# 1                            2             3          4            5            6
# b) plotting the PHs as function of the length of the baseline vector:
# $ft.ph.0.l :
# length of baseline | measured PH | PH error | PH (1.rec) | PH (2.rec) | PH (3.rec) |
# 1                        2             3          4            5            6
rm -f $ft.ph.0.l $ft.ph.0.l.sorted
awk '{ if($1!="#") {ll=sqrt($2^2+$3^2); lll=ll*1000.; print lll, $4, $5. $6, $7, $8;}; }' $ft.ph.0 > $ft.ph.0.l
sort -n -k 1 $ft.ph.0.l > $ft.ph.0.l.sorted


# ./gnupl2.csh
# Input: [name of the ascii file of the data] [text for the measured data]
#        [column of the X axis in the file] [column of the measured data] [column of the errors of the measured data]
#        [text on X axis] [text on Y axis] [name of the ps plot]
#        [text for the fitted data] [column of the fitted data] ...

# to a)
rm -f cpa.ps
$SCRIPTS/gnupl2.csh $ft.ph.0.sorted.0 "measured ft phases" 9 4 5 "number of increasing measured ft phase" "FT phase (deg)" cpa.ps \
"Chi^2/RR = $chi2fta/$ftresa" 6
rm -f cpb.ps
$SCRIPTS/gnupl2.csh $ft.ph.0.sorted.0 "measured ft phases" 9 4 5 "number of increasing measured ft phase" "FT phase (deg)" cpb.ps \
"Chi^2/RR = $chi2ftb/$ftresb" 7
rm -f cpc.ps
$SCRIPTS/gnupl2.csh $ft.ph.0.sorted.0 "measured ft phases" 9 4 5 "number of increasing measured ft phase" "FT phase (deg)" cpc.ps \
"Chi^2/RR = $chi2ftc/$ftresc" 8

# to b)
rm -f cpaa.ps
$SCRIPTS/gnupl3.csh $ft.ph.0.l.sorted "measured ft phases" 1 2 3 "spatial frequency of baseline length (1/arcsec)" "FT phase (deg)" cpaa.ps \
"Chi^2/RR = $chi2fta/$ftresa" 4
rm -f cpbb.ps
$SCRIPTS/gnupl3.csh $ft.ph.0.l.sorted "measured ft phases" 1 2 3 "spatial frequency of baseline length (1/arcsec)" "FT phase (deg)" cpbb.ps \
"Chi^2/RR = $chi2ftb/$ftresb" 5
rm -f cpcc.ps
$SCRIPTS/gnupl3.csh $ft.ph.0.l.sorted "measured ft phases" 1 2 3 "spatial frequency of baseline length (1/arcsec)" "FT phase (deg)" cpcc.ps \
"Chi^2/RR = $chi2ftc/$ftresc" 6

ls -l cpa.ps cpb.ps cpc.ps cpaa.ps cpbb.ps cpcc.ps


cp $SCRIPTS/plot3b.tex .
rm -f plot3b.ps results.1.reconstruction.ps; rm -f t1.eps t2.eps t3.eps; cp gaussudfdda.ps t1.eps; cp visa.ps t2.eps; cp cpa.ps t3.eps; latex plot3b.tex; dvips plot3b.dvi; mv plot3b.ps results.1.reconstruction.ps
rm -f plot3b.ps results.2.reconstruction.ps; rm -f t1.eps t2.eps t3.eps; cp gaussudfdda.ps t1.eps; cp visb.ps t2.eps; cp cpb.ps t3.eps; latex plot3b.tex; dvips plot3b.dvi; mv plot3b.ps results.2.reconstruction.ps
rm -f plot3b.ps results.3.reconstruction.ps; rm -f t1.eps t2.eps t3.eps; cp gaussudfdda.ps t1.eps; cp visc.ps t2.eps; cp cpc.ps t3.eps; latex plot3b.tex; dvips plot3b.dvi; mv plot3b.ps results.3.reconstruction.ps
rm -f plot3b.*

ls -l results.1.reconstruction.ps results.2.reconstruction.ps results.3.reconstruction.ps


echo "b) Abbildung der besten Rek. (ungefaltet/gefaltet) und das Original oder das Startbild & Prior"

	echo "--------------- fdump --------------------- A ---------------------------"

	# quality parameter of the best reconstruction:
         fdump  $fits\[2\] \!rec33.txt QREC -; set qrec = `awk '{if ($0 ~ /  1  /) {printf "%12.4f", $2;}}' rec33.txt`
         fdump  $fits\[2\] \!rec33.txt COST -; set cost = `awk '{if ($0 ~ /  1  /) {printf "%12.4f", $2;}}' rec33.txt`
	 fdump  $fits\[2\] \!rec33.txt CHI2VIS -; set chi2vis = `awk '{if ($0 ~ /  1  /) {printf "%12.4f", $2;}}' rec33.txt`
         fdump  $fits\[2\] \!rec33.txt RRESVIS -; set rresvis = `awk '{if ($0 ~ /  1  /) {printf "%12.4f", $2;}}' rec33.txt`
         fdump  $fits\[2\] \!rec33.txt CHI2AMP -; set chi2amp = `awk '{if ($0 ~ /  1  /) {printf "%12.4f", $2;}}' rec33.txt`
         fdump  $fits\[2\] \!rec33.txt RRESAMP -; set rresamp = `awk '{if ($0 ~ /  1  /) {printf "%12.4f", $2;}}' rec33.txt`
         fdump  $fits\[2\] \!rec33.txt CHI2PHI -; set chi2phi = `awk '{if ($0 ~ /  1  /) {printf "%12.4f", $2;}}' rec33.txt`
         fdump  $fits\[2\] \!rec33.txt RRESPHI -; set rresphi = `awk '{if ($0 ~ /  1  /) {printf "%12.4f", $2;}}' rec33.txt`

        set phqrec = `echo $chi2phi $rresphi | awk '{ Cb=$1; if($1<1.) {Cb=1./$1;}; qrec = ( sqrt((Cb-1.)^2) + sqrt(($2-1.)^2) )/2.; printf "%12.4f", qrec; }'`
	echo "# qrec   |  cost  |  chi2vis |  rresvis | chi2amp | rresamp | chi2phi | rresphi | cpqrec "
	echo  $qrec "  " $cost "  " $chi2vis "   " $rresvis "  " $chi2amp "   " $rresamp "  " $chi2phi "  " $rresphi "  " $phqrec
	echo "--------------- fdump --------------------- E ---------------------------"
        

#        rm -f tkhh; awk -v w=$qrec '{ if(($15 == w",")||($15 == w"0,")||($15 == w"00,")) {print $0;}; }' esorex.log > tkhh
#        set khh0       = `awk '{w=$45;} END{printf "%12.4f", w;}' tkhh`              ; if( x$khh0 == "x" ) set khh = -1; if( x$khh0 != "x" ) set khh = $khh0; echo $khh; rm -f tkhh
        set qrecbest = `awk -v w=$qrec 'BEGIN{mindist=1e10;} { dist=sqrt((w-$15)^2); if(dist < mindist) {mindist=dist; qrecbest=$15;}; } END{print qrecbest;}' esorex.log`
        rm -f tkhh; awk -v w=$qrecbest '{ if($15 == w) {print $0;}; }' esorex.log > tkhh
        set khh0       = `awk 'BEGIN{z=0;} {z=z+1; if(z==1) {w=$45;};} END{printf "%12.3f", w;}' tkhh`              ; if( x$khh0 == "x" ) set khh = -1; if( x$khh0 != "x" ) set khh = $khh0; echo $khh; rm -f tkhh
	set khhm = $khh

	rm -f liste0;  fstruct colinfo=no $fits >> liste0
	set bestrecNr = 0; set bestrecconvNr = `awk '{if ($0 ~ /REC_CONV /) {print $1;}}' liste0`; set uvcoverageNr = `awk '{if ($0 ~ /UV_COVERAGE /) {print $1;}}' liste0`
        set startimageNr = `awk '{if ($0 ~ /START_IMAGE /) {print $1;}}' liste0`; set prioriamgeNr = `awk '{if ($0 ~ /PRIOR_IMAGE /) {print $1;}}' liste0`


	  if( $model != "no" ) then
	    set modelNr = `awk '{if ($0 ~ /MODEL_IMAGE /) {print $1;}}' liste0`; set modelconvNr = `awk '{if ($0 ~ /MODEL_CONV /) {print $1;}}' liste0`
	    rm -f bestrec.fits;  fextract $fits\[$bestrecNr\] \!bestrec.fits\[0\]
	    rm -f bestrecconv.fits;  fextract $fits\[$bestrecconvNr\] \!bestrecconv.fits\[0\]
	    rm -f model.fits;  fextract $fits\[$modelNr\] \!model.fits\[0\]
	    rm -f modelconv.fits;  fextract $fits\[$modelconvNr\] \!modelconv.fits\[0\]
	  else
	    set khhm = -1
	    rm -f bestrec.fits;  fextract $fits\[$bestrecNr\] \!bestrec.fits\[0\]
            rm -f bestrecconv.fits;  fextract $fits\[$bestrecconvNr\] \!bestrecconv.fits\[0\]
	  endif

        if( X$khhm == "X" ) set khhm = 0.000
	echo "khh khhm = " $khh $khhm

	set khhm0 = $khh
	if( ($model != "no") && ($khhm != 0.000) ) set khh = $khhm

	rm -f liste khhm0.ll
        echo " calcVisf0 weightPower convscale = " $calcVisf0  $weightPower $convscale
	echo "-----------------------------------------------------------------"
        echo "# quality parameter of the best reconstruction:" >> liste
	echo "# ----------------------------------------------" >> liste 
        echo "# qrec   |  cost  |  chi2vis |  rresvis | chi2amp | rresamp | chi2phi  | rresphi | dist$convscale || FOV (mas) |Reg-Fct.|oradius step number|mu factor number|calcVisf0|weightPower| npix |startmode|startparam|directory  " >> liste
        echo  $qrec "   " $cost "  " $chi2vis "   " $rresvis "  " $chi2amp "   " $rresamp "   " $chi2phi "  " $rresphi "  " $khh "    " $fov "       " $regFunc "       " $oradiusStart $stepSize $oradiusNumber "     " $muStart $muFactor $muNumber "        " $calcVisf0  "       " $weightPower "      " $npix "     " $startmode "      " $startparam  "   " $cwd:t >> liste
      echo "Eingaben: "
      echo "fov fits convscale costFunc regFunc oradiusStart stepSize oradiusNumber muStart muFactor muNumber npix startmode startparam objname model tgauss tud tfdd weightPower calcVisf0: "
      echo $fov $fits $convscale $costFunc $regFunc $oradiusStart $stepSize $oradiusNumber $muStart $muFactor $muNumber $npix $startmode $startparam $objname $model
      echo $tgauss $tud $tfdd $weightPower $calcVisf0
      echo "-- IRBis.display.Mac.nt.ft.csh ---------------------------------------------------------------"

        echo "# quality parameter of the best reconstruction (dist from esorex.log):" >> khhm0.ll
        echo "# ----------------------------------------------" >> khhm0.ll
        echo "# qrec   |  cost  |  chi2vis |  rresvis | chi2amp | rresamp | chi2phi  | rresphi | dist$convscale || FOV (mas) |Reg-Fct.|oradius step number|mu factor number|calcVisf0|weightPower| npix |startmode|startparam|directory  " >> khhm0.ll
        echo  $qrec "   " $cost "  " $chi2vis "   " $rresvis "  " $chi2amp "   " $rresamp "   " $chi2phi "  " $rresphi "  " $khhm0 "    " $fov "       " $regFunc "       " $oradiusStart $stepSize $oradiusNumber "     " $muStart $muFactor $muNumber "        " $calcVisf0  "       " $weightPower "      " $npix "     " $startmode "      " $startparam  "   " $cwd:t >> khhm0.ll
   

# size0 : Kantenlaenge in Pixel des dargestellten Felder
        set size0 = 256

	# readout of the best reconstruction, unconvolved and convolved, the :
	foreach sc (sqrt lin log)
	  rm -f $fits:r.bestrec.$sc.jpeg $fits:r.bestrec.$sc.pdf $fits:r.bestrec.$sc.ps;  fextract $fits\[$bestrecNr\] \!bestrec.fits\[0\]; \
          set ff = bestrec.fits; ftcopy ''$ff'[-*,-*]' $ff.0; mv $ff.0 $ff; set disp = $sc; if( $sc == "lin" ) set disp = linear; $SCRIPTS/fits2ps.fits2bitmap.csh $ff $disp
	  mv $ff:r.$disp.ps $fits:r.bestrec.$sc.ps; mv $ff:r.$disp.pdf $fits:r.bestrec.$sc.pdf                                       # best rec unconvolved

          rm -f $fits:r.bestrecconv.$sc.jpeg $fits:r.bestrecconv.$sc.pdf $fits:r.bestrecconv.$sc.ps;  fextract $fits\[$bestrecconvNr\] \!bestrecconv.fits\[0\]; \
	  set ff = bestrecconv.fits; ftcopy ''$ff'[-*,-*]' $ff.0; mv $ff.0 $ff; set disp = $sc; if( $sc == "lin" ) set disp = linear; $SCRIPTS/fits2ps.fits2bitmap.csh $ff $disp
	  mv $ff:r.$disp.ps $fits:r.bestrecconv.$sc.ps; mv $ff:r.$disp.pdf $fits:r.bestrecconv.$sc.pdf                               # best rec convolved

	  rm -f $fits:r.uv.$sc.jpeg $fits:r.uv.$sc.pdf $fits:r.uv.$sc.ps;  fextract $fits\[$uvcoverageNr\] \!uv.fits\[0\]; \
	  set ff = uv.fits; ftcopy ''$ff'[-*,-*]' $ff.0; mv $ff.0 $ff; set disp = $sc; if( $sc == "lin" ) set disp = linear; $SCRIPTS/fits2ps.fits2bitmap.csh $ff $disp
	  mv $ff:r.$disp.ps $fits:r.uv.$sc.ps; mv $ff:r.$disp.pdf $fits:r.uv.$sc.pdf                     # uv coverage

	  if( $model == "no" ) then
	    rm -f $fits:r.start.$sc.jpeg $fits:r.start.$sc.pdf $fits:r.start.$sc.ps;  fextract $fits\[$startimageNr\] \!start.fits\[0\]; \
	    set ff = start.fits; ftcopy ''$ff'[-*,-*]' $ff.0; mv $ff.0 $ff; set disp = $sc; if( $sc == "lin" ) set disp = linear; $SCRIPTS/fits2ps.fits2bitmap.csh $ff $disp
            mv $ff:r.$disp.ps $fits:r.start.$sc.ps; mv $ff:r.$disp.pdf $fits:r.start.$sc.pdf                            # start image  unconvolved

	    rm -f $fits:r.prior.$sc.jpeg $fits:r.prior.$sc.pdf $fits:r.prior.$sc.ps;  fextract $fits\[$prioriamgeNr\] \!prior.fits\[0\]; \
	    set ff = prior.fits; ftcopy ''$ff'[-*,-*]' $ff.0; mv $ff.0 $ff; set disp = $sc; if( $sc == "lin" ) set disp = linear; $SCRIPTS/fits2ps.fits2bitmap.csh $ff $disp
            mv $ff:r.$disp.ps $fits:r.prior.$sc.ps; mv $ff:r.$disp.pdf $fits:r.prior.$sc.pdf                            # prior image  unconvolved
          endif

	  if( $model != "no" ) then
	    rm -f $fits:r.model.$sc.jpeg $fits:r.model.$sc.pdf $fits:r.model.$sc.ps;  fextract $fits\[$modelNr\] \!model.fits\[0\]; \
	    set ff = model.fits; ftcopy ''$ff'[-*,-*]' $ff.0; mv $ff.0 $ff; set disp = $sc; if( $sc == "lin" ) set disp = linear; $SCRIPTS/fits2ps.fits2bitmap.csh $ff $disp
            mv $ff:r.$disp.ps $fits:r.model.$sc.ps; mv $ff:r.$disp.pdf  $fits:r.model.$sc.pdf                         # model image  unconvolved

	    rm -f $fits:r.modelconv.$sc.jpeg $fits:r.modelconv.$sc.pdf $fits:r.modelconv.$sc.ps;  fextract $fits\[$modelconvNr\] \!modelconv.fits\[0\]; \
	    set ff = modelconv.fits; ftcopy ''$ff'[-*,-*]' $ff.0; mv $ff.0 $ff; set disp = $sc; if( $sc == "lin" ) set disp = linear; $SCRIPTS/fits2ps.fits2bitmap.csh $ff $disp
            mv $ff:r.$disp.ps $fits:r.modelconv.$sc.ps; mv $ff:r.$disp.pdf $fits:r.modelconv.$sc.pdf                  # model image convolved
	 endif
	end

          if( $model != "no" ) then
            set modelNr = `awk '{if ($0 ~ /MODEL_IMAGE /) {print $1;}}' liste0`; set modelconvNr = `awk '{if ($0 ~ /MODEL_CONV /) {print $1;}}' liste0`
            rm -f bestrec.fits;  fextract $fits\[$bestrecNr\] \!bestrec.fits\[0\]
            rm -f bestrecconv.fits;  fextract $fits\[$bestrecconvNr\] \!bestrecconv.fits\[0\]
            rm -f model.fits;  fextract $fits\[$modelNr\] \!model.fits\[0\]
            rm -f modelconv.fits;  fextract $fits\[$modelconvNr\] \!modelconv.fits\[0\]
          else
            rm -f bestrec.fits;  fextract $fits\[$bestrecNr\] \!bestrec.fits\[0\]
            rm -f bestrecconv.fits;  fextract $fits\[$bestrecconvNr\] \!bestrecconv.fits\[0\]
          endif

#	rm -f tt1.eps tt2.eps tt3.eps; cp gaussudfdda.ps tt1.eps; cp visa.ps tt2.eps; cp cpa.ps tt3.eps
	rm -f tt1.eps tt2.eps tt3.eps; cp gaussudfdda.ps tt1.eps; cp visa.ps tt2.eps; cp cpaa.ps tt3.eps

	foreach sc (sqrt lin log)
	if( $model == "no" ) then
          set i = 0
          foreach f ($fits:r.start.$sc.ps $fits:r.prior.$sc.ps $fits:r.bestrec.$sc.ps $fits:r.bestrecconv.$sc.ps $fits:r.uv.$sc.ps $fits:r.uv.$sc.ps)
	    @ i++; ps2eps $f; rm -f t$i.eps; mv $f:r.eps t$i.eps; ls -l t$i.eps
	  end
	else
          set i = 0
          foreach f ($fits:r.model.$sc.ps $fits:r.modelconv.$sc.ps $fits:r.bestrec.$sc.ps $fits:r.bestrecconv.$sc.ps $fits:r.uv.$sc.ps $fits:r.uv.$sc.ps)
            @ i++; ps2eps $f; rm -f t$i.eps; mv $f:r.eps t$i.eps; ls -l t$i.eps
          end
	endif

rm -f param.txt
echo "FOV $fov mas" >> param.txt
echo "npix $npix" >> param.txt
echo "cost $costFunc" >> param.txt
echo "reg $regFunc" >> param.txt
echo "weightpower $weightPower" >> param.txt
echo "orad $oradiusStart mas / $stepSize mas / $oradiusNumber" >> param.txt
echo "mu $muStart/$muFactor/$muNumber" >> param.txt
echo "startmode $startmode" >> param.txt
echo "startparam $startparam mas" >> param.txt
echo "superres $convscale" >> param.txt
echo "$objname" >> param.txt

rm -f results.txt
set chi2vis00 = `echo $chi2vis | awk '{printf "%12.3f", $1;}'`
set rresvis00 = `echo $rresvis | awk '{printf "%12.3f", $1;}'`
set chi2amp00 = `echo $chi2amp | awk '{printf "%12.3f", $1;}'`
set rresamp00 = `echo $rresamp | awk '{printf "%12.3f", $1;}'`
set chivisa00 = `echo $chivisa | awk '{printf "%12.3f", $1;}'`
set resratioa00 = `echo $resratioa | awk '{printf "%12.3f", $1;}'`
set chi2phi00   = `echo $chi2phi   | awk '{printf "%12.3f", $1;}'`
set rresphi00   = `echo $rresphi   | awk '{printf "%12.3f", $1;}'`
set chi2fta00   = `echo $chi2fta   | awk '{printf "%12.3f", $1;}'`
set ftresa00   = `echo $ftresa   | awk '{printf "%12.3f", $1;}'`
set khh00      = `echo $khh      | awk '{printf "%12.3f", $1;}'`
set qrec00     = `echo $qrec     | awk '{printf "%12.3f", $1;}'`
set phqrec00   = `echo $phqrec   | awk '{printf "%12.3f", $1;}'`
set khhm00     = `echo $khhm     | awk '{printf "%12.3f", $1;}'`
set qreca00    = `echo $qreca    | awk '{printf "%12.3f", $1;}'`
set phqreca00  = `echo $phqreca  | awk '{printf "%12.3f", $1;}'`
echo "Complex-Vis $chi2vis00 $rresvis00" >> results.txt
echo "V $chi2amp00 $rresamp00 $chivisa00 $resratioa00" >> results.txt
echo "PH $chi2phi00 $rresphi00 $chi2fta00 $ftresa00" >> results.txt
echo "dist qrec phqrec $khh00 $qrec00 $phqrec00  $khhm00 $qreca00 $phqreca00" >> results.txt


#    uv coverage plot:
     $SCRIPTS/uvplot.vis.csh; rm -f t66.eps; mv uv.ps t66.eps

        cp $SCRIPTS/plot6b.tex .; latex plot6b.tex; dvips plot6b.dvi; rm -f reks.$sc.ps; mv plot6b.ps reks.$sc.ps
	$SCRIPTS/plot6plus3text.csh; rm -f RR.konv.iterated.rek1.$sc.ps; mv plot6plus3text.ps RR.konv.iterated.rek1.$sc.ps
	rm -f plot6plus3text.* plot6b.* 
	end

        ls -l RR.konv.iterated.rek1.sqrt.ps RR.konv.iterated.rek1.lin.ps RR.konv.iterated.rek1.log.ps
        ls -l results.1.reconstruction.ps results.2.reconstruction.ps results.3.reconstruction.ps

# Erzeugung des ASCII-Files *.vis2 (altes Format fuer Fortran-Code von IRBis):
# Spaltenbelegung von *.vis2:
## set lam0 = `awk 'BEGIN{sum1=0.; sum2=0.;} { if($1!="#") {sum1=sum1+$2; sum2=sum2+1.;}; } END{print (sum1/sum2)*1000000.;}' $vis2`
#  ID     Vis^2   Vis^2Err   u[m]  v[m]


endif

