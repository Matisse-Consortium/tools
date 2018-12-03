#!/bin/tcsh
#
set ft = A.ft.dat
# ft:
# measured complex visibilities vs reconstructed complex visibilities
# nr lambda u v bl amp amperr phi phierr gd ud fdd rec1_amp rec1_phi rec2_amp rec2_phi rec3_amp rec3_phi
# 1  2      3 4 5  6   7      8   9      10 11 12  13       14       15       16       17       18


# set fac = `awk 'BEGIN{min=1e36; max=-1.;} { if(($1!="#")&&($5>0.)) {frad=sqrt($3*$3+$4*$4); fac=($5/frad)/$2; if(fac<min) {min=fac;}; if(fac>max) {max=fac;};}; } END{print fac;}' $vis2`
set fac = `awk 'BEGIN{min=1e36; max=-1.;} { if(($1!="#")&&($5>0.)) {frad=sqrt($3*$3+$4*$4); fac=($5/frad)/$2; if(fac<min) {min=fac;}; if(fac>max) {max=fac;};}; } END{print fac;}' $ft`

set lambdamean = `awk 'BEGIN{sum1=0.;sum2=0.;} { if(($1!="#")&&($5>0.)) {sum1=sum1+$2; sum2=sum2+1.;}; } END{print sum1/sum2;}' $ft`
set lambdamean0 = `echo $lambdamean | awk '{print $1*1000000.;}'`

rm -f uv.txt; awk '{ if($1!="#") {frad=sqrt($3*$3+$4*$4); fac=$5/frad; ux=$3*fac; uy=$4*fac; print ux,uy; print -ux,-uy;}; }' $ft > uv.txt
set maxbl = `awk 'BEGIN{max=-1.;} {rad=sqrt($1^2+$2^2); if(rad>max) {max=rad;};} END{print 1.2*max;}' uv.txt`
echo $maxbl
rm -f uv.ft.ps; $SCRIPTS/gnuplot.uvplott uv.txt uv.vis2.ps $maxbl; ls -l uv.ft.ps

rm -f uv.txt; awk -v lamref=$lambdamean -v fac=$fac '{ if($1!="#") { ux=$3*lamref*fac; uy=$4*lamref*fac; print ux,uy; print -ux,-uy;}; }' $ft > uv.txt
set maxbl = `awk 'BEGIN{max=-1.;} {rad=sqrt($1^2+$2^2); if(rad>max) {max=rad;};} END{print 1.2*max;}' uv.txt`
echo $maxbl
set text = `echo "lam0 = $lambdamean0 mu"`
echo $text
rm -f uv.ps; $SCRIPTS/gnuplot.spatialfreq uv.txt uv.ps $maxbl " $text "; ls -l uv.ps
