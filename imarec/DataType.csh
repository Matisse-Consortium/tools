#!/bin/tcsh
#
#
if($#argv == 0) then
  echo "Goal: estimation of the data types in the oifits file and therefore which image reconstructions are possilble."
  echo
  echo "Input: [ one oifits file]"
  echo
  echo
  echo "algo_mod: 1 = use bispectrum (OI_T3 & OI_VIS2), 2 = use complex visibilities (OI_VIS), "
  echo "          3 = use bispectrum and complex visibilities (OI_T3 & OI_VIS2 & OI_VIS)"
  echo ""
  echo "oit = 1 --> algo_mode = 1, 2, 3: "
  echo "oit = 2 --> algo_mode = 1"
  echo "oit = 3 --> algo_mode = 2"
  echo "oit = 4 --> algo_mode = 2"
  echo "oit = 0 --> no reconstruction possible"
else
  set oifits = $1
  rm -f $oifits:r.txt; fstruct $oifits > $oifits:r.txt
  rm -f $oifits:r.datatypes; awk 'BEGIN{z=0;} { if($2 == "BINTABLE") {if(($3=="OI_T3")||($3=="OI_VIS2")||($3=="OI_VIS")) {z=z+1; print $3;};}; }' $oifits:r.txt > $oifits:r.datatypes
  rm -f $oifits:r.datatypes.s; sort -k 1 $oifits:r.datatypes > $oifits:r.datatypes.s
  rm -f $oifits:r.datatypes.s.t
  awk 'BEGIN{z=0;} { if(z>0) {if(vor != $1) {print $1; vor=$1;}; }; if(z==0) {vor=$1; print $1; z=z+1;}; }' $oifits:r.datatypes.s > $oifits:r.datatypes.s.t
  set datatypes = `awk '{print $1;}' $oifits:r.datatypes.s.t`; rm -f $oifits:r.datatypes.s
  echo "==> the data types stored in the oifits file : "$datatypes
  echo
  set oit = `echo $#datatypes $datatypes | awk '{ if($1==3) {pos=1;}; if($1==2) {if( (($2=="OI_T3")&&($3=="OI_VIS2"))||(($2=="OI_VIS2")&&($3=="OI_T3")) ) {pos=2;}; if( (($2=="OI_T3")&&($3=="OI_VIS"))||(($2=="OI_VIS2")&&($3=="OI_VIS")) ) {pos=3;};}; if($1==1) {pos=0; if($2=="OI_VIS") {pos=4;};}; } END{print pos;}'`
  echo $oit
  echo
  if( $oit == 1 ) echo "==> possible algo_modes: 1, 2, 3"
  if( $oit == 2 ) echo "==> possible algo_modes: 1"
  if( $oit == 3 ) echo "==> possible algo_modes: 2"
  if( $oit == 4 ) echo "==> possible algo_modes: 2"
  if( $oit == 0 ) echo "==> possible algo_modes: none"

  echo "# oit = $oit"
  rm -f *.datatypes* *.txt
endif
