#!/bin/tcsh
#
#
if($#argv == 0) then
  echo "Goal: estimation of the data stored in the oifits file and therefore which image reconstructions are possilble."
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
  rm -f $oifits:r.datatypes ; awk 'BEGIN{z=0;} { if($2 == "BINTABLE") {if(($3=="OI_T3")||($3=="OI_VIS2")||($3=="OI_VIS")) {z=z+1; print $3,$1;};}; }' $oifits:r.txt > $oifits:r.datatypes
  rm -f $oifits:r.datatypes.s; sort -k 1 $oifits:r.datatypes > $oifits:r.datatypes.s
  set index_OI_VIS  = `awk '{ if($1 == "OI_VIS")  {print $2;}; }' $oifits:r.datatypes.s`; echo "-- index_OI_VIS = $index_OI_VIS --"
  set index_OI_VIS2 = `awk '{ if($1 == "OI_VIS2") {print $2;}; }' $oifits:r.datatypes.s`; echo "-- index_OI_VIS2 = $index_OI_VIS2 --"
  set index_OI_T3   = `awk '{ if($1 == "OI_T3")   {print $2;}; }' $oifits:r.datatypes.s`; echo "-- index_OI_T3 = $index_OI_T3 --"
  rm -f $oifits:r.datatypes $oifits:r.txt

  set OIVISja = 1
  if( $#index_OI_VIS > 0 ) then
   rm -f tt.txt0
   set zz = 0
   foreach iinn ($index_OI_VIS)
	   rm -f tt.txt00
	   fdump $oifits+$iinn tt.txt00 FLAG -
	   @ zz++
	   if($zz == 1) then
		   cp tt.txt00 tt.txt0
           else
		   rm -f tt00; cat tt.txt0 tt.txt00 > tt00; mv tt00 tt.txt0
           endif
   end
   rm -f tt.txt0.1; grep F tt.txt0 > tt.txt0.1; set ff = `awk '{ if(($1 == "F")||($2 == "F")) {print $0;}; }' tt.txt0.1`
   # echo $ff; echo $#ff
   if( $#ff == 0 ) then
    set OIVISja = 0
    echo "OI_VIS is in $oifits, but all data are flagged."
   endif
  else
   set OIVISja = 0
   echo "OI_VIS is not in $oifits."
  endif
  if($OIVISja == 0) then
   rm -f tt0; awk '{ if($1 != "OI_VIS") {print $0;}; }' $oifits:r.datatypes.s > tt0; mv tt0 $oifits:r.datatypes.s
  endif
  rm -f tt.txt0 tt.txt0.1

  set OIVIS2ja = 1
  if( $#index_OI_VIS2 > 0 ) then
   rm -f tt.txt0
   set zz = 0
   foreach iinn ($index_OI_VIS2)
	   rm -f tt.txt00
	   fdump $oifits+$iinn tt.txt00 FLAG -
	   @ zz++
	   if($zz == 1) then
		   cp tt.txt00 tt.txt0
           else
		   rm -f tt00; cat tt.txt0 tt.txt00 > tt00; mv tt00 tt.txt0
           endif
   end
   rm -f tt.txt0.1; grep F tt.txt0 > tt.txt0.1; set ff = `awk '{ if(($1 == "F")||($2 == "F")) {print $0;}; }' tt.txt0.1`
   # echo $ff; echo $#ff
   if( $#ff == 0 ) then
    set OIVIS2ja = 0
    echo "OI_VIS2 is in $oifits, but all data are flagged."
   endif
  else
   set OIVIS2ja = 0
   echo "OI_VIS2 is not in $oifits."
  endif
  if($OIVIS2ja == 0) then
   rm -f tt0; awk '{ if($1 != "OI_VIS2") {print $0;}; }' $oifits:r.datatypes.s > tt0; mv tt0 $oifits:r.datatypes.s
  endif
  rm -f tt.txt0 tt.txt0.1

  set OIT3ja = 1
  if( $#index_OI_T3 > 0 ) then
   rm -f tt.txt0
   set zz = 0
   foreach iinn ($index_OI_T3)
	   rm -f tt.txt00
	   fdump $oifits+$iinn tt.txt00 FLAG -
	   @ zz++
	   if($zz == 1) then
		   cp tt.txt00 tt.txt0
           else
		   rm -f tt00; cat tt.txt0 tt.txt00 > tt00; mv tt00 tt.txt0
           endif
   end
   rm -f tt.txt0.1; grep F tt.txt0 > tt.txt0.1; set ff = `awk '{ if(($1 == "F")||($2 == "F")) {print $0;}; }' tt.txt0.1`
   # echo $ff; echo $#ff
   if( $#ff == 0 ) then
    set OIT3ja = 0
    echo "OI_T3 is in $oifits, but all data are flagged."
   endif
  else
   set OIT3ja = 0
   echo "OI_T3 is not in $oifits."
  endif
  if($OIT3ja == 0) then
   rm -f tt0; awk '{ if($1 != "OI_T3") {print $0;}; }' $oifits:r.datatypes.s > tt0; mv tt0 $oifits:r.datatypes.s
  endif
  rm -f tt.txt0 tt.txt0.1

  rm -f $oifits:r.datatypes.s.t
  awk 'BEGIN{z=0;} { if(z>0) {if(vor != $1) {print $1; vor=$1;}; }; if(z==0) {vor=$1; print $1; z=z+1;}; }' $oifits:r.datatypes.s > $oifits:r.datatypes.s.t
  set datatypes = `awk '{print $1;}' $oifits:r.datatypes.s.t`; rm -f $oifits:r.datatypes.s
  echo "==> the data stored in the oifits file : "$datatypes
  echo
  set oit = `echo $#datatypes $datatypes | awk '{ if($1==3) {pos=1;}; if($1==2) {if( (($2=="OI_T3")&&($3=="OI_VIS2"))||(($2=="OI_VIS2")&&($3=="OI_T3")) ) {pos=2;}; if( (($2=="OI_T3")&&($3=="OI_VIS"))||(($2=="OI_VIS2")&&($3=="OI_VIS")) ) {pos=3;};}; if($1==1) {pos=0; if($2=="OI_VIS") {pos=4;};}; } END{print pos;}'`
  echo $oit
  echo
  rm -f $oifits:r.datatypes.s.t
  if( $oit == 1 ) echo "==> possible algo_modes: 1, 2, 3"
  if( $oit == 2 ) echo "==> possible algo_modes: 1"
  if( $oit == 3 ) echo "==> possible algo_modes: 2"
  if( $oit == 4 ) echo "==> possible algo_modes: 2"
  if( $oit == 0 ) echo "==> possible algo_modes: none"

  echo "# oit = $oit"
#  rm -f *.datatypes* *.txt
endif
