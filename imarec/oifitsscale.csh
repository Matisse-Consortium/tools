#
#
if($#argv == 0) then
  echo "Tasks: scaling of interferometric data in an oifits file: "
  echo "       a) scale factor >0: scaling the complex amplitudes V, the V^2 values and T3 amplitudes (and their errors) with factor, factor^2 and factor^3, respectively."
  echo "       b) scale factor <0: scaling only the errors of the complex amplitudes, V^2 values and T3 values with |factor|." 
  echo "  (getestet am 24-10-18  --> funktioniert!) (b) hinzugefuegt am 12-3-21 --> funktioniert!)i (a) getestet nur fuer V^2 --> funktioniert!)"
  echo 
  echo "Input: [oifits file] [scaling factor: >0 for the complex visibility; <0 for errors only] [mit/ohne OI_VIS (1/0)]"
  echo 
  echo "Output: the oifits file *.s.#scale.fits with the complex visibility, squared visibility, and triple amplitude scaled."
else
  set origfits = $1
  set scale0   = $2
  set joivis   = $3

  set scale = `echo $scale0 | awk '{print sqrt($1*$1);}'`; echo $scale
  set marke = `echo $scale0 | awk '{m=1; if($1<0) {m=0;}; print m;}'`; echo $marke   # marke = 1 ==> a) scale factor >0: scaling the complex amplitudes V, the V^2 values and T3 amplitudes (and their errors) with factor, factor^2 and factor^3, respectively.
  										     # marke = 0 ==> b) scale factor <0: scaling only the errors of the complex amplitudes, V^2 values and T3 values with |factor|.

  set resfits  = $origfits:r.s.$scale0.fits
  rm -f $resfits


  set scale2   = `echo $scale | awk '{print $1*$1;}'`; echo $scale2
  set scale3   = `echo $scale | awk '{print $1*$1*$1;}'`; echo $scale3

  rm -f out1.fits; cp $origfits out1.fits

# OI_VIS2 ----------- A ------------------
  if($marke == 1) then
  rm -f out2.fits
  ftcalc <<EOF
  out1.fits[OI_VIS2]
  out2.fits
  VIS2DATA
  VIS2DATA * $scale2
EOF
  rm -f out1.fits; mv out2.fits out1.fits

  rm -f out2.fits
  ftcalc <<EOF
  out1.fits[OI_VIS2]
  out2.fits
  VIS2ERR
  VIS2ERR * $scale2
EOF
  rm -f out1.fits; mv out2.fits out1.fits
  endif

  if($marke == 0) then
  rm -f out2.fits
  ftcalc <<EOF
  out1.fits[OI_VIS2]
  out2.fits
  VIS2ERR
  VIS2ERR * $scale
EOF
  rm -f out1.fits; mv out2.fits out1.fits
  endif
# OI_VIS2 ----------- E ------------------

if( $joivis == 1 ) then
# OI_VIS ----------- A ------------------
  if($marke == 1) then
  rm -f out2.fits
  ftcalc <<EOF
  out1.fits[OI_VIS]
  out2.fits
  VISAMP
  VISAMP * $scale
EOF
  rm -f out1.fits; mv out2.fits out1.fits
  endif

  rm -f out2.fits
  ftcalc <<EOF
  out1.fits[OI_VIS]
  out2.fits
  VISAMPERR
  VISAMPERR * $scale
EOF
  rm -f out1.fits; mv out2.fits out1.fits

  if($marke == 0) then
  rm -f out2.fits
  ftcalc <<EOF
  out1.fits[OI_VIS]
  out2.fits
  VISPHIERR
  VISPHIERR * $scale
EOF
  rm -f out1.fits; mv out2.fits out1.fits
  endif
# OI_VIS ----------- E ------------------
endif

# OI_T3  ----------- A ------------------
  if($marke == 1) then
  rm -f out2.fits
  ftcalc <<EOF
  out1.fits[OI_T3]
  out2.fits
  T3AMP
  T3AMP * $scale3
EOF
  rm -f out1.fits; mv out2.fits out1.fits

  rm -f out2.fits
  ftcalc <<EOF
  out1.fits[OI_T3]
  out2.fits
  T3AMPERR
  T3AMPERR * $scale3
EOF
  rm -f out1.fits; mv out2.fits out1.fits
  endif

  if($marke == 0) then
  rm -f out2.fits
  ftcalc <<EOF
  out1.fits[OI_T3]
  out2.fits
  T3PHIERR
  T3PHIERR * $scale
EOF
  rm -f out1.fits; mv out2.fits out1.fits

  rm -f out2.fits
  ftcalc <<EOF
  out1.fits[OI_T3]
  out2.fits
  T3AMPERR
  T3AMPERR * $scale
EOF
  rm -f out1.fits; mv out2.fits out1.fits
  endif
# OI_T3  ----------- E ------------------

  echo "---------------- Ende von ftcalc ------------------------"
  echo
  echo "-------------------------------"
  dir out1.fits
  mv out1.fits $resfits
  dir $resfits

endif

