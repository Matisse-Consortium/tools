#
#
if($#argv == 0) then
  echo "Task: in an oifits file: scaling the complex amplitude with a factor, and scaling the V^2 values and T3 amplitude with factor^2 and factor^3, respectively."
  echo "  (getestet am 24-10-18  --> funktioniert!)"
  echo 
  echo "Input: [oifits file] [scaling factor for the complex visibility] [mit/ohne OI_VIS (1/0)]"
  echo 
  echo "Output: the oifits file *.sfactor.fits with the complex visibility, squared visibility, and triple amplitude scaled."
else
  set origfits = $1
  set scale    = $2
  set joivis   = $3

  set resfits  = $origfits:r.s$scale.fits
  rm -f $resfits

  set scale2   = `echo $scale | awk '{print $1*$1;}'`; echo $scale2
  set scale3   = `echo $scale | awk '{print $1*$1*$1;}'`; echo $scale3

# OI_VIS2 ----------- A ------------------
  rm -f out1.fits
  ftcalc <<EOF
  ${origfits}[OI_VIS2]
  out1.fits
  VIS2DATA
  VIS2DATA * $scale2
EOF

  rm -f out2.fits
  ftcalc <<EOF
  out1.fits[OI_VIS2]
  out2.fits
  VIS2ERR
  VIS2ERR * $scale2
EOF
  rm -f out1.fits; mv out2.fits out1.fits
# OI_VIS2 ----------- E ------------------

if( $joivis == 1 ) then
# OI_VIS ----------- A ------------------
  rm -f out2.fits
  ftcalc <<EOF
  out1.fits[OI_VIS]
  out2.fits
  VISAMP
  VISAMP * $scale
EOF
  rm -f out1.fits; mv out2.fits out1.fits

  rm -f out2.fits
  ftcalc <<EOF
  out1.fits[OI_VIS]
  out2.fits
  VISAMPERR
  VISAMPERR * $scale
EOF
  rm -f out1.fits; mv out2.fits out1.fits
# OI_VIS ----------- E ------------------
endif

# OI_T3  ----------- A ------------------
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
# OI_T3  ----------- E ------------------

  echo "---------------- Ende von ftcalc ------------------------"
  echo
  echo "-------------------------------"
  dir out1.fits
  mv out1.fits $resfits
  dir $resfits

endif

