#!/bin/tcsh
#
set marke = $#argv
if($marke == 0) then
  echo "Goal: convert a fits image into pdf image with heat lookuptable and linear/sqrt/log display."
  echo
  echo "input: [fits image] [linear/sqrt/log]"
else
  set fits = $1
  set disp = $2

##  set colormap = "hot"
#  set colormap = "afmhot"
  set colormap = "gist_heat"
##  set colormap = "autumn"
#  set colormap = "afmhot"

	echo "------------------ fits2ps.fits2bitmap.csh --------- 1 -----------"
  set pdf  = $fits:r.$disp.pdf
  set png  = $fits:r.png
  set ps   = $fits:r.$disp.ps
  if( -f $pdf ) rm -f $pdf
  if( -f $ps  ) rm -f $ps
  if( -f $png ) rm -f $png

	echo "------------------ fits2ps.fits2bitmap.csh --------- 2 -----------"
  set flag = `uname -a | awk '{ flag=0; if($0 ~ "Ubuntu") {flag=1;}; print flag; }'`
#  if( $flag == 1 ) then
   if( 2 == 1 ) then
#   Ubuntu
    fits2bitmap --cmap $colormap --scale "$disp" $fits
	echo "------------------ fits2ps.fits2bitmap.csh --------- 3 -----------"
    if($status != 0) then
      fits2bitmap --cmap $colormap --stretch "$disp" $fits
    endif
    convert $png $ps; convert $png $pdf; rm -f $png
  endif
	echo "------------------ fits2ps.fits2bitmap.csh --------- 4 -----------"
#  if( $flag == 0 ) then
   if( 1 == 1 ) then
#   Fedora or other
    rm -f {$fits:r}mr180.fits; imrot -l -r 180 $fits
    rm -f $fits.ORIG; mv $fits $fits.ORIG; mv {$fits:r}mr180.fits $fits
    fits2bitmap --cmap $colormap --stretch "$disp" $fits
    if($status != 0) then
      rm -f $fits; mv $fits.ORIG $fits
      fits2bitmap --cmap $colormap --scale "$disp" $fits
    endif
    if(-f $fits.ORIG) then
      rm -f $fits; mv $fits.ORIG $fits
    endif
    convert $png $ps; convert $png $pdf; rm -f $png
	echo "------------------ fits2ps.fits2bitmap.csh --------- 5 -----------"
  endif
  ls -l $pdf $ps
endif
