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

  set pdf  = $fits:r.$disp.pdf
  set png  = $fits:r.png
  set ps   = $fits:r.$disp.ps
  if( -f $pdf ) rm -f $pdf
  if( -f $ps  ) rm -f $ps
  if( -f $png ) rm -f $png

  set flag = `uname -a | awk '{ flag=0; if($0 ~ "Ubuntu") {flag=1;}; print flag; }'`
  if( $flag == 1 ) then
#   Ubuntu
    fits2bitmap --cmap $colormap --scale "$disp" $fits; convert $png $ps; convert $png $pdf; rm -f $png
  endif
  if( $flag == 0 ) then
#   Fedora or other
    fits2bitmap --cmap $colormap $fits; convert $png $ps; convert $png $pdf; rm -f $png
  endif
  ls -l $pdf $ps
endif
