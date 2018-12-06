#!/bin/tcsh
#
set marke = $#argv
if($marke == 0) then
  echo "Goal: plot of measured and fitted data"
  echo
  echo "Input: [name of the ascii file of the data] [text for the measured data]"
  echo "       [column of the X axis in the file] [column of the measured data] [column of the errors of the measured data]"
  echo "       [text on X axis] [text on Y axis] [name of the ps plot]"
  echo "       [text for the fitted data] [column of the fitted data] ..."
else
  set file   = $argv[1]
  set textm  = "$argv[2]"
  set cx     = $argv[3]
  set cym    = $argv[4]
  set cem    = $argv[5]
  set textx  = "$argv[6]"
  set texty  = "$argv[7]"
  set psplot = $argv[8]
  if( $#argv == 10 ) then
    set fitteddata = ("$argv[9]" "$argv[10]")
  endif
  if( $#argv == 12 ) then
    set fitteddata = ("$argv[9]" "$argv[10]" "$argv[11]" "$argv[12]")
  endif
  if( $#argv == 14 ) then
    set fitteddata = ("$argv[9]" "$argv[10]" "$argv[11]" "$argv[12]" "$argv[13]" "$argv[14]")
  endif
  if( $#argv == 16 ) then
    set fitteddata = ("$argv[9]" "$argv[10]" "$argv[11]" "$argv[12]" "$argv[13]" "$argv[14]" "$argv[15]" "$argv[16]")
  endif
  set anzahl = `echo $#fitteddata | awk '{ z=$1/2; print z;}'`

#  echo "Input: " $file $textm $textf $cx $cym $cem $cyf $textx $texty $psplot

  rm -f $psplot

  rm -f skript0
  echo "# " >> skript0
  echo "# " >> skript0
  echo "gnuplot <<EOF" >> skript0
  #  set term postscript enh color solid \"Helvetica\" 20 >> skript0
  #  set term postscript portrait enh color solid \"Helvetica\" 18 >> skript0
  # echo set term postscript landscape enh color solid \"Helvetica\" 24 >> skript0
  echo set term postscript landscape enh color solid \"Helvetica\" 20 >> skript0
  #  set size 1.2,0.6 >> skript0
  #  set size 1.0,0.5   # portait >> skript0
  echo "set size 1.0,1.0    # landscape" >> skript0
  echo set output \"$psplot\" >> skript0
  #  set title "{/Helvetica=16 {$title}}" >> skript0
  #  set xla "{/Symbol l} [{/Symbol m}m]" >> skript0
  echo set xla \"$textx\" >> skript0
  echo set yla \"$texty\" >> skript0
  #  set yrange [vmin : vmax] >> skript0
  #  set xrange [xmin : xmax] >> skript0
  echo set mytics 10 >> skript0
#  echo set mxtics 10 >> skript0
  echo set mxtics 100 >> skript0
#  echo set format x \"%.1f\" >> skript0
  echo set format x \"%.0f\" >> skript0
  echo set format y \"%.1f\" >> skript0
  echo set bars small >> skript0
  echo set pointsize 0.5 >> skript0
  echo set grid >> skript0
  echo null\(x\) \= 0. >> skript0
  echo plot "\" >> skript0
  echo \"$file\" u $cx\:$cym\:$cem tit \"$textm\" wi errorbars lt 7 lw 2 pt 7, "\" >> skript0
  echo \"$file\" u $cx\:$cym notit wi lin lt 7, "\" >> skript0
  set i = 1
  while( $i <= $anzahl )
    set zweite = `echo $i | awk '{print $1*2;}'`
    set erste  = `echo $zweite | awk '{print $1-1;}'`
    set textf  = "$fitteddata[$erste]"
    set cyf    = $fitteddata[$zweite]
    if( $i == $anzahl ) then
      echo \"$file\" u $cx\:$cyf tit \"$textf\" wi lin lt $i >> skript0
    else
      echo \"$file\" u $cx\:$cyf tit \"$textf\" wi lin lt $i, "\" >> skript0
    endif
  @ i++
  end
  echo EOF >> skript0

  ls -l skript0
#  chmod +x skript0; ./skript0 rm -f skript0
  chmod +x skript0; ./skript0
 
  ls -l $psplot
endif
