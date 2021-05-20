#!/bin/tcsh
#
# set echo
set psfiles = ()
set sqpsfiles = ()
set logpsfiles = ()

if( $#argv == 0 ) then
#  the reconstructions are displayed with increasing nr in E.nr 
 set i = 1
 while( $i > 0 )
  set dir0 = E.$i
  if( -d $dir0 ) then
    set psfiles = ($psfiles $dir0/E.lin.ps)
    set sqpsfiles = ($sqpsfiles $dir0/E.sqrt.ps)
    set logpsfiles = ($logpsfiles $dir0/E.log.ps)
  else
    goto Label111
  endif
  @ i++
 end
else
# the reconstructions as sorted in E.liste.00.sorted (according to increasing qrec) or E.liste.03.sorted (according to increasing cpqrec)
 set Eliste = $1
 set i = 1
 while( $i > 0 )
  set dir0 = `awk -v i=$i 'BEGIN{z=0;} { z=z+1; if(z==i) {print $23;}; }' $Eliste`
  if( X$dir0 == "X" ) set dir0 = E.0
  if( -d $dir0 ) then
    set psfiles = ($psfiles $dir0/E.lin.ps)
    set sqpsfiles = ($sqpsfiles $dir0/E.sqrt.ps)
    set logpsfiles = ($logpsfiles $dir0/E.log.ps)
  else
    goto Label111
  endif
 @ i++
 end
endif

Label111:
echo $psfiles
echo "$#psfiles Verzeichnisse sind vorhanden."


set plot = $cwd:t.lin.ps
#rm -f $plot; psmerge -o$plot $psfiles
rm -f $plot; psjoin $psfiles > $plot
set plotsq = $cwd:t.sqrt.ps
#rm -f $plotsq; psmerge -o$plotsq $sqpsfiles
rm -f $plotsq; psjoin $sqpsfiles > $plotsq
set plotlog = $cwd:t.log.ps
#rm -f $plotlog; psmerge -o$plotlog $logpsfiles
rm -f $plotlog; psjoin $logpsfiles > $plotlog
dir $plot $plotsq $plotlog
