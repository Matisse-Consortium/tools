#!/bin/tcsh
#
#
# removing all *.ps *.eps and *.pdf, but not E.*.ps
set files = (E.*.ps)
#
foreach f (*.ps)
  if( ($f != $files[1]) && ($f != $files[2]) && ($f != $files[3]) ) then
    rm -f $f
  endif
  rm *.pdf *.eps
end

