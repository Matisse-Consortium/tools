#!/bin/tcsh
#
#
# removing all *.ps *.eps and *.pdf, but not E.*.ps
set files = (E.*.ps bestrec.*.ps bestrecconv.*.ps)
mkdir container; mv $files container/
#
foreach f (*.ps)
#  if( ($f != $files[1]) && ($f != $files[2]) && ($f != $files[3]) ) then
    rm -f $f
#  endif
  rm *.pdf *.eps
end
mv container/* .; rmdir container/

