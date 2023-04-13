#!/bin/bash
#

SRCDIR="../../src"

for ffile in `ls ${SRCDIR}/*/*.[fF]90` ; do
 modname=$(basename $(basename $ffile .f90) .F90)
 if [ $# -eq 1 ];  then
   if [ $1 = $modname ]; then
     echo "Proposing signature file for: $modname"
     f2py --overwrite-signature -m ${modname} -h ${modname}.pyf-proposed $ffile
   fi
 else
   echo "Proposing signature file for: $modname"
   f2py --overwrite-signature -m ${modname} -h ${modname}.pyf-proposed $ffile
 fi
done
