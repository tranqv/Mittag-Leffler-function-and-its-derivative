#!/bin/bash 

   mod=mod_mlf_garrappa
   pro=tst_main_v02

   com=gfortran
#  com=ifort

   iffl="-O3 -ipo -fno-alias -static -zero"
#  iffl="-O3"
   ifex="ifc" 

   gffl="-O3"
   gfex="gfc" 



   if [ $# -gt 0 ]
   then 
      com=$1
      if [ "$com" != "ifort" ] && [ "$com" != "gfortran" ] 
      then 
         echo "Support compliers: ifort, gfortran. Default: gfortran"
         echo "Syntax:    $0 ifort"
         echo "or,        $0 gfortran"
         echo "or,        $0 "
         echo "or,        $0 with_some_silly_name"
         com="gfortran"
      fi 
    fi 

   if [ $com = "ifort" ]
   then 
      flag=$iffl
      exte=".$ifex"
   else
      flag=$gffl
      exte=".$gfex"
   fi 

   echo "$com ${flag} -c ${mod}.f90"
   echo "$com ${flag} ${pro}.f90 ${mod}.o -o ${pro}.exe${exte}"

   rm -rf ${pro}.exe *.mod *.o 

   $com ${flag} -c ${mod}.f90                            \
&& $com ${flag} ${pro}.f90 ${mod}.o -o ${pro}.exe${exte} \
&& (echo "cp  ${pro}.exe${exte}  ../" ; cp ${pro}.exe${exte} ../ ) \
&& cat << OUT
      module : $mod 
      compler: $com 

      Example: to run testcase c40 with data reference file tt_mlfm_c40.txt

      ./${pro}.exe${exte} cas=c40 eep=6 
      ./${pro}.exe${exte} cas=c40 eep=8 
      ./${pro}.exe${exte} cas=c40 eep=10
      ./${pro}.exe${exte} cas=c40 eep=15
      ./${pro}.exe${exte} cas=c40
OUT
