# Mittag-Leffler-function-and-its-derivative

Package (version 0.01) includes

+ mod_mlf_garrappa.f90  is a portable Fortran module for calculation of Mittag Leffler function and its derivative.
+ tst_main_01.f90 is the first main program for more than 70 test caces (uploading ...)
+ tst_main_01.sh is the bash (Linux) script that supports users to compile the package. 

In case the script tst_main_01.sh makes you to be confused, read inside the files mod_mlf_garrappa.f90 and tst_main_01.f90 to know how to compile. It is very simple. Here is an example with gfortran compiler (version 4.9.2)

Compiling:

+ gfortran  -O3  -c  mod_mlf_garrappa.f90
+ gfortran  -O3  tst_main_01.f90  mod_mlf_garrappa.o  -o  tst_main_01.exe.gfc

Running: with test case 01 (ref. file: tcases/tt_mlfm_c01.txt), 

+ ./tst_main_01.exe.gfc cas=c01 eep=6
+ ./tst_main_01.exe.gfc cas=c01 eep=8
+ ./tst_main_01.exe.gfc cas=c01 eep=10
+ ./tst_main_01.exe.gfc cas=c01 

where cas=c01 is to select the test case c01, and eep=6 is to perform with error estimate about 10^(-6)

Cheers.
