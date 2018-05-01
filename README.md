This is a Fortran module based on the Matlab code of Prof. Robert Garrappa with a minor modification

+ the Matlab code: 
https://www.mathworks.com/matlabcentral/fileexchange/48154-the-mittag-leffler-function

I would like to thank him for allowing me to reuse his code so that this Fortran module can be made and to be released on public domain. Feel free to use or modify this module. In doing so, please cite references included inside the source codes.

--------------------------------

The package includes

+ mod_mlf_garrappa.f90 is a portable Fortran module for calculation of Mittag Leffler function and its derivative.

+ tst_main_vXX.f90 is the first main program for more than 70 test caces which all are maintained in the folder tcases. I'm going to update more testcases.

+ tst_main_vXX.sh is the bash (Linux) script that supports users to compile the package.

where _vXX marks the current version, e.g. _v01, _v02, etc. The newest version of the module is stayed the same name: mod_mlf_garrappa.f90, while older version are marked by the _vXX.

In case the script tst_main_vXX.sh makes you to be confused, read inside the files mod_mlf_garrappa.f90 and tst_main_vXX.f90 to know how to compile. It is very simple. Here is an example with gfortran compiler (version 4.9.2)

Compiling:

+ gfortran -O3 -c mod_mlf_garrappa.f90

+ gfortran -O3 tst_main_vXX.f90 mod_mlf_garrappa.o -o tst_main_vXX.exe.gfc

Running: with test case 01 (ref. file data: tcases/tt_mlfm_c01.txt),

+ ./tst_main_vXX.exe.gfc cas=c01 eep=6
+ ./tst_main_vXX.exe.gfc cas=c01 eep=8
+ ./tst_main_vXX.exe.gfc cas=c01 eep=10
+ ./tst_main_vXX.exe.gfc cas=c01

where cas=c01 is to select the test case c01, and eep=6 is to perform with error estimate about 10^(-6). 

Details of the reference file data can be read at the README file inside the folder tcase.

--------------------------------

Please test this Fortran code by yourself. Let me know if you can find some regular case that it is FAILED. 

I hope that this package is helpful for those who work on this topic.

Every report on this code is welcomed at the emails tranquocviet@tdt.edu.vn or viet204@mail.com.

Cheers.
