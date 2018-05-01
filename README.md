This is a Fortran module based on the Matlab code of Prof. Robert Garrappa with a minor modification

+ the Matlab code: 
https://www.mathworks.com/matlabcentral/fileexchange/48154-the-mittag-leffler-function

I would like to thank him for allowing me to reuse his code so that this Fortran module can be made and to be released on public domain. Feel free to use or modify this module. In doing so, please cite references included inside the source codes.

--------------------------------

The package includes

+ mod_mlf_garrappa.f90 is a portable Fortran module for calculation of Mittag Leffler function and its derivative.

+ tst_main_vXX.f90 is the first main program for more than 70 test caces which all are maintained in the folder tcases. I'm going to update more testcases.

+ tst_main_vXX.sh is the bash (Linux) script that supports users to compile the package.

where _vXX marks the current version, e.g. _v01, _v02, etc. The newest version of the module is stayed the same name: mod_mlf_garrappa.f90, while older versions are marked by the _vXX.

In case the script tst_main_vXX.sh makes you to be confused, read inside the files mod_mlf_garrappa.f90 and tst_main_vXX.f90 to know how to compile. It is very simple. Here is an example with gfortran compiler (version 4.9.2)

Compiling:

+ gfortran -O3 -c mod_mlf_garrappa.f90

+ gfortran -O3 tst_main_vXX.f90 mod_mlf_garrappa.o -o tst_main_vXX.exe.gfc

Running: To run testcase c01 with data reference file tcases/tt_mlfm_c01.txt,

+ ./tst_main_vXX.exe.gfc cas=c01 eep=6
+ ./tst_main_vXX.exe.gfc cas=c01 eep=8
+ ./tst_main_vXX.exe.gfc cas=c01 eep=10
+ ./tst_main_vXX.exe.gfc cas=c01

where cas=c01 is to select the test case c01, and eep=6 is to perform with error estimate about 10^(-6). Details of the reference data files can be read at README.txt in the folder tcases (extract tcases.zip).

To calculate the generalized Mittag-Leffler function for three parameters alpha, beta and gamma, and for each complex number z, perform:

+ ./tst_main_v02.exe.gfc afa=0.75  bta=1.0 gma=1.2  z='( -10.0, 1.0 )'

which should give the output:

alpha = 7.5000000E-01, beta = 1.0000000E+00, gamma = 1.2000000E+00
  z  = (-1.00000000000000000E+01, 1.00000000000000000E+00),
E(z) = ( 8.72265651199321222E-03, 1.30548998870915031E-03),  time =  2.00000E-04 (s) (sec)

Make sure that there is no white-space around the character "=" as above.

--------------------------------
Syntax:

+  e0d = genmlf ( afa, bta, gma, z0d )
+  e1d = genmlf ( afa, bta, gma, z1d )
+  e2d = genmlf ( afa, bta, gma, z2d )
+  e3d = genmlf ( afa, bta, gma, z3d )

for *0d, *1d, *2d, and *3d stand for scalar (0d) and arrays in 1d, 2d, 3d, respectively.

--------------------------------

Please test this Fortran code by yourself. Let me know if you can find some regular case that it is FAILED. 

I hope that this package is helpful for those who work on this topic.

Every report on this code is welcomed at the emails tranquocviet@tdt.edu.vn or viet204@mail.com.

Cheers.
