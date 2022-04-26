This is a Fortran module translated from the Matlab code of Prof. Roberto Garrappa:

+ https://www.mathworks.com/matlabcentral/fileexchange/48154-the-mittag-leffler-function


--------------------------------

The package includes

+ **mod_mlf_garrappa.f90** is a portable Fortran module for calculation of Mittag Leffler function and its derivative with double precision (i.e. real(8) and complex(8)).

+ **tst_main_v02.f90** is a main program for more than 70 test caces which all are maintained in the folder **tcases**. I'm going to update more test cases.

+ **tst_main_v02.sh** is the bash (Linux) script that supports users to compile the package.

+ **tcases** includes more than 70 test caces for which I compute the Mittag-Leffler function and its derivative directly using the **FM** package of Prof. David M. Smith (link: http://dmsmith.lmu.build/)

Herein **v02** means **version 02**, i.e., the current version. Older stuffs are moved to **obsolete/**

If **tst_main_v02.sh** confuses, read comments inside **mod_mlf_garrappa.f90** and **tst_main_v02.f90** to know how to compile. 

Example with the gfortran compiler (4.9.2):

Step 1) Compiling:
```
  gfortran -O3 -c mod_mlf_garrappa.f90
  gfortran -O3 tst_main_v02.f90 mod_mlf_garrappa.o -o tst_main_v02.exe.gfc
```
Step 2) Running: To run testcase c01 with data reference file **tcases/tt_mlfm_c01.txt**,
```
 ./tst_main_v02.exe.gfc cas=c01 eep=6
 ./tst_main_v02.exe.gfc cas=c01 eep=8
 ./tst_main_v02.exe.gfc cas=c01 eep=10
 ./tst_main_v02.exe.gfc cas=c01
```
where cas=**c01** is to select the test case **c01**, and eep=**6** is to perform the computation with the tolerance **10^(-6)**. Details of the data files (all the test cases) can be found in **tcases/README.md**.

The above example was in Linux. It should work in Windows as well.

To calculate the generalized Mittag-Leffler function for three parameters _alpha_ (afa), _beta_ (bta) and _gamma_ (gma), and for each complex number z, perform the command line:
```
 ./tst_main_v02.exe.gfc afa=0.75  bta=1.0 gma=1.2  z='( -10.0, 1.0 )'
```
where _alpha = 0.75_, _beta = 1.0_, _gamma = 1.2_ and _z = -10 + i_. Make sure that there is no white-space around the character "=", as above. 

Then the outcome should be:
```
alpha = 7.5000000E-01, beta = 1.0000000E+00, gamma = 1.2000000E+00
  z  = (-1.00000000000000000E+01, 1.00000000000000000E+00),
E(z) = ( 8.72265651199321222E-03, 1.30548998870915031E-03),  time =  2.00000E-04 (s) (sec)
```
where _E(z)_ is the numerical value of the Mittag-Leffler function.

--------------------------------
Fortran syntax in use:
```
  e0d = genmlf ( afa, bta, gma, z0d )
  e1d = genmlf ( afa, bta, gma, z1d )
  e2d = genmlf ( afa, bta, gma, z2d )
  e3d = genmlf ( afa, bta, gma, z3d )
```
where the symbols 0d, 1d, 2d, and 3d stand for input when it given as scalar (0d), or arrays in 1d, 2d, 3d, respectively. Read **tst_main_v02.f90** for details. (There are lots of examples for input as scalar or array with various dimensions.)

--------------------------------

Contact: **viet204@mail.com**
