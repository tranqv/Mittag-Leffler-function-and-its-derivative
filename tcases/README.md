# Description of the test cases (roughly)

These files are output from direct computations of Mittag-Leffler function and its derivative
by adopting the high precision package FM (http://dmsmith.lmu.build/)

The sums of ML function and its derivative were accumulated directly with more than 1000 significant digits and with 
a relative error tolerance about 10^(-50).

This folder includes more than 70 test cases. For each the test cases, ML function and its derivative are computed on a line in complex
plane with a fixed angle, i.e. arg(z) = constant. 

So, what are inside the files? For example, let us open the file tt_mlfm_c01.txt, 

```
#   1001  0.850000000000000E+00   0.100000000000000E+01   0.000000000000000E+00
-4.500000000000000E+001  0.000000000000000E+000  3.693082757843946E-003  0.000000000000000E+000  8.486317294263145E-005  0.000000000000000E+000
....
```

What are there?

The 1st line: 
```
#   1001  0.850000000000000E+00   0.100000000000000E+01   0.000000000000000E+00
     N    alpha                   beta                    arg(z) (fixed)
```
the next line:
```
-4.500000000000000E+001  0.000000000000000E+000  3.693082757843946E-003  0.000000000000000E+000  8.486317294263145E-005  0.000000000000000E+000

  z_j    E(alpha,beta,z_j)   and  dE/dz -> 6 colums for 3 complex numbers.
```
and so on for j = 1,N, where N=1001.
