cf_expm
=======

`cf_expm` calculates the residuals and poles needed for the rational function (partial fraction) approximation to the matrix exponential. It does this using the Carathéodory-Fejér method, and it is based on the MATLAB code in L.N. Trefethen, J.A.C. Weideman, T. Schmelzer, "Talbot quadratures and rational approximations," BIT Numer. Math. 46 (2006) 653–670. [doi:10.1007/s10543-006-0077-9](http://dx.doi.org/10.1007/s10543-006-0077-9)

It requires LAPACK and [FFTW](http://www.fftw.org/).

Usage
-------

Compile and link using `make`. Avoid using the compiler flag `-ffast-math`, as it seems to cause slightly incorrect results.

Run using `./cf_expm n` where `n` is the type (n, n) of approximation (i.e., number of terms):

    $ ./cf_expm 10

It can also be run without specifying n for the default (10):

    $ ./cf_expm

License
-------

`cf_expm` is released under the modified BSD license, see LICENSE for details.

Author
------

Created by [Kyle Niemeyer](http://kyleniemeyer.com). Email address: [kyle.niemeyer@gmail.com](mailto:kyle.niemeyer@gmail.com)