cf_expm
=======

`cf_expm` calculates the residuals and poles needed for the rational function (partial fraction) approximation to the matrix exponential.

It requires LAPACK and [FFTW](http://www.fftw.org/).

Usage
-------

Compile and link using `make`.

Run using `./cf_expm n` where `n` is the type (n, n) of approximation (i.e., number of terms):

    $ ./cf_expm 10

It can also be run without specifying n for the default (10):

    $ ./cf_expm

License
-------

`cf_expm` is released under the modified BSD license, see LICENSE for details.