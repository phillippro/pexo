This directory contains Chebyshev polynomial coefficients
for the JPL Planetary ephemeris DE430. DE430 includes 572
constants used as the basis for integration. This number exceeds
the original limit on the number of parameters in the program
testeph.f provided as an example for reading the ephemerides.

In addition, this version of DE430 includes Chebyshev polynomial
coefficients fit to the integrated value of TT-TDB evaluated
at the geocenter. Older versions of software will likely ignore
these coefficients.

A new version of testeph.f for reading this set of coefficients
will be released shortly.

DE430 documentation is available at
http://ipnpr.jpl.nasa.gov/progress_report/42-196/196C.pdf
