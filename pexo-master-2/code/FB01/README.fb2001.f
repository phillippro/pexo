fb2001.f has been written by Alan W. Irwin <irwin@beluga.phys.uvic.ca> to
read fb2001.dat (the definitive mass-corrected version of the coefficients
for the Fairhead and Bretagnon analytical series for the time ephemeris) and
from that series calculate the time ephemeris as a function of teph.  The
routine is compiled and run as follows:

gfortran -O2 fb2001.f -o fb2001
./fb2001 <fb2001.in >fb2001.out

fb2001.out contains the independent variable (teph expressed as a Julian
date), the time-ephemeris (without linear term) in days, and its derivative
(dimensionless).

The independent variable of this routine, teph is related to TCB via
equation 1 of Irwin and Fukushima, 1999 A&A 348, 642. That reference also
presents the definitive time ephemeris, TE405, whose maximum errors are
estimated to be only 0.1 ns from 1600 to 2200 and which is available at
ftp://astroftp.phys.uvic.ca/pub/irwin/tephemeris/ as an extremely fast
routine to evaluate TE405 from its Chebyshev form.

It turns out the present series whose coefficients are given by fb2001.dat
produces results that are in close numerical agreement with FB3C (a
non-definitive version of the FB series that had been corrected by Irwin and
Fukushima for mass errors).  That is, if you calculate a time ephemeris
using fb2001 and compare with TE405, you get a plot that is almost
indistinguishable from the last "FB3C" panel in Figure 3 of Irwin and
Fukushima with maximum errors of 15 ns from 1600 to 2200 and substantially
less than that for smaller epoch ranges near 2000.

Alan W. Irwin 2003-11-13
