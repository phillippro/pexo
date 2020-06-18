# PEXO v2.0
[![DOI](https://zenodo.org/badge/210655784.svg)](https://zenodo.org/badge/latestdoi/210655784)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/phillippro/pexo/binder) <sup>[What is Binder?](https://mybinder.readthedocs.io/en/latest/)</sup>

Compared with previous models and packages, PEXO is general enough to account for binary motion and stellar reflex motions induced by planetary companions. PEXO is precise enough to treat various relativistic effects both in the Solar System and in the target system (Roemer, Shapiro, and Einstein delays).

PEXO is able to model timing to a precision of 1 ns, astrometry to a precision of 1 μas, and radial velocity to a precision of 1 μm/s. There are [pdf](https://github.com/phillippro/pexo/blob/master/docs/manual.pdf) and [html](http://rpubs.com/Fabo/pexo) versions of the manual available for instructions of how to use PEXO. The fitting mode and a python wrapper are in development and expected to be released soon.

The relevant paper was published by ApJS. If you use PEXO in your work, please cite the paper:
```
@article{Feng_2019,
	doi = {10.3847/1538-4365/ab40b6},
	url = {https://doi.org/10.3847%2F1538-4365%2Fab40b6},
	year = 2019,
	month = {oct},
	publisher = {American Astronomical Society},
	volume = {244},
	number = {2},
	pages = {39},
	author = {Fabo Feng and Maksym Lisogorskyi and Hugh R. A. Jones and Sergei M. Kopeikin and R. Paul Butler and Guillem Anglada-Escud{\'{e}} and Alan P. Boss},
	title = {{PEXO}: A Global Modeling Framework for Nanosecond Timing, Microarcsecond Astrometry, and $\upmu$m s-1 Radial Velocities},
	journal = {The Astrophysical Journal Supplement Series}
}
```

