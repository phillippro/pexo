Msun <- 1#solar mass
##Earth mass in unit of solar mass
Mearth <- MEMsol <- 3.003e-6
##Mercury mass
Mmercury <- 0.0553*MEMsol
##Venus mass
Mvenus <- 0.815*MEMsol
##Mars mass
Mmars <- 0.11*MEMsol
##Jupiter mass
Mjupiter <- 0.000954791915#solar mass
##Saturn mass
Msaturn <- 95.2*MEMsol
##Uranus mass
Muranus <- 14.6*MEMsol
##Neptune mass
Mneptune <- 17.2*MEMsol
##Moon
Mmoon <- Mearth/81.30056
##all masses of solar system star and planets (ssp)
Mssp <- c('Sun'=Msun,'Mercury'=Mmercury,'Venus'=Mvenus,'Earth'=Mearth,'Moon'=Mmoon,'Mars'=Mmars,'Jupiter'=Mjupiter,'Saturn'=Msaturn,'Uranus'=Muranus,'Neptune'=Mneptune)

###constants related to clock correction
IFTE.TEPH0 <- -65.564518e-6
IFTE.MJD0 <- 43144.0003725
####SOFA constants
Rsun.km <- 6.957e5#km
## Pi */
DPI <- 3.141592653589793238462643
## 2Pi */
D2PI <- 6.283185307179586476925287
## Radians to degrees */
rad2deg <- DR2D <- 57.29577951308232087679815
## Degrees to radians */
DD2R <- 1.745329251994329576923691e-2
## Radians to arcseconds */
DR2AS <- 206264.8062470963551564734
## Arcseconds to radians */
DAS2R <- 4.848136811095359935899141e-6
## Seconds of time to radians */
DS2R <- 7.272205216643039903848712e-5
## Arcseconds in a full circle */
TURNAS <- 1296000.0
## Milliarcseconds to radians */
pxConv <- DMAS2R <- DAS2R / 1e3
## Length of tropical year B1900 (days) */
DTY <- 365.242198781
## Seconds per day. */
DAYSEC <- 86400.0
## Days per Julian year */
DJY <- 365.25
## Seconds per Julian year
SJY <- DAYSEC*DJY
## Days per Julian century */
DJC <- 36525.0
## Days per Julian millennium */
DJM <- 365250.0
## Reference epoch (J2000.0), Julian Date */
DJ00 <- 2451545.0
## Julian Date of Modified Julian Date zero */
DJM0 <- 2400000.5
## Reference epoch (J2000.0), Modified Julian Date */
DJM00 <- 51544.5
## 1977 Jan 1.0 as MJD */
DJM77 <- 43144.0
## TT minus TAI (s) */
TTMTAI <- 32.184
## Astronomical unit (m, IAU 2012) */
DAU <- 149597870.7e3
au2km <- 149597870.7
## Speed of light (m/s) */
CMPS <- 299792458.0
## Speed of light (km/s) */
CKMPS <- CMPS*1e-3
## Light time for 1 au (s) */
AULT <- DAU/CMPS
## Number of Light seconds in one AU
AULTSC <- DAU/CMPS
## Speed of light (au per day) */
DC <- DAYSEC/AULT
## Speed of light (au per year) */
Cauyr <- YC <- DAYSEC/AULT*DJY
## Speed of light (au per s)
Caus <- SC <- CKMPS/au2km
## L_G = 1 − d(TT)/d(TCG) */
ELG <- 6.969290134e-10
## L_B = 1 − d(TDB)/d(TCB), and TDB (s) at TAI 1977/1/1.0 */
ELB <- 1.55051976772e-8
##tempo2 defination
IFTE.KM1 <- ELB/(1-ELB)
IFTE.K <- IFTE.KM1+1
IFTE.LC <- 1.48082686742e-8
TDB0 <- -6.55e-5
##https://www.iers.org/SharedDocs/Publikationen/EN/IERS/Publications/tn/TechnNote36/tn36_151.pdf?__blob=publicationFile&v=1
TT0 <- 2443144.5003725
## K
## Schwarzschild radius of the Sun (au) */
## = 2 * 1.32712440041e20 / (2.99792458e8)^2 / 1.49597870700e11 */
SRS <- 1.97412574336e-8#au
###Light time for crossing the Schwarzschild radius of the Sun (s)
SRSLT <- SRS*AULT
################Other constants
##au/yr to km/s
auyr2kms <- 4.740470463533347889
##Mass of Sun (kg)
Msun.kg <- 1.98892e30
##Gravitational constant
G <- 6.674e-11
##gauss constant
g <- 4*pi^2
##pc to au
pc2au <- 206264.8062470963551564734
pc2km <- pc2au*au2km
##pc to m
PCM <- 3.08568025e16
##pc to km
PCKM <- 3.08568025e13
##solar radius in m
Rsun.m <- 6.95700e8
##solar radius in au
Rsun.au <- 0.00465047
##escape velocity on the Sun's surface
Ve.sun <- 617.5#km/s
##G*Msun/Rsun; (km/s)^2
Phi.sun <- -G*Msun/Rsun.m/1000
###
#Tsun <- G*Msun/CMPS^3#second
Tsun <- 4.925490947e-6#second
Tsun.day <- Tsun/DAYSEC#day
Tsun.yr <- Tsun/DAYSEC/DJY#yr
SUNMASS <- Tsun#since c=g=1
