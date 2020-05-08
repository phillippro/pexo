sofa_Tttcg <- function(tt){
####################################
##  Time scale transformation:  Terrestrial Time, TT, to Geocentric
##  Coordinate Time, TCG.
##
##  This function is part of the International Astronomical Union's
##  SOFA (Standards of Fundamental Astronomy) software collection.
##
##  Status:  canonical.
##
##  Given:
##     tt1,tt2    double    TT as a 2-part Julian Date
##
##  Returned:
##     tcg1,tcg2  double    TCG as a 2-part Julian Date
##
##  Returned (function value):
##                int       status:  0 = OK
##
##  Note:
##
##     tt1+tt2 is Julian Date, apportioned in any convenient way between
##     the two arguments, for example where tt1 is the Julian Day Number
##     and tt2 is the fraction of a day.  The returned tcg1,tcg2 follow
##     suit.
##
##  References:
##
##     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
##     IERS Technical Note No. 32, BKG (2004)
##
##     IAU 2000 Resolution B1.9
##
##  This revision:  2019 June 20
##
##  SOFA release 2019-07-22
##
##  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
####################################
    t77t <- DJM77 + TTMTAI/DAYSEC

    ## TT to TCG rate
    elgg <- ELG/(1.0-ELG)

    tt1 <- tt[,1]
    tt2 <- tt[,2]
    ##Result, safeguarding precision.
    tcg1 <- tt1
    tcg2 <- tt2 + ( ( tt1 - DJM0 ) + ( tt2 - t77t ) ) * elgg
    cbind(tcg1,tcg2)
}

sofa_Tcbtdb <- function(tcb){
####################################
##  Time scale transformation:  Barycentric Coordinate Time, TCB, to
##  Barycentric Dynamical Time, TDB.
##
##  This function is part of the International Astronomical Union's
##  SOFA (Standards of Fundamental Astronomy) software collection.
##
##  Status:  canonical.
##
##  Given:
##     tcb1,tcb2  double    TCB as a 2-part Julian Date
##
##  Returned:
##     tdb1,tdb2  double    TDB as a 2-part Julian Date
##
##  Returned (function value):
##                int       status:  0 = OK
##
##  Notes:
##
##  1) tcb1+tcb2 is Julian Date, apportioned in any convenient way
##     between the two arguments, for example where tcb1 is the Julian
##     Day Number and tcb2 is the fraction of a day.  The returned
##     tdb1,tdb2 follow suit.
##
##  2) The 2006 IAU General Assembly introduced a conventional linear
##     transformation between TDB and TCB.  This transformation
##     compensates for the drift between TCB and terrestrial time TT,
##     and keeps TDB approximately centered on TT.  Because the
##     relationship between TT and TCB depends on the adopted solar
##     system ephemeris, the degree of alignment between TDB and TT over
##     long intervals will vary according to which ephemeris is used.
##     Former definitions of TDB attempted to avoid this problem by
##     stipulating that TDB and TT should differ only by periodic
##     effects.  This is a good description of the nature of the
##     relationship but eluded precise mathematical formulation.  The
##     conventional linear relationship adopted in 2006 sidestepped
##     these difficulties whilst delivering a TDB that in practice was
##     consistent with values before that date.
##
##  3) TDB is essentially the same as Teph, the time argument for the
##     JPL solar system ephemerides.
##
##  Reference:
##
##     IAU 2006 Resolution B3
##
##  This revision:  2019 June 20
##
##  SOFA release 2019-07-22
##
##  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
####################################
## 1977 Jan 1 00:00:32.184 TT, as two-part JD
   t77td <- DJM0 + DJM77
   t77tf <- TTMTAI/DAYSEC

## TDB (days) at TAI 1977 Jan 1.0
   tdb0 <- TDB0/DAYSEC
   tcb1 <- tcb[,1]
   tcb2 <- tcb[,2]
##Result, safeguarding precision.
   d <- tcb1 - t77td
   tdb1 <- tcb1
   tdb2 <- tcb2 + tdb0 - ( d + ( tcb2 - t77tf ) ) * ELB
   cbind(tdb1,tdb2)
}

sofa_Tttdb <- function(tt){
####################################
##  Time scale transformation:  Terrestrial Time, TT, to Barycentric
##  Dynamical Time, TDB.
##
##  This function is part of the International Astronomical Union's
##  SOFA (Standards of Fundamental Astronomy) software collection.
##
##  Status:  canonical.
##
##  Given:
##     tt1,tt2    double    TT as a 2-part Julian Date
##     dtr        double    TDB-TT in seconds
##
##  Returned:
##     tdb1,tdb2  double    TDB as a 2-part Julian Date
##
##  Returned (function value):
##                int       status:  0 = OK
##
##  Notes:
##
##  1) tt1+tt2 is Julian Date, apportioned in any convenient way between
##     the two arguments, for example where tt1 is the Julian Day Number
##     and tt2 is the fraction of a day.  The returned tdb1,tdb2 follow
##     suit.
##
##  2) The argument dtr represents the quasi-periodic component of the
##     GR transformation between TT and TCB.  It is dependent upon the
##     adopted solar-system ephemeris, and can be obtained by numerical
##     integration, by interrogating a precomputed time ephemeris or by
##     evaluating a model such as that implemented in the SOFA function
##     iauDtdb.   The quantity is dominated by an annual term of 1.7 ms
##     amplitude.
##
##  3) TDB is essentially the same as Teph, the time argument for the JPL
##     solar system ephemerides.
##
##  References:
##
##     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
##     IERS Technical Note No. 32, BKG (2004)
##
##     IAU 2006 Resolution 3
##
##  This revision:  2019 June 20
##
##  SOFA release 2019-07-22
##
##  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
####################################
##Result, safeguarding precision.
    tt1 <- tt[1]
    tt2 <- tt[2]
###only an approximation; absolute value should be calculated according to JPL ephemeris
    g <- 6.24 + 0.017202*(sum(tt)-2451545)
    dtr <- 0.001657*sin(g)#s
    dtrd = dtr / DAYSEC
    if ( tt1 > tt2 ) {
        tdb1 = tt1
        tdb2 = tt2 + dtrd
    } else {
        tdb1 <- tt1 + dtrd
        tdb2 <- tt2
    }
    c(tdb1,tdb2)
}

##Delta(AT) (=TAI-UTC) for a given UTC date
sofa_Dat <- function(g){
####################################
##  For a given UTC date, calculate Delta(AT) = TAI-UTC.
##
##     :------------------------------------------:
##     :                                          :
##     :                 IMPORTANT                :
##     :                                          :
##     :  A new version of this function must be  :
##     :  produced whenever a new leap second is  :
##     :  announced.  There are four items to     :
##     :  change on each such occasion:           :
##     :                                          :
##     :  1) A new line must be added to the set  :
##     :     of statements that initialize the    :
##     :     array "changes".                     :
##     :                                          :
##     :  2) The constant IYV must be set to the  :
##     :     current year.                        :
##     :                                          :
##     :  3) The "Latest leap second" comment     :
##     :     below must be set to the new leap    :
##     :     second date.                         :
##     :                                          :
##     :  4) The "This revision" comment, later,  :
##     :     must be set to the current date.     :
##     :                                          :
##     :  Change (2) must also be carried out     :
##     :  whenever the function is re-issued,     :
##     :  even if no leap seconds have been       :
##     :  added.                                  :
##     :                                          :
##     :  Latest leap second:  2016 December 31   :
##     :                                          :
##     :__________________________________________:
##
##  This function is part of the International Astronomical Union's
##  SOFA (Standards Of Fundamental Astronomy) software collection.
##
##  Status:  user-replaceable support function.
##
##  Given:
##     iy     int      UTC:  year (Notes 1 and 2)
##     im     int            month (Note 2)
##     id     int            day (Notes 2 and 3)
##     fd     double         fraction of day (Note 4)
##
##  Returned:
##     deltat double   TAI minus UTC, seconds
##
##  Returned (function value):
##            int      status (Note 5):
##                       1 = dubious year (Note 1)
##                       0 = OK
##                      -1 = bad year
##                      -2 = bad month
##                      -3 = bad day (Note 3)
##                      -4 = bad fraction (Note 4)
##                      -5 = internal error (Note 5)
##
##  Notes:
##
##  1) UTC began at 1960 January 1.0 (JD 2436934.5) and it is improper
##     to call the function with an earlier date.  If this is attempted,
##     zero is returned together with a warning status.
##
##     Because leap seconds cannot, in principle, be predicted in
##     advance, a reliable check for dates beyond the valid range is
##     impossible.  To guard against gross errors, a year five or more
##     after the release year of the present function (see the constant
##     IYV) is considered dubious.  In this case a warning status is
##     returned but the result is computed in the normal way.
##
##     For both too-early and too-late years, the warning status is +1.
##     This is distinct from the error status -1, which signifies a year
##     so early that JD could not be computed.
##
##  2) If the specified date is for a day which ends with a leap second,
##     the TAI-UTC value returned is for the period leading up to the
##     leap second.  If the date is for a day which begins as a leap
##     second ends, the TAI-UTC returned is for the period following the
##     leap second.
##
##  3) The day number must be in the normal calendar range, for example
##     1 through 30 for April.  The "almanac" convention of allowing
##     such dates as January 0 and December 32 is not supported in this
##     function, in order to avoid confusion near leap seconds.
##
##  4) The fraction of day is used only for dates before the
##     introduction of leap seconds, the first of which occurred at the
##     end of 1971.  It is tested for validity (0 to 1 is the valid
##     range) even if not used;  if invalid, zero is used and status -4
##     is returned.  For many applications, setting fd to zero is
##     acceptable;  the resulting error is always less than 3 ms (and
##     occurs only pre-1972).
##
##  5) The status value returned in the case where there are multiple
##     errors refers to the first error detected.  For example, if the
##     month and day are 13 and 32 respectively, status -2 (bad month)
##     will be returned.  The "internal error" status refers to a
##     case that is impossible but causes some compilers to issue a
##     warning.
##
##  6) In cases where a valid result is not available, zero is returned.
##
##  References:
##
##  1) For dates from 1961 January 1 onwards, the expressions from the
##     file ftp://maia.usno.navy.mil/ser7/tai-utc.dat are used.
##
##  2) The 5ms timestep at 1961 January 1 is taken from 2.58.1 (p87) of
##     the 1992 Explanatory Supplement.
##
##  Called:
##     iauCal2jd    Gregorian calendar to JD
##
##  This revision:  2019 July 5
##
##  SOFA release 2019-07-22
##
##  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
####################################
    IYV <- 2019
    iy <- g[,1]
    im <- g[,2]
    id <- g[,3]
    fd <- g[,4]
##Reference dates (MJD) and drift rates (s/day), pre leap seconds
   drift <-rbind(
      c(37300.0, 0.0012960 ),
      c( 37300.0, 0.0012960 ),
      c( 37300.0, 0.0012960 ),
      c( 37665.0, 0.0011232 ),
      c( 37665.0, 0.0011232 ),
      c( 38761.0, 0.0012960 ),
      c( 38761.0, 0.0012960 ),
      c( 38761.0, 0.0012960 ),
      c( 38761.0, 0.0012960 ),
      c( 38761.0, 0.0012960 ),
      c( 38761.0, 0.0012960 ),
      c( 38761.0, 0.0012960 ),
      c( 39126.0, 0.0025920 ),
      c( 39126.0, 0.0025920 )
   )

##Number of Delta(AT) expressions before leap seconds were introduced */
   NERA1 <- nrow(drift)

##Number of Delta(AT) changes
   NDAT <- nrow(delat)

## Miscellaneous local variables
#   int j, i, m;
#   double da, djm0, djm;

##Initialize the result to zero.
   Dt <- rep(0,nrow(g))

##Convert the date into an MJD.
   jd <- time_Cal2JD(cbind(iy, im, id))
   djm0 <- DJM0
   djm <- (jd[,1]-DJM0)+jd[,2]

#If pre-UTC year, set warning status and give up.
   if(any(iy < delat[1,1])){
       cat('pre-UTC year ',delat[1,1],'!\n')
   }

# If suspiciously late year, set warning status but proceed.
   if(any(iy > IYV + 5)){
       cat('Warning: suspiciously late years after ',IYV + 5,'and TAI-UTC cannot be accurately determined!\n')
   }

#Combine year and month to form a date-ordered integer...
#the original SOFA routine only use integer months, I add fraction months to account for the recent dates in the delat table.
   m <- 12*iy + im

# Get the Delta(AT)
    ms <- 12*delat[,1] + delat[,2]
    mmax <- max(ms)
    mmin <- min(ms)
    Dt[m>mmax] <- delat[nrow(delat),3]#later UTC-TAI does not change

#...and use it to find the preceding table entry.
    ind0 <- which(m>=mmin & m<=mmax)#input epoch index
    inds <- unlist(sapply(m[ind0], function(x) which(x<ms)[1]-1))#list index
    inds[is.na(inds)] <- length(ms)#index for late dates
    Dt[ind0] <- delat[inds,3]#UTC-TAI for 1960-IYV

# If pre-1972, adjust for drift
# m: ind0; ms: inds
    ind1 <- which(m[ind0]<ms[NERA1+1] & m[ind0]>=mmin)
    if(length(ind1)>0){
        ind.m <- ind0[ind1]#input epoch index
        ind.ms <- inds[ind1]#drift index
        Dt[ind.m] <- Dt[ind.m] + (djm[ind.m] + fd[ind.m] - drift[ind.ms,1])*drift[ind.ms,2]
    }
# for pre-1960 dates, zero difference is used.
    index <- which(m<mmin)
    if(length(index)>0){
        cat('Warning: ',length(index),' epochs earlier than 1960 January 1.0 (JD 2436934.5) when UTC began and zero UTC-TAI drift is used.\n')
        if(length(Dt)>length(index)){
            Dt[index] <- min(Dt[-index])
        }else{
            Dt[index] <- delat[1,3]+(time_Cal2JD(cbind(1960, 1, 1))[,2]-drift[1,1])*drift[1,2]
        }
    }
   return(Dt)
}

sofa_Tttai <- function(utc,tt,Par){
####################################
##  Time scale transformation:  Terrestrial Time, TT, to International
##  Atomic Time, TAI.
##
##  This function is part of the International Astronomical Union's
##  SOFA (Standards of Fundamental Astronomy) software collection.
##
##  Status:  canonical.
##
##  Given:
##     tt1,tt2    double    TT as a 2-part Julian Date
##
##  Returned:
##     tai1,tai2  double    TAI as a 2-part Julian Date
##
##  Returned (function value):
##                int       status:  0 = OK
##
##  Note:
##
##     tt1+tt2 is Julian Date, apportioned in any convenient way between
##     the two arguments, for example where tt1 is the Julian Day Number
##     and tt2 is the fraction of a day.  The returned tai1,tai2 follow
##     suit.
##
##  References:
##
##     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
##     IERS Technical Note No. 32, BKG (2004)
##
##     Explanatory Supplement to the Astronomical Almanac,
##     P. Kenneth Seidelmann (ed), University Science Books (1992)
##
##  This revision:  2019 June 20
##
##  SOFA release 2019-07-22
##
##  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
####################################
##tt
    if(Par$TtType=='BIPM'){
        dt <- time_BipmCorr(utc,BIPM=Par$BIPM)
    }else{
        dt <- 0
    }
##TT minus TAI (days).
    dtat <- TTMTAI/DAYSEC

##Result, safeguarding precision.
    tai1 <- tt1
    tai2 <- tt2 - dtat-dt

    return(cbind(tai1,tai2))
}

###Taitt is a modified version of iauTaitt in SOFA to account for the BIPMXX realization of TT.
sofa_Taitt <- function(utc,tai,Par){
####################################
##  Time scale transformation:  International Atomic Time, TAI, to
##  Terrestrial Time, TT.
##
##  This function is part of the International Astronomical Union's
##  SOFA (Standards of Fundamental Astronomy) software collection.
##
##  Status:  canonical.
##
##  Given:
##     tai1,tai2  double    TAI as a 2-part Julian Date
##
##  Returned:
##     tt1,tt2    double    TT as a 2-part Julian Date
##
##  Returned (function value):
##                int       status:  0 = OK
##
##  Note:
##
##     tai1+tai2 is Julian Date, apportioned in any convenient way
##     between the two arguments, for example where tai1 is the Julian
##     Day Number and tai2 is the fraction of a day.  The returned
##     tt1,tt2 follow suit.
##
##  References:
##
##     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
##     IERS Technical Note No. 32, BKG (2004)
##
##     Explanatory Supplement to the Astronomical Almanac,
##     P. Kenneth Seidelmann (ed), University Science Books (1992)
##
##  This revision:  2019 June 20
##
##  SOFA release 2019-07-22
##
##  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
####################################
    if(Par$TtType=='BIPM'){
        dt <- time_BipmCorr((utc[,1]-DJM0)+utc[,2],BIPM=Par$BIPM)/DAYSEC
    }else{
        dt <- 0
    }
    dtat <-  TTMTAI/DAYSEC
#Result, safeguarding precision.
    cbind(tai[,1],tai[,2] + dtat+dt)
}

sofa_Utctai <- function(utc,TaiType='scale'){
####################################
    ## This function is adapted from iauUtctai in the SOFA libraries
    ## While iauUtctai only implement the case of TaiType=scale (treat JD day as 86399 or 86401 days if leap seconds happended), sofa_Utctai implements both scaled and instant cases.
    ## ref: http://www.iausofa.org/2019_0722_C/sofa/utctai.c
    ##  Time scale transformation:  Coordinated Universal Time, UTC, to
    ##  International Atomic Time, TAI.
    ##
    ##  This function is part of the International Astronomical Union's
    ##  SOFA (Standards of Fundamental Astronomy) software collection.
    ##
    ##  Status:  canonical.
    ##
    ##  Given:
    ##     utc1,utc2  double   UTC as a 2-part quasi Julian Date (Notes 1-4)
    ##
    ##  Returned:
    ##     tai1,tai2  double   TAI as a 2-part Julian Date (Note 5)
    ##
    ##  Returned (function value):
    ##                int      status: +1 = dubious year (Note 3)
    ##                                  0 = OK
    ##                                 -1 = unacceptable date
    ##
    ##  Notes:
    ##
    ##  1) utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
    ##     convenient way between the two arguments, for example where utc1
    ##     is the Julian Day Number and utc2 is the fraction of a day.
    ##
    ##  2) JD cannot unambiguously represent UTC during a leap second unless
    ##     special measures are taken.  The convention in the present
    ##     function is that the JD day represents UTC days whether the
    ##     length is 86399, 86400 or 86401 SI seconds.  In the 1960-1972 era
    ##     there were smaller jumps (in either direction) each time the
    ##     linear UTC(TAI) expression was changed, and these "mini-leaps"
    ##     are also included in the SOFA convention.
    ##
    ##  3) The warning status "dubious year" flags UTCs that predate the
    ##     introduction of the time scale or that are too far in the future
    ##     to be trusted.  See iauDat for further details.
    ##
    ##  4) The function iauDtf2d converts from calendar date and time of day
    ##     into 2-part Julian Date, and in the case of UTC implements the
    ##     leap-second-ambiguity convention described above.
    ##
    ##  5) The returned TAI1,TAI2 are such that their sum is the TAI Julian
    ##     Date.
    ##
    ##  Called:
    ##     iauJd2cal    JD to Gregorian calendar
    ##     iauDat       delta(AT) = TAI-UTC
    ##     iauCal2jd    Gregorian calendar to JD
    ##
    ##  References:
    ##
    ##     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    ##     IERS Technical Note No. 32, BKG (2004)
    ##
    ##     Explanatory Supplement to the Astronomical Almanac,
    ##     P. Kenneth Seidelmann (ed), University Science Books (1992)
    ##
    ##  This revision:  2019 June 20
    ##
    ##  SOFA release 2019-07-22
    ##
    ##  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
####################################
                                        #utc1/u1 should always be larger than utc2/u2
    j <- time_Jd2cal(utc)
    u1 <- utc[,1]
    u2 <- utc[,2]
    iy <- j[,1]
    im <- j[,2]
    id <- j[,3]
    fd <- j[,4]
    if(TaiType=='instant'){
        dat0 <- sofa_Dat(j[,1:4,drop=FALSE])
        dat.past <- sofa_Dat(cbind(j[,1:3,drop=FALSE],j[,4]-1e-6))
        dat.future <- sofa_Dat(cbind(j[,1:3,drop=FALSE],j[,4]+1e-6))
        dlod <- 2.0 * (dat0 - dat.past)
        dleap <- dat.future - (dat.past + dlod)
        fd <- fd+dleap/DAYSEC
        fd <- fd+dlod/DAYSEC
    }else{
        ## Get TAI-UTC at 0h at given epochs.
        dat0 <- sofa_Dat(cbind(j[,1:3,drop=FALSE],t(t(rep(0,nrow(j))))))

        ##Get TAI-UTC at 12h (to detect drift).
        dat12 <- sofa_Dat(cbind(j[,1:3,drop=FALSE],t(t(rep(0.5,nrow(j))))))

###Get TAI-UTC at 0h tomorrow (to detect jumps).
###Note that the original sofa code at http://www.iausofa.org/2018_0130_C/sofa/utctai.c has a typo the following commented function would yield 12h tomorrow
        ##    j <- time_Jd2cal(cbind(u1+1.5, u2-fd))
        j <- time_Jd2cal(cbind(u1+1.5, u2-fd))
        iyt <- j[,1]
        imt <- j[,2]
        idt <- j[,3]
        w <- j[,4]
        dat24 <- sofa_Dat(cbind(j[,1:3,drop=FALSE],t(t(rep(0,nrow(j))))))

        ## Separate TAI-UTC change into per-day (DLOD) and any jump (DLEAP). */
        dlod <- 2.0 * (dat12 - dat0)
        dleap <- dat24 - (dat0 + dlod)

        ##Remove any scaling applied to spread leap into preceding day.
        fd <- fd*(DAYSEC+dleap)/DAYSEC

        ## Scale from (pre-1972) UTC seconds to SI seconds.
        fd <- fd*(DAYSEC+dlod)/DAYSEC
    }
    ##given epoch calendar date to 2-part JD.
    tmp <- time_Cal2JD(cbind(iy, im, id))
    z1 <- tmp[,1]
    z2 <- tmp[,2]

    ## Assemble the TAI result, preserving the UTC split and order.
    a2 <- z1 - u1
    a2 <- a2+z2
    a2 <- a2+fd + dat0/DAYSEC
    tai <- cbind(u1,a2)
    return(list(tai=tai,leap=dleap))
}

sofa_Taiut1 <- function(tai1,tai2,dta){
####################################
##  Time scale transformation:  International Atomic Time, TAI, to
##  Universal Time, UT1.
##
##  This function is part of the International Astronomical Union's
##  SOFA (Standards of Fundamental Astronomy) software collection.
##
##  Status:  canonical.
##
##  Given:
##     tai1,tai2  double    TAI as a 2-part Julian Date
##     dta        double    UT1-TAI in seconds
##
##  Returned:
##     ut11,ut12  double    UT1 as a 2-part Julian Date
##
##  Returned (function value):
##                int       status:  0 = OK
##
##  Notes:
##
##  1) tai1+tai2 is Julian Date, apportioned in any convenient way
##     between the two arguments, for example where tai1 is the Julian
##     Day Number and tai2 is the fraction of a day.  The returned
##     UT11,UT12 follow suit.
##
##  2) The argument dta, i.e. UT1-TAI, is an observed quantity, and is
##     available from IERS tabulations.
##
##  Reference:
##
##     Explanatory Supplement to the Astronomical Almanac,
##     P. Kenneth Seidelmann (ed), University Science Books (1992)
##
##  This revision:  2019 June 20
##
##  SOFA release 2019-07-22
##
##  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
####################################
   dtad <- dta/DAYSEC
   ut11 <- tai1
   ut12 <- tai2 + dtad
   return(cbind(ut11,ut12))
}

sofa_Taiutc <- function(tai,TaiType){
####################################
##  Time scale transformation:  International Atomic Time, TAI, to
##  Coordinated Universal Time, UTC.
##
##  This function is part of the International Astronomical Union's
##  SOFA (Standards of Fundamental Astronomy) software collection.
##
##  Status:  canonical.
##
##  Given:
##     tai1,tai2  double   TAI as a 2-part Julian Date (Note 1)
##
##  Returned:
##     utc1,utc2  double   UTC as a 2-part quasi Julian Date (Notes 1-3)
##
##  Returned (function value):
##                int      status: +1 = dubious year (Note 4)
##                                  0 = OK
##                                 -1 = unacceptable date
##
##  Notes:
##
##  1) tai1+tai2 is Julian Date, apportioned in any convenient way
##     between the two arguments, for example where tai1 is the Julian
##     Day Number and tai2 is the fraction of a day.  The returned utc1
##     and utc2 form an analogous pair, except that a special convention
##     is used, to deal with the problem of leap seconds - see the next
##     note.
##
##  2) JD cannot unambiguously represent UTC during a leap second unless
##     special measures are taken.  The convention in the present
##     function is that the JD day represents UTC days whether the
##     length is 86399, 86400 or 86401 SI seconds.  In the 1960-1972 era
##     there were smaller jumps (in either direction) each time the
##     linear UTC(TAI) expression was changed, and these "mini-leaps"
##     are also included in the SOFA convention.
##
##  3) The function iauD2dtf can be used to transform the UTC quasi-JD
##     into calendar date and clock time, including UTC leap second
##     handling.
##
##  4) The warning status "dubious year" flags UTCs that predate the
##     introduction of the time scale or that are too far in the future
##     to be trusted.  See iauDat for further details.
##
##  Called:
##     iauUtctai    UTC to TAI
##
##  References:
##
##     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
##     IERS Technical Note No. 32, BKG (2004)
##
##     Explanatory Supplement to the Astronomical Almanac,
##     P. Kenneth Seidelmann (ed), University Science Books (1992)
##
##  This revision:  2019 June 20
##
##  SOFA release 2019-07-22
##
##  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
####################################
## Put the two parts of the TAI into big-first order.
    a1 <- tai[,1]
    a2 <- tai[,2]
##Initial guess for UTC.
   u1 <- a1
   u2 <- a2

##Iterate (though in most cases just once is enough).
   for(i in 1:3){
##Guessed UTC to TAI.
      g <- sofa_Utctai(cbind(u1, u2),TaiType=TaiType)
##Adjust guessed UTC.
      u2 <- u2+a1 - g[,1]
      u2 <- u2+a2 - g[,2]
   }

##Return the UTC result, preserving the TAI order.
      utc1 <- u1
      utc2 <- u2
##return
   return(c(u1,u2))
}

sofa_Utcut1 <- function(utc,tempo=TRUE){
####################################
##  Time scale transformation:  Coordinated Universal Time, UTC, to
##  Universal Time, UT1.
##
##  This function is part of the International Astronomical Union's
##  SOFA (Standards of Fundamental Astronomy) software collection.
##
##  Status:  canonical.
##
##  Given:
##     utc1,utc2  double   UTC as a 2-part quasi Julian Date (Notes 1-4)
##     dut1       double   Delta UT1 = UT1-UTC in seconds (Note 5)
##
##  Returned:
##     ut11,ut12  double   UT1 as a 2-part Julian Date (Note 6)
##
##  Returned (function value):
##                int      status: +1 = dubious year (Note 3)
##                                  0 = OK
##                                 -1 = unacceptable date
##
##  Notes:
##
##  1) utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
##     convenient way between the two arguments, for example where utc1
##     is the Julian Day Number and utc2 is the fraction of a day.
##
##  2) JD cannot unambiguously represent UTC during a leap second unless
##     special measures are taken.  The convention in the present
##     function is that the JD day represents UTC days whether the
##     length is 86399, 86400 or 86401 SI seconds.
##
##  3) The warning status "dubious year" flags UTCs that predate the
##     introduction of the time scale or that are too far in the future
##     to be trusted.  See iauDat for further details.
##
##  4) The function iauDtf2d converts from calendar date and time of
##     day into 2-part Julian Date, and in the case of UTC implements
##     the leap-second-ambiguity convention described above.
##
##  5) Delta UT1 can be obtained from tabulations provided by the
##     International Earth Rotation and Reference Systems Service.
##     It is the caller's responsibility to supply a dut1 argument
##     containing the UT1-UTC value that matches the given UTC.
##
##  6) The returned ut11,ut12 are such that their sum is the UT1 Julian
##     Date.
##
##  References:
##
##     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
##     IERS Technical Note No. 32, BKG (2004)
##
##     Explanatory Supplement to the Astronomical Almanac,
##     P. Kenneth Seidelmann (ed), University Science Books (1992)
##
##  Called:
##     iauJd2cal    JD to Gregorian calendar
##     iauDat       delta(AT) = TAI-UTC
##     iauUtctai    UTC to TAI
##     iauTaiut1    TAI to UT1
##
##  This revision:  2013 August 12
##
##  SOFA release 2019-07-22
##
##  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
####################################
    utc1 <- utc[,1]
    utc2 <- utc[,2]
    ##Look up TAI-UTC.
    tmp <- time_Jd2cal(utc)
    iy <- tmp[,1]
    im <- tmp[,2]
    id <- tmp[,3]
    w <- tmp[,4]
    dat <- sofa_Dat(cbind(iy, im, id, 0.0))

    ## Form UT1-TAI.
    mjd <- (utc1-DJM0)+utc2
    if(tempo){
        dut1 <- calDut1dot(mjd)$dut1
    }else{
        dut1 <- time_Ut1mUtc(mjd)
    }
    dta <- dut1 - dat

    ## UTC to TAI to UT1. */
    tai <- sofa_Utctai(utc,TaiType=Par$TaiType)
    Taiut1(tai[,1], tai[,2], dta)
}

sofa_Eform <- function(n){
####################################
##  Earth reference ellipsoids.
##
##  This function is part of the International Astronomical Union's
##  SOFA (Standards of Fundamental Astronomy) software collection.
##
##  Status:  canonical.
##
##  Given:
##     n    int         ellipsoid identifier (Note 1)
##
##  Returned:
##     a    double      equatorial radius (meters, Note 2)
##     f    double      flattening (Note 2)
##
##  Returned (function value):
##          int         status:  0 = OK
##                              -1 = illegal identifier (Note 3)
##
##  Notes:
##
##  1) The identifier n is a number that specifies the choice of
##     reference ellipsoid.  The following are supported:
##
##        n    ellipsoid
##
##        1     WGS84
##        2     GRS80
##        3     WGS72
##
##     The n value has no significance outside the SOFA software.  For
##     convenience, symbols WGS84 etc. are defined in sofam.h.
##
##  2) The ellipsoid parameters are returned in the form of equatorial
##     radius in meters (a) and flattening (f).  The latter is a number
##     around 0.00335, i.e. around 1/298.
##
##  3) For the case where an unsupported n value is supplied, zero a and
##     f are returned, as well as error status.
##
##  References:
##
##     Department of Defense World Geodetic System 1984, National
##     Imagery and Mapping Agency Technical Report 8350.2, Third
##     Edition, p3-2.
##
##     Moritz, H., Bull. Geodesique 66-2, 187 (1992).
##
##     The Department of Defense World Geodetic System 1972, World
##     Geodetic System Committee, May 1974.
##
##     Explanatory Supplement to the Astronomical Almanac,
##     P. Kenneth Seidelmann (ed), University Science Books (1992),
##     p220.
##
##  This revision:  2013 June 18
##
##  SOFA release 2019-07-22
##
##  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
####################################
    if(n==1){
        a <- 6378137.0
        f <- 1.0 / 298.257223563
    }

    if(n==2){
        a = 6378137.0
        f = 1.0 / 298.257222101
    }

    if(n==3){
        a = 6378135.0
        f = 1.0 / 298.26
    }
    return(c(a,f))
   }

sofa_Gd2gce <- function(a,f,elong, phi, height){
####################################
##  Transform geodetic coordinates to geocentric for a reference
##  ellipsoid of specified form.
##
##  This function is part of the International Astronomical Union's
##  SOFA (Standards of Fundamental Astronomy) software collection.
##
##  Status:  support function.
##
##  Given:
##     a       double     equatorial radius (Notes 1,4)
##     f       double     flattening (Notes 2,4)
##     elong   double     longitude (radians, east +ve)
##     phi     double     latitude (geodetic, radians, Note 4)
##     height  double     height above ellipsoid (geodetic, Notes 3,4)
##
##  Returned:
##     xyz     double[3]  geocentric vector (Note 3)
##
##  Returned (function value):
##             int        status:  0 = OK
##                                -1 = illegal case (Note 4)
##  Notes:
##
##  1) The equatorial radius, a, can be in any units, but meters is
##     the conventional choice.
##
##  2) The flattening, f, is (for the Earth) a value around 0.00335,
##     i.e. around 1/298.
##
##  3) The equatorial radius, a, and the height, height, must be
##     given in the same units, and determine the units of the
##     returned geocentric vector, xyz.
##
##  4) No validation is performed on individual arguments.  The error
##     status -1 protects against (unrealistic) cases that would lead
##     to arithmetic exceptions.  If an error occurs, xyz is unchanged.
##
##  5) The inverse transformation is performed in the function
##     iauGc2gde.
##
##  6) The transformation for a standard ellipsoid (such as WGS84) can
##     more conveniently be performed by calling iauGd2gc,  which uses a
##     numerical code to identify the required a and f values.
##
##  References:
##
##     Green, R.M., Spherical Astronomy, Cambridge University Press,
##     (1985) Section 4.5, p96.
##
##     Explanatory Supplement to the Astronomical Almanac,
##     P. Kenneth Seidelmann (ed), University Science Books (1992),
##     Section 4.22, p202.
##
##  This revision:  2013 June 18
##
##  SOFA release 2019-07-22
##
##  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
####################################
   sp = sin(phi)
   cp = cos(phi)
   w = 1.0 - f
   w = w * w
   d = cp*cp + w*sp*sp
   if( d <= 0.0 ) return(NULL)
   ac = a / sqrt(d)
   as = w * ac

## Geocentric vector.
    r = (ac + height)*cp
    x = r * cos(elong)
    y = r * sin(elong)
    z = (as + height)*sp
    return(cbind(x,y,z))
}

sofa_Gd2gc <- function(n, elong, phi, height){
####################################
##  Transform geodetic coordinates to geocentric using the specified
##  reference ellipsoid.
##
##  This function is part of the International Astronomical Union's
##  SOFA (Standards of Fundamental Astronomy) software collection.
##
##  Status:  canonical transformation.
##
##  Given:
##     n       int        ellipsoid identifier (Note 1)
##     elong   double     longitude (radians, east +ve)
##     phi     double     latitude (geodetic, radians, Note 3)
##     height  double     height above ellipsoid (geodetic, Notes 2,3)
##
##  Returned:
##     xyz     double[3]  geocentric vector (Note 2)
##
##  Returned (function value):
##             int        status:  0 = OK
##                                -1 = illegal identifier (Note 3)
##                                -2 = illegal case (Note 3)
##
##  Notes:
##
##  1) The identifier n is a number that specifies the choice of
##     reference ellipsoid.  The following are supported:
##
##        n    ellipsoid
##
##        1     WGS84
##        2     GRS80
##        3     WGS72
##
##     The n value has no significance outside the SOFA software.  For
##     convenience, symbols WGS84 etc. are defined in sofam.h.
##
##  2) The height (height, given) and the geocentric vector (xyz,
##     returned) are in meters.
##
##  3) No validation is performed on the arguments elong, phi and
##     height.  An error status -1 means that the identifier n is
##     illegal.  An error status -2 protects against cases that would
##     lead to arithmetic exceptions.  In all error cases, xyz is set
##     to zeros.
##
##  4) The inverse transformation is performed in the function iauGc2gd.
##
##  Called:
##     iauEform     Earth reference ellipsoids
##     iauGd2gce    geodetic to geocentric transformation, general
##     iauZp        zero p-vector
##
##  This revision:  2013 June 18
##
##  SOFA release 2019-07-22
##
##  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
####################################

# Obtain reference ellipsoid parameters.
    af <- sofa_Eform(n)

# transform longitude, geodetic latitude, height to x,y,z.
    sofa_Gd2gce(af[1], af[2], elong, phi, height)
}
sofa_Tdbtcb <- function(tdb){
####################################
##  Time scale transformation:  Barycentric Dynamical Time, TDB, to
##  Barycentric Coordinate Time, TCB.
##
##  This function is part of the International Astronomical Union's
##  SOFA (Standards of Fundamental Astronomy) software collection.
##
##  Status:  canonical.
##
##  Given:
##     tdb1,tdb2  double    TDB as a 2-part Julian Date
##
##  Returned:
##     tcb1,tcb2  double    TCB as a 2-part Julian Date
##
##  Returned (function value):
##                int       status:  0 = OK
##
##  Notes:
##
##  1) tdb1+tdb2 is Julian Date, apportioned in any convenient way
##     between the two arguments, for example where tdb1 is the Julian
##     Day Number and tdb2 is the fraction of a day.  The returned
##     tcb1,tcb2 follow suit.
##
##  2) The 2006 IAU General Assembly introduced a conventional linear
##     transformation between TDB and TCB.  This transformation
##     compensates for the drift between TCB and terrestrial time TT,
##     and keeps TDB approximately centered on TT.  Because the
##     relationship between TT and TCB depends on the adopted solar
##     system ephemeris, the degree of alignment between TDB and TT over
##     long intervals will vary according to which ephemeris is used.
##     Former definitions of TDB attempted to avoid this problem by
##     stipulating that TDB and TT should differ only by periodic
##     effects.  This is a good description of the nature of the
##     relationship but eluded precise mathematical formulation.  The
##     conventional linear relationship adopted in 2006 sidestepped
##     these difficulties whilst delivering a TDB that in practice was
##     consistent with values before that date.
##
##  3) TDB is essentially the same as Teph, the time argument for the
##     JPL solar system ephemerides.
##
##  Reference:
##
##     IAU 2006 Resolution B3
##
##  This revision:  2019 June 20
##
##  SOFA release 2019-07-22
##
##  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
####################################

    ##1977 Jan 1 00:00:32.184 TT, as two-part JD
    t77td = DJM0 + DJM77
    t77tf = TTMTAI/DAYSEC

    ##TDB (days) at TAI 1977 Jan 1.0
    tdb0 = TDB0/DAYSEC

    ## TDB to TCB rate
    elbb = ELB/(1.0-ELB)

    tdb1 <- tdb[,1]
    tdb2 <- tdb[,2]
    ##Result, preserving date format but safeguarding precision.
    d <- t77td - tdb1
    f <- tdb2 - tdb0
    tcb1 <- tdb1
    tcb2 <- f - ( d - ( f - t77tf ) ) * elbb
    cbind(tcb1,tcb2)
}

sofa_Fal03 <- function(t){
####################################
##  Fundamental argument, IERS Conventions (2003):
##  mean anomaly of the Moon.
##
##  This function is part of the International Astronomical Union's
##  SOFA (Standards Of Fundamental Astronomy) software collection.
##
##  Status:  canonical model.
##
##  Given:
##     t     double    TDB, Julian centuries since J2000.0 (Note 1)
##
##  Returned (function value):
##           double    l, radians (Note 2)
##
##  Notes:
##
##  1) Though t is strictly TDB, it is usually more convenient to use
##     TT, which makes no significant difference.
##
##  2) The expression used is as adopted in IERS Conventions (2003) and
##     is from Simon et al. (1994).
##
##  References:
##
##     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
##     IERS Technical Note No. 32, BKG (2004)
##
##     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
##     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
##
##  This revision:  2013 June 18
##
##  SOFA release 2019-07-22
##
##  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
####################################
          ((485868.249036  +
             t * ( 1717915923.2178 +
             t * (         31.8792 +
             t * (          0.051635 +
             t * (        -0.00024470 ) ) ) ))%%TURNAS ) * DAS2R
}

sofa_Falp03 <- function(t){
####################################
##  Fundamental argument, IERS Conventions (2003):
##  mean anomaly of the Sun.
##
##  This function is part of the International Astronomical Union's
##  SOFA (Standards Of Fundamental Astronomy) software collection.
##
##  Status:  canonical model.
##
##  Given:
##     t     double    TDB, Julian centuries since J2000.0 (Note 1)
##
##  Returned (function value):
##           double    l', radians (Note 2)
##
##  Notes:
##
##  1) Though t is strictly TDB, it is usually more convenient to use
##     TT, which makes no significant difference.
##
##  2) The expression used is as adopted in IERS Conventions (2003) and
##     is from Simon et al. (1994).
##  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
####################################
        ( ( 1287104.793048 +
             t * ( 129596581.0481 +
             t * (       -0.5532 +
             t * (         0.000136 +
             t * (       -0.00001149 ) ) ) ))%%TURNAS)*DAS2R
}

sofa_Fad03 <- function(t){
####################################
##  Fundamental argument, IERS Conventions (2003):
##  mean elongation of the Moon from the Sun.
##
##  This function is part of the International Astronomical Union's
##  SOFA (Standards Of Fundamental Astronomy) software collection.
##
##  Status:  canonical model.
##
##  Given:
##     t     double    TDB, Julian centuries since J2000.0 (Note 1)
##
##  Returned (function value):
##           double    D, radians (Note 2)
##
##  Notes:
##
##  1) Though t is strictly TDB, it is usually more convenient to use
##     TT, which makes no significant difference.
##
##  2) The expression used is as adopted in IERS Conventions (2003) and
##     is from Simon et al. (1994).
##
##  References:
##
##     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
##     IERS Technical Note No. 32, BKG (2004)
##
##     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
##     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
##
##  This revision:  2013 June 18
##
##  SOFA release 2019-07-22
##
##  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
####################################

          ((1072260.703692 +
             t * ( 1602961601.2090 +
             t * (        - 6.3706 +
             t * (          0.006593 +
             t * (        - 0.00003169 ) ) ) ))%%TURNAS)*DAS2R
}

sofa_Faf03 <- function(t){
####################################
##  Fundamental argument, IERS Conventions (2003):
##  mean longitude of the Moon minus mean longitude of the ascending
##  node.
##
##  This function is part of the International Astronomical Union's
##  SOFA (Standards Of Fundamental Astronomy) software collection.
##
##  Status:  canonical model.
##
##  Given:
##     t     double    TDB, Julian centuries since J2000.0 (Note 1)
##
##  Returned (function value):
##           double    F, radians (Note 2)
##
##  Notes:
##
##  1) Though t is strictly TDB, it is usually more convenient to use
##     TT, which makes no significant difference.
##
##  2) The expression used is as adopted in IERS Conventions (2003) and
##     is from Simon et al. (1994).
##
##  References:
##
##     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
##     IERS Technical Note No. 32, BKG (2004)
##
##     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
##     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
##
##  This revision:  2013 June 18
##
##  SOFA release 2019-07-22
##
##  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
####################################
    ((335779.526232 +
             t * ( 1739527262.8478 +
             t * (       - 12.7512 +
             t * (        - 0.001037 +
             t * (          0.00000417 ) ) ) ))%%TURNAS)*DAS2R
}

sofa_Fave03 <- function(t){
####################################
##  Fundamental argument, IERS Conventions (2003):
##  mean longitude of Venus.
##
##  This function is part of the International Astronomical Union's
##  SOFA (Standards Of Fundamental Astronomy) software collection.
##
##  Status:  canonical model.
##
##  Given:
##     t     double    TDB, Julian centuries since J2000.0 (Note 1)
##
##  Returned (function value):
##           double    mean longitude of Venus, radians (Note 2)
##
##  Notes:
##
##  1) Though t is strictly TDB, it is usually more convenient to use
##     TT, which makes no significant difference.
##
##  2) The expression used is as adopted in IERS Conventions (2003) and
##     comes from Souchay et al. (1999) after Simon et al. (1994).
##
##  References:
##
##     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
##     IERS Technical Note No. 32, BKG (2004)
##
##     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
##     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
##
##     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
##     Astron.Astrophys.Supp.Ser. 135, 111
##
##  This revision:  2013 June 18
##
##  SOFA release 2019-07-22
##
##  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
####################################
    (3.176146697 + 1021.3285546211 * t)%%D2PI
}

sofa_Faom03 <- function(t){
####################################
##  Fundamental argument, IERS Conventions (2003):
##  mean longitude of the Moon's ascending node.
##
##  This function is part of the International Astronomical Union's
##  SOFA (Standards Of Fundamental Astronomy) software collection.
##
##  Status:  canonical model.
##
##  Given:
##     t     double    TDB, Julian centuries since J2000.0 (Note 1)
##
##  Returned (function value):
##           double    Omega, radians (Note 2)
##
##  Notes:
##
##  1) Though t is strictly TDB, it is usually more convenient to use
##     TT, which makes no significant difference.
##
##  2) The expression used is as adopted in IERS Conventions (2003) and
##     is from Simon et al. (1994).
##
##  References:
##
##     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
##     IERS Technical Note No. 32, BKG (2004)
##
##     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
##     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
##
##  This revision:  2013 June 18
##
##  SOFA release 2019-07-22
##
##  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
####################################
    (( 450160.398036 +
             t * ( - 6962890.5431 +
             t * (         7.4722 +
             t * (         0.007702 +
             t * (       - 0.00005939 ) ) ) ))%%TURNAS ) * DAS2R;
}

sofa_Fae03 <- function(t){
####################################
##  Fundamental argument, IERS Conventions (2003):
##  mean longitude of Earth.
##
##  This function is part of the International Astronomical Union's
##  SOFA (Standards Of Fundamental Astronomy) software collection.
##
##  Status:  canonical model.
##
##  Given:
##     t     double    TDB, Julian centuries since J2000.0 (Note 1)
##
##  Returned (function value):
##           double    mean longitude of Earth, radians (Note 2)
##
##  Notes:
##
##  1) Though t is strictly TDB, it is usually more convenient to use
##     TT, which makes no significant difference.
##
##  2) The expression used is as adopted in IERS Conventions (2003) and
##     comes from Souchay et al. (1999) after Simon et al. (1994).
##
##  References:
##
##     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
##     IERS Technical Note No. 32, BKG (2004)
##
##     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
##     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
##
##     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
##     Astron.Astrophys.Supp.Ser. 135, 111
##
##  This revision:  2013 June 18
##
##  SOFA release 2019-07-22
##
##  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
####################################
    (1.753470314 + 628.3075849991 * t)%%D2PI
}

sofa_Fapa03 <- function(t){
####################################
##  Fundamental argument, IERS Conventions (2003):
##  general accumulated precession in longitude.
##
##  This function is part of the International Astronomical Union's
##  SOFA (Standards Of Fundamental Astronomy) software collection.
##
##  Status:  canonical model.
##
##  Given:
##     t     double    TDB, Julian centuries since J2000.0 (Note 1)
##
##  Returned (function value):
##           double    general precession in longitude, radians (Note 2)
##
##  Notes:
##
##  1) Though t is strictly TDB, it is usually more convenient to use
##     TT, which makes no significant difference.
##
##  2) The expression used is as adopted in IERS Conventions (2003).  It
##     is taken from Kinoshita & Souchay (1990) and comes originally
##     from Lieske et al. (1977).
##
##  References:
##
##     Kinoshita, H. and Souchay J. 1990, Celest.Mech. and Dyn.Astron.
##     48, 187
##
##     Lieske, J.H., Lederle, T., Fricke, W. & Morando, B. 1977,
##     Astron.Astrophys. 58, 1-16
##
##     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
##     IERS Technical Note No. 32, BKG (2004)
##
##  This revision:  2013 June 18
##
##  SOFA release 2019-07-22
##
##  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
####################################
    (0.024381750 + 0.00000538691 * t) * t
}
##
sofa_S06 <- function(tt,x,y){
####################################
##  The CIO locator s, positioning the Celestial Intermediate Origin on
##  the equator of the Celestial Intermediate Pole, given the CIP's X,Y
##  coordinates.  Compatible with IAU 2006/2000A precession-nutation.
##
##  This function is part of the International Astronomical Union's
##  SOFA (Standards Of Fundamental Astronomy) software collection.
##
##  Status:  canonical model.
##
##  Given:
##     date1,date2   double    TT as a 2-part Julian Date (Note 1)
##     x,y           double    CIP coordinates (Note 3)
##
##  Returned (function value):
##                   double    the CIO locator s in radians (Note 2)
##
##  Notes:
##
##  1) The TT date date1+date2 is a Julian Date, apportioned in any
##     convenient way between the two arguments.  For example,
##     JD(TT)=2450123.7 could be expressed in any of these ways,
##     among others:
##
##            date1          date2
##
##         2450123.7           0.0       (JD method)
##         2451545.0       -1421.3       (J2000 method)
##         2400000.5       50123.2       (MJD method)
##         2450123.5           0.2       (date & time method)
##
##     The JD method is the most natural and convenient to use in
##     cases where the loss of several decimal digits of resolution
##     is acceptable.  The J2000 method is best matched to the way
##     the argument is handled internally and will deliver the
##     optimum resolution.  The MJD method and the date & time methods
##     are both good compromises between resolution and convenience.
##
##  2) The CIO locator s is the difference between the right ascensions
##     of the same point in two systems:  the two systems are the GCRS
##     and the CIP,CIO, and the point is the ascending node of the
##     CIP equator.  The quantity s remains below 0.1 arcsecond
##     throughout 1900-2100.
##
##  3) The series used to compute s is in fact for s+XY/2, where X and Y
##     are the x and y components of the CIP unit vector;  this series
##     is more compact than a direct series for s would be.  This
##     function requires X,Y to be supplied by the caller, who is
##     responsible for providing values that are consistent with the
##     supplied date.
##
##  4) The model is consistent with the "P03" precession (Capitaine et
##     al. 2003), adopted by IAU 2006 Resolution 1, 2006, and the
##     IAU 2000A nutation (with P03 adjustments).
##
##  Called:
##     iauFal03     mean anomaly of the Moon
##     iauFalp03    mean anomaly of the Sun
##     iauFaf03     mean argument of the latitude of the Moon
##     iauFad03     mean elongation of the Moon from the Sun
##     iauFaom03    mean longitude of the Moon's ascending node
##     iauFave03    mean longitude of Venus
##     iauFae03     mean longitude of Earth
##     iauFapa03    general accumulated precession in longitude
##
##  References:
##
##     Capitaine, N., Wallace, P.T. & Chapront, J., 2003, Astron.
##     Astrophys. 432, 355
##
##     McCarthy, D.D., Petit, G. (eds.) 2004, IERS Conventions (2003),
##     IERS Technical Note No. 32, BKG
##
##  This revision:  2019 June 23
##
##  SOFA release 2019-07-22
##
##  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
####################################

    sp <-   c(94.00e-6,  3808.65e-6,  -122.68e-6,-72574.11e-6, 27.98e-6,15.62e-6)
    s0 <-  rbind(
        c( 0,  0,  0,  0,  1,  0,  0,  0, -2640.73e-6,   0.39e-6 ),
        c( 0,  0,  0,  0,  2,  0,  0,  0,   -63.53e-6,   0.02e-6 ),
        c( 0,  0,  2, -2,  3,  0,  0,  0,   -11.75e-6,  -0.01e-6 ),
        c( 0,  0,  2, -2,  1,  0,  0,  0,   -11.21e-6,  -0.01e-6 ),
        c( 0,  0,  2, -2,  2,  0,  0,  0,     4.57e-6,   0.00e-6 ),
        c( 0,  0,  2,  0,  3,  0,  0,  0,    -2.02e-6,   0.00e-6 ),
        c( 0,  0,  2,  0,  1,  0,  0,  0,    -1.98e-6,   0.00e-6 ),
        c( 0,  0,  0,  0,  3,  0,  0,  0,     1.72e-6,   0.00e-6 ),
        c( 0,  1,  0,  0,  1,  0,  0,  0,     1.41e-6,   0.01e-6 ),
        c( 0,  1,  0,  0, -1,  0,  0,  0,     1.26e-6,   0.01e-6 ),
        c( 1,  0,  0,  0, -1,  0,  0,  0,     0.63e-6,   0.00e-6 ),
        c( 1,  0,  0,  0,  1,  0,  0,  0,     0.63e-6,   0.00e-6 ),
        c( 0,  1,  2, -2,  3,  0,  0,  0,    -0.46e-6,   0.00e-6 ),
        c( 0,  1,  2, -2,  1,  0,  0,  0,    -0.45e-6,   0.00e-6 ),
        c( 0,  0,  4, -4,  4,  0,  0,  0,    -0.36e-6,   0.00e-6 ),
        c( 0,  0,  1, -1,  1, -8, 12,  0,     0.24e-6,   0.12e-6 ),
        c( 0,  0,  2,  0,  0,  0,  0,  0,    -0.32e-6,   0.00e-6 ),
        c( 0,  0,  2,  0,  2,  0,  0,  0,    -0.28e-6,   0.00e-6 ),
        c( 1,  0,  2,  0,  3,  0,  0,  0,    -0.27e-6,   0.00e-6 ),
        c( 1,  0,  2,  0,  1,  0,  0,  0,    -0.26e-6,   0.00e-6 ),
        c( 0,  0,  2, -2,  0,  0,  0,  0,     0.21e-6,   0.00e-6 ),
        c( 0,  1, -2,  2, -3,  0,  0,  0,    -0.19e-6,   0.00e-6 ),
        c( 0,  1, -2,  2, -1,  0,  0,  0,    -0.18e-6,   0.00e-6 ),
        c( 0,  0,  0,  0,  0,  8,-13, -1,     0.10e-6,  -0.05e-6 ),
        c( 0,  0,  0,  2,  0,  0,  0,  0,    -0.15e-6,   0.00e-6 ),
        c( 2,  0, -2,  0, -1,  0,  0,  0,     0.14e-6,   0.00e-6 ),
        c( 0,  1,  2, -2,  2,  0,  0,  0,     0.14e-6,   0.00e-6 ),
        c( 1,  0,  0, -2,  1,  0,  0,  0,    -0.14e-6,   0.00e-6 ),
        c( 1,  0,  0, -2, -1,  0,  0,  0,    -0.14e-6,   0.00e-6 ),
        c( 0,  0,  4, -2,  4,  0,  0,  0,    -0.13e-6,   0.00e-6 ),
        c( 0,  0,  2, -2,  4,  0,  0,  0,     0.11e-6,   0.00e-6 ),
        c( 1,  0, -2,  0, -3,  0,  0,  0,    -0.11e-6,   0.00e-6 ),
        c( 1,  0, -2,  0, -1,  0,  0,  0,    -0.11e-6,   0.00e-6 ))

    s1 <- rbind(c( 0,  0,  0,  0,  2,  0,  0,  0,    -0.07e-6,   3.57e-6 ),
                c( 0,  0,  0,  0,  1,  0,  0,  0, 1.73e-6,  -0.03e-6 ),
                c( 0,  0,  2, -2,  3,  0,  0,  0,     0.00e-6,   0.48e-6 ))

    s2 <- rbind(
        c( 0,  0,  0,  0,  1,  0,  0,  0,   743.52e-6,  -0.17e-6 ),
        c( 0,  0,  2, -2,  2,  0,  0,  0,    56.91e-6,   0.06e-6 ),
        c( 0,  0,  2,  0,  2,  0,  0,  0,     9.84e-6,  -0.01e-6 ),
        c( 0,  0,  0,  0,  2,  0,  0,  0,    -8.85e-6,   0.01e-6 ),
        c( 0,  1,  0,  0,  0,  0,  0,  0,    -6.38e-6,  -0.05e-6 ),
        c( 1,  0,  0,  0,  0,  0,  0,  0,    -3.07e-6,   0.00e-6 ),
        c( 0,  1,  2, -2,  2,  0,  0,  0,     2.23e-6,   0.00e-6 ),
        c( 0,  0,  2,  0,  1,  0,  0,  0,     1.67e-6,   0.00e-6 ),
        c( 1,  0,  2,  0,  2,  0,  0,  0,     1.30e-6,   0.00e-6 ),
        c( 0,  1, -2,  2, -2,  0,  0,  0,     0.93e-6,   0.00e-6 ),
        c( 1,  0,  0, -2,  0,  0,  0,  0,     0.68e-6,   0.00e-6 ),
        c( 0,  0,  2, -2,  1,  0,  0,  0,    -0.55e-6,   0.00e-6 ),
        c( 1,  0, -2,  0, -2,  0,  0,  0,     0.53e-6,   0.00e-6 ),
        c( 0,  0,  0,  2,  0,  0,  0,  0,    -0.27e-6,   0.00e-6 ),
        c( 1,  0,  0,  0,  1,  0,  0,  0,    -0.27e-6,   0.00e-6 ),
        c( 1,  0, -2, -2, -2,  0,  0,  0,    -0.26e-6,   0.00e-6 ),
        c( 1,  0,  0,  0, -1,  0,  0,  0,    -0.25e-6,   0.00e-6 ),
        c( 1,  0,  2,  0,  1,  0,  0,  0,     0.22e-6,   0.00e-6 ),
        c( 2,  0,  0, -2,  0,  0,  0,  0,    -0.21e-6,   0.00e-6 ),
        c( 2,  0, -2,  0, -1,  0,  0,  0,     0.20e-6,   0.00e-6 ),
        c( 0,  0,  2,  2,  2,  0,  0,  0,     0.17e-6,   0.00e-6 ),
        c( 2,  0,  2,  0,  2,  0,  0,  0,     0.13e-6,   0.00e-6 ),
        c( 2,  0,  0,  0,  0,  0,  0,  0,    -0.13e-6,   0.00e-6 ),
        c( 1,  0,  2, -2,  2,  0,  0,  0,    -0.12e-6,   0.00e-6 ),
        c( 0,  0,  2,  0,  0,  0,  0,  0,    -0.11e-6,   0.00e-6 ))

    s3 <- rbind(c( 0,  0,  0,  0,  1,  0,  0,  0,     0.30e-6, -23.42e-6 ),
                c( 0,  0,  2, -2,  2,  0,  0,  0,    -0.03e-6,  -1.46e-6 ),
                c( 0,  0,  2,  0,  2,  0,  0,  0,    -0.01e-6,  -0.25e-6 ),
                c( 0,  0,  0,  0,  2,  0,  0,  0,     0.00e-6,   0.23e-6 ))

    s4 <- rbind(c( 0,  0,  0,  0,  1,  0,  0,  0,    -0.26e-6,  -0.01e-6 ))

    NS0 <- nrow(s0)
    NS1 <- nrow(s1)
    NS2 <- nrow(s2)
    NS3 <- nrow(s3)
    NS4 <- nrow(s4)

    t <- ((tt[,1] - DJ00) + tt[,2]) / DJC
    ## Fundamental Arguments (from IERS Conventions 2003) */

    fa <- array(NA,dim=c(length(t),8))
    ## Mean anomaly of the Moon. */
    fa[,1] = sofa_Fal03(t)

    ## Mean anomaly of the Sun. */
    fa[,2] = sofa_Falp03(t)

    ## Mean longitude of the Moon minus that of the ascending node. */
    fa[,3] = sofa_Faf03(t)

    ## Mean elongation of the Moon from the Sun. */
    fa[,4] = sofa_Fad03(t)

    ## Mean longitude of the ascending node of the Moon. */
    fa[,5] = sofa_Faom03(t)

    ## Mean longitude of Venus. */
    fa[,6] = sofa_Fave03(t)

    ## Mean longitude of Earth. */
    fa[,7] = sofa_Fae03(t)

    ##General precession in longitude. */
    fa[,8] = sofa_Fapa03(t)

    ## Evaluate s. */
    w0 <- sp[1]
    w1 <- sp[2]
    w2 <- sp[3]
    w3 <- sp[4]
    w4 <- sp[5]
    w5 <- sp[6]

    for(i in NS0:1){
        a <- 0
        for(j in 1:8){
            a <- a+ s0[i,j]*fa[,j]
        }
        w0 <- w0+s0[i,9]*sin(a) + s0[i,10] * cos(a)
    }

    for (i in NS1:1){
        a <- 0.0
        for (j in 1:8) {
            a <- a + s1[i,j]*fa[,j]
        }
        w1 <- w1+s1[i,9] * sin(a) + s1[i,10]*cos(a)
    }

    for (i in NS2:1){
        a <- 0.0
        for(j in 1:8){
            a <- a+s2[i,j] * fa[,j]
        }
        w2 <- w2+s2[i,9]*sin(a) + s2[i,10]*cos(a)
    }

    for (i in NS3:1){
        a <- 0.0
        for(j in 1:8){
            a <- a+s3[i,j]*fa[,j]
        }
        w3 <- w3+s3[i,9]*sin(a)+s3[i,10] * cos(a);
    }

    for (i in NS4:1){
        a <- 0.0
        for(j in 1:8){
            a <- a+s4[i,j]*fa[,j]
        }
        w4 <- w4+s4[i,9]*sin(a) + s4[i,10]*cos(a)
    }
    s <- (w0 + (w1 + (w2 + (w3 + (w4 + w5 * t) * t) * t) * t) * t)*DAS2R - x*y/2.0
    return(s)
}

sofa_Pom00 <- function(xp, yp, sp){
####################################
## Construct the matrix W in https://www.iers.org/SharedDocs/Publikationen/EN/IERS/Publications/tn/TechnNote32/tn32.pdf?__blob=publicationFile&v=1
##  Form the matrix of polar motion for a given date, IAU 2000.
##
##  This function is part of the International Astronomical Union's
##  SOFA (Standards Of Fundamental Astronomy) software collection.
##
##  Status:  support function.
##
##  Given:
##     xp,yp    double    coordinates of the pole (radians, Note 1)
##     sp       double    the TIO locator s' (radians, Note 2)
##
##  Returned:
##     rpom     double[3][3]   polar-motion matrix (Note 3)
##
##  Notes:
##
##  1) The arguments xp and yp are the coordinates (in radians) of the
##     Celestial Intermediate Pole with respect to the International
##     Terrestrial Reference System (see IERS Conventions 2003),
##     measured along the meridians to 0 and 90 deg west respectively.
##
##  2) The argument sp is the TIO locator s', in radians, which
##     positions the Terrestrial Intermediate Origin on the equator.  It
##     is obtained from polar motion observations by numerical
##     integration, and so is in essence unpredictable.  However, it is
##     dominated by a secular drift of about 47 microarcseconds per
##     century, and so can be taken into account by using s' = -47*t,
##     where t is centuries since J2000.0.  The function iauSp00
##     implements this approximation.
##
##  3) The matrix operates in the sense V(TRS) = rpom * V(CIP), meaning
##     that it is the final rotation when computing the pointing
##     direction to a celestial source.
##
##  Called:
##     iauIr        initialize r-matrix to identity
##     iauRz        rotate around Z-axis
##     iauRy        rotate around Y-axis
##     iauRx        rotate around X-axis
##
##  Reference:
##
##     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
##     IERS Technical Note No. 32, BKG (2004)
##
##  This revision:  2013 June 18
##
##  SOFA release 2019-07-22
##
##  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
####################################
    rpom <- array(dim=c(3,3,length(xp)))
    for(j in 1:length(xp)){
        rpom[,,j] <- gen_R1(-yp[j])%*%gen_R2(-xp[j])%*%gen_R3(sp[j])
    }
    return(rpom)
}

sofa_Era00 <- function(ut){
####################################
##  Earth rotation angle (IAU 2000 model).
##
##  This function is part of the International Astronomical Union's
##  SOFA (Standards Of Fundamental Astronomy) software collection.
##
##  Status:  canonical model.
##
##  Given:
##     dj1,dj2   double    UT1 as a 2-part Julian Date (see note)
##
##  Returned (function value):
##               double    Earth rotation angle (radians), range 0-2pi
##
##  Notes:
##
##  1) The UT1 date dj1+dj2 is a Julian Date, apportioned in any
##     convenient way between the arguments dj1 and dj2.  For example,
##     JD(UT1)=2450123.7 could be expressed in any of these ways,
##     among others:
##
##             dj1            dj2
##
##         2450123.7           0.0       (JD method)
##         2451545.0       -1421.3       (J2000 method)
##         2400000.5       50123.2       (MJD method)
##         2450123.5           0.2       (date & time method)
##
##     The JD method is the most natural and convenient to use in
##     cases where the loss of several decimal digits of resolution
##     is acceptable.  The J2000 and MJD methods are good compromises
##     between resolution and convenience.  The date & time method is
##     best matched to the algorithm used:  maximum precision is
##     delivered when the dj1 argument is for 0hrs UT1 on the day in
##     question and the dj2 argument lies in the range 0 to 1, or vice
##     versa.
##
##  2) The algorithm is adapted from Expression 22 of Capitaine et al.
##     2000.  The time argument has been expressed in days directly,
##     and, to retain precision, integer contributions have been
##     eliminated.  The same formulation is given in IERS Conventions
##     (2003), Chap. 5, Eq. 14.
##
##  Called:
##     iauAnp       normalize angle into range 0 to 2pi
##
##  References:
##
##     Capitaine N., Guinot B. and McCarthy D.D, 2000, Astron.
##     Astrophys., 355, 398-405.
##
##     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
##     IERS Technical Note No. 32, BKG (2004)
##
##  This revision:  2013 June 18
##
##  SOFA release 2019-07-22
##
##  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
####################################

## Days since fundamental epoch.
   t = (ut[,1]-DJ00) + ut[,2]
## Fractional part of T (days).
   f = rowSums(ut)%%1

##Earth rotation angle at this UT1.
   theta <- 2*pi* ((f + 0.7790572732640 + 0.00273781191135448 * t)%%1)
   return(theta)
}

sofa_C2t00b <- function(tt, ut, xp, yp,rpom){
####################################
##  Form the celestial to terrestrial matrix given the date, the UT1 and
##  the polar motion, using the IAU 2000B nutation model.
##
##  This function is part of the International Astronomical Union's
##  SOFA (Standards Of Fundamental Astronomy) software collection.
##
##  Status:  support function.
##
##  Given:
##     tta,ttb  double         TT as a 2-part Julian Date (Note 1)
##     uta,utb  double         UT1 as a 2-part Julian Date (Note 1)
##     xp,yp    double         coordinates of the pole (radians, Note 2)
##
##  Returned:
##     rc2t     double[3][3]   celestial-to-terrestrial matrix (Note 3)
##
##  Notes:
##
##  1) The TT and UT1 dates tta+ttb and uta+utb are Julian Dates,
##     apportioned in any convenient way between the arguments uta and
##     utb.  For example, JD(UT1)=2450123.7 could be expressed in any of
##     these ways, among others:
##
##             uta            utb
##
##         2450123.7           0.0       (JD method)
##         2451545.0       -1421.3       (J2000 method)
##         2400000.5       50123.2       (MJD method)
##         2450123.5           0.2       (date & time method)
##
##     The JD method is the most natural and convenient to use in
##     cases where the loss of several decimal digits of resolution is
##     acceptable.  The J2000 and MJD methods are good compromises
##     between resolution and convenience.  In the case of uta,utb, the
##     date & time method is best matched to the Earth rotation angle
##     algorithm used:  maximum precision is delivered when the uta
##     argument is for 0hrs UT1 on the day in question and the utb
##     argument lies in the range 0 to 1, or vice versa.
##
##  2) The arguments xp and yp are the coordinates (in radians) of the
##     Celestial Intermediate Pole with respect to the International
##     Terrestrial Reference System (see IERS Conventions 2003),
##     measured along the meridians to 0 and 90 deg west respectively.
##
##  3) The matrix rc2t transforms from celestial to terrestrial
##     coordinates:
##
##        [TRS] = RPOM * R_3(ERA) * RC2I * [CRS]
##
##              = rc2t * [CRS]
##
##     where [CRS] is a vector in the Geocentric Celestial Reference
##     System and [TRS] is a vector in the International Terrestrial
##     Reference System (see IERS Conventions 2003), RC2I is the
##     celestial-to-intermediate matrix, ERA is the Earth rotation
##     angle and RPOM is the polar motion matrix.
##
##  4) The present function is faster, but slightly less accurate (about
##     1 mas), than the iauC2t00a function.
##
##  Called:
##     iauC2i00b    celestial-to-intermediate matrix, IAU 2000B
##     iauEra00     Earth rotation angle, IAU 2000
##     iauPom00     polar motion matrix
##     iauC2tcio    form CIO-based celestial-to-terrestrial matrix
##
##  Reference:
##
##     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
##     IERS Technical Note No. 32, BKG (2004)
##
##  This revision:  2013 June 18
##
##  SOFA release 2019-07-22
##
##  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
####################################
    ## Form the celestial-to-intermediate matrix for this TT (IAU 2000B).
    rc2i <- C2i00b(tt)
    ## Predict the Earth rotation angle for this UT1.
    era <- Era00(ut)
    ## Form the polar motion matrix (neglecting s').
#    rpom <- sofa_Pom00(xp, yp, rep(0,length(xp)))
    ## Combine to form the celestial-to-terrestrial matrix.
    rt2c <- rc2t <- array(NA,dim=c(3,3,nrow(tt)))
    for(j in 1:length(xp)){
        val <- rpom[,,j]%*%R3(era[j])%*%rc2i[,,j]
        rc2t[,,j] <- val
        rt2c[,,j] <- t(val)
    }
#    rc2t
    list(Mc2t=rc2t,Mt2c=rt2c)
#    return(list(Mcio=rc2i,Mc2t=rc2t,era=era,rpom=rpom))
}

##
sofa_C2i00a <- function(tt){
####################################
##  Form the celestial-to-intermediate matrix for a given date using the
##  IAU 2000A precession-nutation model.
##
##  This function is part of the International Astronomical Union's
##  SOFA (Standards Of Fundamental Astronomy) software collection.
##
##  Status:  support function.
##
##  Given:
##     date1,date2 double       TT as a 2-part Julian Date (Note 1)
##
##  Returned:
##     rc2i        double[3][3] celestial-to-intermediate matrix (Note 2)
##
##  Notes:
##
##  1) The TT date date1+date2 is a Julian Date, apportioned in any
##     convenient way between the two arguments.  For example,
##     JD(TT)=2450123.7 could be expressed in any of these ways,
##     among others:
##
##            date1          date2
##
##         2450123.7           0.0       (JD method)
##         2451545.0       -1421.3       (J2000 method)
##         2400000.5       50123.2       (MJD method)
##         2450123.5           0.2       (date & time method)
##
##     The JD method is the most natural and convenient to use in
##     cases where the loss of several decimal digits of resolution
##     is acceptable.  The J2000 method is best matched to the way
##     the argument is handled internally and will deliver the
##     optimum resolution.  The MJD method and the date & time methods
##     are both good compromises between resolution and convenience.
##
##  2) The matrix rc2i is the first stage in the transformation from
##     celestial to terrestrial coordinates:
##
##        [TRS]  =  RPOM * R_3(ERA) * rc2i * [CRS]
##
##               =  rc2t * [CRS]
##
##     where [CRS] is a vector in the Geocentric Celestial Reference
##     System and [TRS] is a vector in the International Terrestrial
##     Reference System (see IERS Conventions 2003), ERA is the Earth
##     Rotation Angle and RPOM is the polar motion matrix.
##
##  3) A faster, but slightly less accurate result (about 1 mas), can be
##     obtained by using instead the iauC2i00b function.
##
##  Called:
##     iauPnm00a    classical NPB matrix, IAU 2000A
##     iauC2ibpn    celestial-to-intermediate matrix, given NPB matrix
##
##  References:
##
##     "Expressions for the Celestial Intermediate Pole and Celestial
##     Ephemeris Origin consistent with the IAU 2000A precession-
##     nutation model", Astron.Astrophys. 400, 1145-1154
##     (2003)
##
##     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
##          intermediate origin" (CIO) by IAU 2006 Resolution 2.
##
##     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
##     IERS Technical Note No. 32, BKG (2004)
##
##  This revision:  2013 June 18
##
##  SOFA release 2019-07-22
##
##  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
####################################
    rbpn <- Pnm00a(tt)
    rc2i <- C2ibpn(tt, rbpn)
}

sofa_C2ibpn <- function(tt,rbpn){
####################################
##  Form the celestial-to-intermediate matrix for a given date given
##  the bias-precession-nutation matrix.  IAU 2000.
##
##  This function is part of the International Astronomical Union's
##  SOFA (Standards Of Fundamental Astronomy) software collection.
##
##  Status:  support function.
##
##  Given:
##     date1,date2 double       TT as a 2-part Julian Date (Note 1)
##     rbpn        double[3][3] celestial-to-true matrix (Note 2)
##
##  Returned:
##     rc2i        double[3][3] celestial-to-intermediate matrix (Note 3)
##
##  Notes:
##
##  1) The TT date date1+date2 is a Julian Date, apportioned in any
##     convenient way between the two arguments.  For example,
##     JD(TT)=2450123.7 could be expressed in any of these ways,
##     among others:
##
##            date1          date2
##
##         2450123.7           0.0       (JD method)
##         2451545.0       -1421.3       (J2000 method)
##         2400000.5       50123.2       (MJD method)
##         2450123.5           0.2       (date & time method)
##
##     The JD method is the most natural and convenient to use in
##     cases where the loss of several decimal digits of resolution
##     is acceptable.  The J2000 method is best matched to the way
##     the argument is handled internally and will deliver the
##     optimum resolution.  The MJD method and the date & time methods
##     are both good compromises between resolution and convenience.
##
##  2) The matrix rbpn transforms vectors from GCRS to true equator (and
##     CIO or equinox) of date.  Only the CIP (bottom row) is used.
##
##  3) The matrix rc2i is the first stage in the transformation from
##     celestial to terrestrial coordinates:
##
##        [TRS] = RPOM * R_3(ERA) * rc2i * [CRS]
##
##              = RC2T * [CRS]
##
##     where [CRS] is a vector in the Geocentric Celestial Reference
##     System and [TRS] is a vector in the International Terrestrial
##     Reference System (see IERS Conventions 2003), ERA is the Earth
##     Rotation Angle and RPOM is the polar motion matrix.
##
##  4) Although its name does not include "00", This function is in fact
##     specific to the IAU 2000 models.
##
##  Called:
##     iauBpn2xy    extract CIP X,Y coordinates from NPB matrix
##     iauC2ixy     celestial-to-intermediate matrix, given X,Y
##
##  References:
##     "Expressions for the Celestial Intermediate Pole and Celestial
##     Ephemeris Origin consistent with the IAU 2000A precession-
##     nutation model", Astron.Astrophys. 400, 1145-1154 (2003)
##
##     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
##          intermediate origin" (CIO) by IAU 2006 Resolution 2.
##
##     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
##     IERS Technical Note No. 32, BKG (2004)
##
##  This revision:  2013 June 18
##
##  SOFA release 2019-07-22
##
##  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
####################################
##Extract the X,Y coordinates.
    x <- rbpn[3,1,]
    y <- rbpn[3,2,]
##Form the celestial-to-intermediate matrix (n.b. IAU 2000 specific).
    rc2i <- C2ixy(tt, x, y)
    return(rc2i)
}

sofa_Pn00b <- function(tt){
####################################
##  Precession-nutation, IAU 2000B model:  a multi-purpose function,
##  supporting classical (equinox-based) use directly and CIO-based
##  use indirectly.
##
##  This function is part of the International Astronomical Union's
##  SOFA (Standards Of Fundamental Astronomy) software collection.
##
##  Status:  support function.
##
##  Given:
##     date1,date2  double          TT as a 2-part Julian Date (Note 1)
##
##  Returned:
##     dpsi,deps    double          nutation (Note 2)
##     epsa         double          mean obliquity (Note 3)
##     rb           double[3][3]    frame bias matrix (Note 4)
##     rp           double[3][3]    precession matrix (Note 5)
##     rbp          double[3][3]    bias-precession matrix (Note 6)
##     rn           double[3][3]    nutation matrix (Note 7)
##     rbpn         double[3][3]    GCRS-to-true matrix (Notes 8,9)
##
##  Notes:
##
##  1)  The TT date date1+date2 is a Julian Date, apportioned in any
##      convenient way between the two arguments.  For example,
##      JD(TT)=2450123.7 could be expressed in any of these ways,
##      among others:
##
##             date1          date2
##
##          2450123.7           0.0       (JD method)
##          2451545.0       -1421.3       (J2000 method)
##          2400000.5       50123.2       (MJD method)
##          2450123.5           0.2       (date & time method)
##
##      The JD method is the most natural and convenient to use in
##      cases where the loss of several decimal digits of resolution
##      is acceptable.  The J2000 method is best matched to the way
##      the argument is handled internally and will deliver the
##      optimum resolution.  The MJD method and the date & time methods
##      are both good compromises between resolution and convenience.
##
##  2)  The nutation components (luni-solar + planetary, IAU 2000B) in
##      longitude and obliquity are in radians and with respect to the
##      equinox and ecliptic of date.  For more accurate results, but
##      at the cost of increased computation, use the iauPn00a function.
##      For the utmost accuracy, use the iauPn00  function, where the
##      nutation components are caller-specified.
##
##  3)  The mean obliquity is consistent with the IAU 2000 precession.
##
##  4)  The matrix rb transforms vectors from GCRS to J2000.0 mean
##      equator and equinox by applying frame bias.
##
##  5)  The matrix rp transforms vectors from J2000.0 mean equator and
##      equinox to mean equator and equinox of date by applying
##      precession.
##
##  6)  The matrix rbp transforms vectors from GCRS to mean equator and
##      equinox of date by applying frame bias then precession.  It is
##      the product rp x rb.
##
##  7)  The matrix rn transforms vectors from mean equator and equinox
##      of date to true equator and equinox of date by applying the
##      nutation (luni-solar + planetary).
##
##  8)  The matrix rbpn transforms vectors from GCRS to true equator and
##      equinox of date.  It is the product rn x rbp, applying frame
##      bias, precession and nutation in that order.
##
##  9)  The X,Y,Z coordinates of the IAU 2000B Celestial Intermediate
##      Pole are elements (3,1-3) of the GCRS-to-true matrix,
##      i.e. rbpn[2][0-2].
##
##  10) It is permissible to re-use the same array in the returned
##      arguments.  The arrays are filled in the stated order.
##
##  Called:
##     iauNut00b    nutation, IAU 2000B
##     iauPn00      bias/precession/nutation results, IAU 2000
##
##  Reference:
##
##     Capitaine, N., Chapront, J., Lambert, S. and Wallace, P.,
##     "Expressions for the Celestial Intermediate Pole and Celestial
##     Ephemeris Origin consistent with the IAU 2000A precession-
##     nutation model", Astron.Astrophys. 400, 1145-1154 (2003).
##
##     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
##          intermediate origin" (CIO) by IAU 2006 Resolution 2.
##
##  This revision:  2013 November 13
##
##  SOFA release 2019-07-22
##
##  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
####################################

##Nutation.
    tmp1 <- Nut00b(tt)

##Remaining results.
    tmp2 <- Pn00(tt,tmp1$dpsi,tmp1$deps)
    return(c(tmp1,tmp2))
}

sofa_C2i00b <- function(tt){
####################################
##  Form the celestial-to-intermediate matrix for a given date using the
##  IAU 2000B precession-nutation model.
##
##  This function is part of the International Astronomical Union's
##  SOFA (Standards Of Fundamental Astronomy) software collection.
##
##  Status:  support function.
##
##  Given:
##     date1,date2 double       TT as a 2-part Julian Date (Note 1)
##
##  Returned:
##     rc2i        double[3][3] celestial-to-intermediate matrix (Note 2)
##
##  Notes:
##
##  1) The TT date date1+date2 is a Julian Date, apportioned in any
##     convenient way between the two arguments.  For example,
##     JD(TT)=2450123.7 could be expressed in any of these ways,
##     among others:
##
##            date1          date2
##
##         2450123.7           0.0       (JD method)
##         2451545.0       -1421.3       (J2000 method)
##         2400000.5       50123.2       (MJD method)
##         2450123.5           0.2       (date & time method)
##
##     The JD method is the most natural and convenient to use in
##     cases where the loss of several decimal digits of resolution
##     is acceptable.  The J2000 method is best matched to the way
##     the argument is handled internally and will deliver the
##     optimum resolution.  The MJD method and the date & time methods
##     are both good compromises between resolution and convenience.
##
##  2) The matrix rc2i is the first stage in the transformation from
##     celestial to terrestrial coordinates:
##
##        [TRS]  =  RPOM * R_3(ERA) * rc2i * [CRS]
##
##               =  rc2t * [CRS]
##
##     where [CRS] is a vector in the Geocentric Celestial Reference
##     System and [TRS] is a vector in the International Terrestrial
##     Reference System (see IERS Conventions 2003), ERA is the Earth
##     Rotation Angle and RPOM is the polar motion matrix.
##
##  3) The present function is faster, but slightly less accurate (about
##     1 mas), than the iauC2i00a function.
##
##  Called:
##     iauPnm00b    classical NPB matrix, IAU 2000B
##     iauC2ibpn    celestial-to-intermediate matrix, given NPB matrix
##
##  References:
##
##     "Expressions for the Celestial Intermediate Pole and Celestial
##     Ephemeris Origin consistent with the IAU 2000A precession-
##     nutation model", Astron.Astrophys. 400, 1145-1154
##     (2003)
##
##     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
##          intermediate origin" (CIO) by IAU 2006 Resolution 2.
##
##     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
##     IERS Technical Note No. 32, BKG (2004)
##
##  This revision:  2013 June 18
##
##  SOFA release 2019-07-22
##
##  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
####################################
##Obtain the celestial-to-true matrix (IAU 2000B).
#   rbpn <- Pnm00b(tt)
   rbpn <- Pn00b(tt)$rbpn

##Form the celestial-to-intermediate matrix.
   rc2i <- C2ibpn(tt, rbpn)

   return(rc2i)
}

sofa_Bi00 <- function(){
####################################
##  Frame bias components of IAU 2000 precession-nutation models (part
##  of MHB2000 with additions).
##
##  This function is part of the International Astronomical Union's
##  SOFA (Standards Of Fundamental Astronomy) software collection.
##
##  Status:  canonical model.
##
##  Returned:
##     dpsibi,depsbi  double  longitude and obliquity corrections
##     dra            double  the ICRS RA of the J2000.0 mean equinox
##
##  Notes:
##
##  1) The frame bias corrections in longitude and obliquity (radians)
##     are required in order to correct for the offset between the GCRS
##     pole and the mean J2000.0 pole.  They define, with respect to the
##     GCRS frame, a J2000.0 mean pole that is consistent with the rest
##     of the IAU 2000A precession-nutation model.
##
##  2) In addition to the displacement of the pole, the complete
##     description of the frame bias requires also an offset in right
##     ascension.  This is not part of the IAU 2000A model, and is from
##     Chapront et al. (2002).  It is returned in radians.
##
##  3) This is a supplemented implementation of one aspect of the IAU
##     2000A nutation model, formally adopted by the IAU General
##     Assembly in 2000, namely MHB2000 (Mathews et al. 2002).
##
##  References:
##
##     Chapront, J., Chapront-Touze, M. & Francou, G., Astron.
##     Astrophys., 387, 700, 2002.
##
##     Mathews, P.M., Herring, T.A., Buffet, B.A., "Modeling of nutation
##     and precession   New nutation series for nonrigid Earth and
##     insights into the Earth's interior", J.Geophys.Res., 107, B4,
##     2002.  The MHB2000 code itself was obtained on 9th September 2002
##     from ftp://maia.usno.navy.mil/conv2000/chapter5/IAU2000A.
##
##  This revision:  2013 June 18
##
##  SOFA release 2019-07-22
##
##  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
####################################
    DPBIAS = -0.041775  * DAS2R
    DEBIAS = -0.0068192 * DAS2R
    DRA0 = -0.0146 * DAS2R
    dpsibi = DPBIAS
    depsbi = DEBIAS
    dra = DRA0
    return(c(dpsibi, depsbi, dra))
}

sofa_Bp00 <- function(tt){
####################################
##  Frame bias and precession, IAU 2000.
##
##  This function is part of the International Astronomical Union's
##  SOFA (Standards Of Fundamental Astronomy) software collection.
##
##  Status:  canonical model.
##
##  Given:
##     date1,date2  double         TT as a 2-part Julian Date (Note 1)
##
##  Returned:
##     rb           double[3][3]   frame bias matrix (Note 2)
##     rp           double[3][3]   precession matrix (Note 3)
##     rbp          double[3][3]   bias-precession matrix (Note 4)
##
##  Notes:
##
##  1) The TT date date1+date2 is a Julian Date, apportioned in any
##     convenient way between the two arguments.  For example,
##     JD(TT)=2450123.7 could be expressed in any of these ways,
##     among others:
##
##             date1         date2
##
##         2450123.7           0.0       (JD method)
##         2451545.0       -1421.3       (J2000 method)
##         2400000.5       50123.2       (MJD method)
##         2450123.5           0.2       (date & time method)
##
##     The JD method is the most natural and convenient to use in
##     cases where the loss of several decimal digits of resolution
##     is acceptable.  The J2000 method is best matched to the way
##     the argument is handled internally and will deliver the
##     optimum resolution.  The MJD method and the date & time methods
##     are both good compromises between resolution and convenience.
##
##  2) The matrix rb transforms vectors from GCRS to mean J2000.0 by
##     applying frame bias.
##
##  3) The matrix rp transforms vectors from J2000.0 mean equator and
##     equinox to mean equator and equinox of date by applying
##     precession.
##
##  4) The matrix rbp transforms vectors from GCRS to mean equator and
##     equinox of date by applying frame bias then precession.  It is
##     the product rp x rb.
##
##  5) It is permissible to re-use the same array in the returned
##     arguments.  The arrays are filled in the order given.
##
##  Called:
##     iauBi00      frame bias components, IAU 2000
##     iauPr00      IAU 2000 precession adjustments
##     iauIr        initialize r-matrix to identity
##     iauRx        rotate around X-axis
##     iauRy        rotate around Y-axis
##     iauRz        rotate around Z-axis
##     iauCr        copy r-matrix
##     iauRxr       product of two r-matrices
##
##  Reference:
##     "Expressions for the Celestial Intermediate Pole and Celestial
##     Ephemeris Origin consistent with the IAU 2000A precession-
##     nutation model", Astron.Astrophys. 400, 1145-1154 (2003)
##
##     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
##          intermediate origin" (CIO) by IAU 2006 Resolution 2.
##
##  This revision:  2013 August 21
##
##  SOFA release 2019-07-22
##
##  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
####################################

##J2000.0 obliquity (Lieske et al. 1977)
   EPS0 = 84381.448 * DAS2R
   Nt <- nrow(tt)

##Interval between fundamental epoch J2000.0 and current date (JC).
   t = ((tt[,1] - DJ00) + tt[,2]) / DJC

##Frame bias.
   tmp <- Bi00()
   dpsibi <- tmp[1]
   depsbi <- tmp[2]
   dra0 <- tmp[3]

##Precession angles (Lieske et al. 1977)
   psia77 = (5038.7784 + (-1.07259 + (-0.001147) * t) * t) * t * DAS2R
   oma77  =       EPS0 + ((0.05127 + (-0.007726) * t) * t) * t * DAS2R
   chia   = (  10.5526 + (-2.38064 + (-0.001125) * t) * t) * t * DAS2R

##Apply IAU 2000 precession corrections.
   out <- Pr00(tt)
   psia = psia77 + out$dpsipr
   oma  = oma77  + out$depspr

   rb <- rp <- rbp <- array(NA,dim=c(3,3,Nt))
   for(j in 1:Nt){
       ##Frame bias matrix: GCRS to J2000.0.
       rbw <- R1(-depsbi)%*%R2(dpsibi*sin(EPS0))%*%R3(dra0)
       rb[,,j] <- rbw

       ##Precession matrix: J2000.0 to mean of date.
       rp[,,j] <- R3(chia[j])%*%R1(-oma[j])%*%R3(-psia[j])%*%R1(EPS0)

       ##Bias-precession matrix: GCRS to mean of date.
       rbp[,,j] <- rp[,,j]%*%rbw
   }
##return
   return(list(rb=rb,rp=rp,rbp=rbp))
}

sofa_Pr00 <- function(tt){
####################################
##  Precession-rate part of the IAU 2000 precession-nutation models
##  (part of MHB2000).
##
##  This function is part of the International Astronomical Union's
##  SOFA (Standards Of Fundamental Astronomy) software collection.
##
##  Status:  canonical model.
##
##  Given:
##     date1,date2    double  TT as a 2-part Julian Date (Note 1)
##
##  Returned:
##     dpsipr,depspr  double  precession corrections (Notes 2,3)
##
##  Notes:
##
##  1) The TT date date1+date2 is a Julian Date, apportioned in any
##     convenient way between the two arguments.  For example,
##     JD(TT)=2450123.7 could be expressed in any of these ways,
##     among others:
##
##            date1          date2
##
##         2450123.7           0.0       (JD method)
##         2451545.0       -1421.3       (J2000 method)
##         2400000.5       50123.2       (MJD method)
##         2450123.5           0.2       (date & time method)
##
##     The JD method is the most natural and convenient to use in
##     cases where the loss of several decimal digits of resolution
##     is acceptable.  The J2000 method is best matched to the way
##     the argument is handled internally and will deliver the
##     optimum resolution.  The MJD method and the date & time methods
##     are both good compromises between resolution and convenience.
##
##  2) The precession adjustments are expressed as "nutation
##     components", corrections in longitude and obliquity with respect
##     to the J2000.0 equinox and ecliptic.
##
##  3) Although the precession adjustments are stated to be with respect
##     to Lieske et al. (1977), the MHB2000 model does not specify which
##     set of Euler angles are to be used and how the adjustments are to
##     be applied.  The most literal and straightforward procedure is to
##     adopt the 4-rotation epsilon_0, psi_A, omega_A, xi_A option, and
##     to add dpsipr to psi_A and depspr to both omega_A and eps_A.
##
##  4) This is an implementation of one aspect of the IAU 2000A nutation
##     model, formally adopted by the IAU General Assembly in 2000,
##     namely MHB2000 (Mathews et al. 2002).
##
##  References:
##
##     Lieske, J.H., Lederle, T., Fricke, W. & Morando, B., "Expressions
##     for the precession quantities based upon the IAU (1976) System of
##     Astronomical Constants", Astron.Astrophys., 58, 1-16 (1977)
##
##     Mathews, P.M., Herring, T.A., Buffet, B.A., "Modeling of nutation
##     and precession   New nutation series for nonrigid Earth and
##     insights into the Earth's interior", J.Geophys.Res., 107, B4,
##     2002.  The MHB2000 code itself was obtained on 9th September 2002
##     from ftp://maia.usno.navy.mil/conv2000/chapter5/IAU2000A.
##
##     Wallace, P.T., "Software for Implementing the IAU 2000
##     Resolutions", in IERS Workshop 5.1 (2002).
##
##  This revision:  2013 June 18
##
##  SOFA release 2019-07-22
##
##  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
####################################
    ##Precession and obliquity corrections (radians per century)
    PRECOR = -0.29965 * DAS2R
    OBLCOR = -0.02524 * DAS2R

    ##Interval between fundamental epoch J2000.0 and given date (JC).
    t = ((tt[,1] - DJ00) + tt[,2]) / DJC

    ##Precession rate contributions with respect to IAU 1976/80.
    dpsipr <- PRECOR * t
    depspr <- OBLCOR * t
    return(list(dpsipr=dpsipr,depspr=depspr))
}

sofa_Pn00 <- function(tt,dpsi,deps){
####################################
##  Precession-nutation, IAU 2000 model:  a multi-purpose function,
##  supporting classical (equinox-based) use directly and CIO-based
##  use indirectly.
##
##  This function is part of the International Astronomical Union's
##  SOFA (Standards Of Fundamental Astronomy) software collection.
##
##  Status:  support function.
##
##  Given:
##     date1,date2  double          TT as a 2-part Julian Date (Note 1)
##     dpsi,deps    double          nutation (Note 2)
##
##  Returned:
##     epsa         double          mean obliquity (Note 3)
##     rb           double[3][3]    frame bias matrix (Note 4)
##     rp           double[3][3]    precession matrix (Note 5)
##     rbp          double[3][3]    bias-precession matrix (Note 6)
##     rn           double[3][3]    nutation matrix (Note 7)
##     rbpn         double[3][3]    GCRS-to-true matrix (Note 8)
##
##  Notes:
##
##  1) The TT date date1+date2 is a Julian Date, apportioned in any
##     convenient way between the two arguments.  For example,
##     JD(TT)=2450123.7 could be expressed in any of these ways,
##     among others:
##
##            date1          date2
##
##         2450123.7           0.0       (JD method)
##         2451545.0       -1421.3       (J2000 method)
##         2400000.5       50123.2       (MJD method)
##         2450123.5           0.2       (date & time method)
##
##     The JD method is the most natural and convenient to use in
##     cases where the loss of several decimal digits of resolution
##     is acceptable.  The J2000 method is best matched to the way
##     the argument is handled internally and will deliver the
##     optimum resolution.  The MJD method and the date & time methods
##     are both good compromises between resolution and convenience.
##
##  2) The caller is responsible for providing the nutation components;
##     they are in longitude and obliquity, in radians and are with
##     respect to the equinox and ecliptic of date.  For high-accuracy
##     applications, free core nutation should be included as well as
##     any other relevant corrections to the position of the CIP.
##
##  3) The returned mean obliquity is consistent with the IAU 2000
##     precession-nutation models.
##
##  4) The matrix rb transforms vectors from GCRS to J2000.0 mean
##     equator and equinox by applying frame bias.
##
##  5) The matrix rp transforms vectors from J2000.0 mean equator and
##     equinox to mean equator and equinox of date by applying
##     precession.
##
##  6) The matrix rbp transforms vectors from GCRS to mean equator and
##     equinox of date by applying frame bias then precession.  It is
##     the product rp x rb.
##
##  7) The matrix rn transforms vectors from mean equator and equinox of
##     date to true equator and equinox of date by applying the nutation
##     (luni-solar + planetary).
##
##  8) The matrix rbpn transforms vectors from GCRS to true equator and
##     equinox of date.  It is the product rn x rbp, applying frame
##     bias, precession and nutation in that order.
##
##  9) It is permissible to re-use the same array in the returned
##     arguments.  The arrays are filled in the order given.
##
##  Called:
##     iauPr00      IAU 2000 precession adjustments
##     iauObl80     mean obliquity, IAU 1980
##     iauBp00      frame bias and precession matrices, IAU 2000
##     iauCr        copy r-matrix
##     iauNumat     form nutation matrix
##     iauRxr       product of two r-matrices
##
##  Reference:
##
##     Capitaine, N., Chapront, J., Lambert, S. and Wallace, P.,
##     "Expressions for the Celestial Intermediate Pole and Celestial
##     Ephemeris Origin consistent with the IAU 2000A precession-
##     nutation model", Astron.Astrophys. 400, 1145-1154 (2003)
##
##     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
##          intermediate origin" (CIO) by IAU 2006 Resolution 2.
##
##  This revision:  2013 June 18
##
##  SOFA release 2019-07-22
##
##  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
####################################
    Nt <- nrow(tt)
    ## IAU 2000 precession-rate adjustments. , &dpsipr, &depspr)
    tmp <- Pr00(tt)
    dpsipr <- tmp$dpsipr
    depspr <- tmp$depspr

    ##Mean obliquity, consistent with IAU 2000 precession-nutation.
    epsa <- sofa_Obl80(tt) + depspr

    ##Frame bias and precession matrices and their product.
    rbpw <- Bp00(tt)$rbp
    rb <- Bp00(tt)$rb
    rp <- Bp00(tt)$rp
    rbp <- rbpw

    ##Nutation matrix.
    rn <- rbpn <- array(NA,dim=c(3,3,Nt))
    for(j in 1:Nt){
        rnw <- R1(-(epsa[j]+deps[j]))%*%R3(-dpsi[j])%*%R1(epsa[j])
        rn[,,j] <- rnw
        ##Bias-precession-nutation matrix (classical).
        rbpn[,,j] <- rnw%*%rbpw[,,j]
    }
    return(list(rb=rb,rp=rp,rbp=rbp,rn=rn,rbpn=rbpn))
}

sofa_Obl80 <- function(tt){
####################################
##  Mean obliquity of the ecliptic, IAU 1980 model.
##
##  This function is part of the International Astronomical Union's
##  SOFA (Standards Of Fundamental Astronomy) software collection.
##
##  Status:  canonical model.
##
##  Given:
##     date1,date2   double    TT as a 2-part Julian Date (Note 1)
##
##  Returned (function value):
##                   double    obliquity of the ecliptic (radians, Note 2)
##
##  Notes:
##
##  1) The TT date date1+date2 is a Julian Date, apportioned in any
##     convenient way between the two arguments.  For example,
##     JD(TT)=2450123.7 could be expressed in any of these ways,
##     among others:
##
##            date1          date2
##
##         2450123.7           0.0       (JD method)
##         2451545.0       -1421.3       (J2000 method)
##         2400000.5       50123.2       (MJD method)
##         2450123.5           0.2       (date & time method)
##
##     The JD method is the most natural and convenient to use in
##     cases where the loss of several decimal digits of resolution
##     is acceptable.  The J2000 method is best matched to the way
##     the argument is handled internally and will deliver the
##     optimum resolution.  The MJD method and the date & time methods
##     are both good compromises between resolution and convenience.
##
##  2) The result is the angle between the ecliptic and mean equator of
##     date date1+date2.
##
##  Reference:
##
##     Explanatory Supplement to the Astronomical Almanac,
##     P. Kenneth Seidelmann (ed), University Science Books (1992),
##     Expression 3.222-1 (p114).
##
##  This revision:  2013 June 18
##
##  SOFA release 2019-07-22
##
##  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
####################################
##Interval between fundamental epoch J2000.0 and given date (JC).
   t <- ((tt[,1] - DJ00) + tt[,2])/DJC

##Mean obliquity of date.
   eps0 <- DAS2R * (84381.448  +
                  (-46.8150   +
                  (-0.00059   +
                  ( 0.001813) * t) * t) * t)
   return(eps0)
}

sofa_Pnm06a <- function(tt){
####################################
##  Form the matrix of precession-nutation for a given date (including
##  frame bias), IAU 2006 precession and IAU 2000A nutation models.
##
##  This function is part of the International Astronomical Union's
##  SOFA (Standards Of Fundamental Astronomy) software collection.
##
##  Status:  support function.
##
##  Given:
##     date1,date2 double       TT as a 2-part Julian Date (Note 1)
##
##  Returned:
##     rnpb        double[3][3] bias-precession-nutation matrix (Note 2)
##
##  Notes:
##
##  1) The TT date date1+date2 is a Julian Date, apportioned in any
##     convenient way between the two arguments.  For example,
##     JD(TT)=2450123.7 could be expressed in any of these ways,
##     among others:
##
##            date1          date2
##
##         2450123.7           0.0       (JD method)
##         2451545.0       -1421.3       (J2000 method)
##         2400000.5       50123.2       (MJD method)
##         2450123.5           0.2       (date & time method)
##
##     The JD method is the most natural and convenient to use in
##     cases where the loss of several decimal digits of resolution
##     is acceptable.  The J2000 method is best matched to the way
##     the argument is handled internally and will deliver the
##     optimum resolution.  The MJD method and the date & time methods
##     are both good compromises between resolution and convenience.
##
##  2) The matrix operates in the sense V(date) = rnpb * V(GCRS), where
##     the p-vector V(date) is with respect to the true equatorial triad
##     of date date1+date2 and the p-vector V(GCRS) is with respect to
##     the Geocentric Celestial Reference System (IAU, 2000).
##
##  Called:
##     iauPfw06     bias-precession F-W angles, IAU 2006
##     iauNut06a    nutation, IAU 2006/2000A
##     iauFw2m      F-W angles to r-matrix
##
##  Reference:
##
##     Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855.
##
##  This revision:  2013 June 18
##
##  SOFA release 2019-07-22
##
##  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
####################################
##Fukushima-Williams angles for frame bias and precession.
   gppe <- Pfw06(tt)

##Nutation components.
   dpe <- Nut06a(tt)

##Equinox based nutation x precession x bias matrix. return rnpb
   R1(-gppe[4])%*%R3(-(gppe[3] + dpe[1]))%*%R1(gppe[2])%*%R3(gppe[1])
}

sofa_Pfw06 <- function(tt){
####################################
##  Precession angles, IAU 2006 (Fukushima-Williams 4-angle formulation).
##
##  This function is part of the International Astronomical Union's
##  SOFA (Standards Of Fundamental Astronomy) software collection.
##
##  Status:  canonical model.
##
##  Given:
##     date1,date2  double   TT as a 2-part Julian Date (Note 1)
##
##  Returned:
##     gamb         double   F-W angle gamma_bar (radians)
##     phib         double   F-W angle phi_bar (radians)
##     psib         double   F-W angle psi_bar (radians)
##     epsa         double   F-W angle epsilon_A (radians)
##
##  Notes:
##
##  1) The TT date date1+date2 is a Julian Date, apportioned in any
##     convenient way between the two arguments.  For example,
##     JD(TT)=2450123.7 could be expressed in any of these ways,
##     among others:
##
##            date1          date2
##
##         2450123.7           0.0       (JD method)
##         2451545.0       -1421.3       (J2000 method)
##         2400000.5       50123.2       (MJD method)
##         2450123.5           0.2       (date & time method)
##
##     The JD method is the most natural and convenient to use in
##     cases where the loss of several decimal digits of resolution
##     is acceptable.  The J2000 method is best matched to the way
##     the argument is handled internally and will deliver the
##     optimum resolution.  The MJD method and the date & time methods
##     are both good compromises between resolution and convenience.
##
##  2) Naming the following points:
##
##           e = J2000.0 ecliptic pole,
##           p = GCRS pole,
##           E = mean ecliptic pole of date,
##     and   P = mean pole of date,
##
##     the four Fukushima-Williams angles are as follows:
##
##        gamb = gamma_bar = epE
##        phib = phi_bar = pE
##        psib = psi_bar = pEP
##        epsa = epsilon_A = EP
##
##  3) The matrix representing the combined effects of frame bias and
##     precession is:
##
##        PxB = R_1(-epsa).R_3(-psib).R_1(phib).R_3(gamb)
##
##  4) The matrix representing the combined effects of frame bias,
##     precession and nutation is simply:
##
##        NxPxB = R_1(-epsa-dE).R_3(-psib-dP).R_1(phib).R_3(gamb)
##
##     where dP and dE are the nutation components with respect to the
##     ecliptic of date.
##
##  Reference:
##
##     Hilton, J. et al., 2006, Celest.Mech.Dyn.Astron. 94, 351
##
##  Called:
##     iauObl06     mean obliquity, IAU 2006
##
##  This revision:  2013 June 18
##
##  SOFA release 2019-07-22
##
##  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
####################################
    t <- ((tt[1] - DJ00) + tt[2]) / DJC

    ##P03 bias+precession angles.
    gamb <- (    -0.052928     +
                 (    10.556378     +
                      (     0.4932044    +
                           (    -0.00031238   +
                                (    -0.000002788  +
                                     (     0.0000000260 )
                                 * t) * t) * t) * t) * t) * DAS2R
    phib <- ( 84381.412819     +
                 (   -46.811016     +
                      (     0.0511268    +
                           (     0.00053289   +
                                (    -0.000000440  +
                                     (    -0.0000000176 )
                                 * t) * t) * t) * t) * t) * DAS2R
    psib <- (    -0.041775     +
                 (  5038.481484     +
                      (     1.5584175    +
                           (    -0.00018522   +
                                (    -0.000026452  +
                                     (    -0.0000000148 )
                                 * t) * t) * t) * t) * t) * DAS2R
    epsa <- iauObl06(date1, date2)
    return(c(gamb,phib,psib,epsa))
}

##calculate classical NPB matrix
sofa_Pnm00a <- function(tt){
####################################
##  Form the matrix of precession-nutation for a given date (including
##  frame bias), equinox-based, IAU 2000A model.
##
##  This function is part of the International Astronomical Union's
##  SOFA (Standards Of Fundamental Astronomy) software collection.
##
##  Status:  support function.
##
##  Given:
##     date1,date2  double     TT as a 2-part Julian Date (Note 1)
##
##  Returned:
##     rbpn         double[3][3]    classical NPB matrix (Note 2)
##
##  Notes:
##
##  1) The TT date date1+date2 is a Julian Date, apportioned in any
##     convenient way between the two arguments.  For example,
##     JD(TT)=2450123.7 could be expressed in any of these ways,
##     among others:
##
##            date1          date2
##
##         2450123.7           0.0       (JD method)
##         2451545.0       -1421.3       (J2000 method)
##         2400000.5       50123.2       (MJD method)
##         2450123.5           0.2       (date & time method)
##
##     The JD method is the most natural and convenient to use in
##     cases where the loss of several decimal digits of resolution
##     is acceptable.  The J2000 method is best matched to the way
##     the argument is handled internally and will deliver the
##     optimum resolution.  The MJD method and the date & time methods
##     are both good compromises between resolution and convenience.
##
##  2) The matrix operates in the sense V(date) = rbpn * V(GCRS), where
##     the p-vector V(date) is with respect to the true equatorial triad
##     of date date1+date2 and the p-vector V(GCRS) is with respect to
##     the Geocentric Celestial Reference System (IAU, 2000).
##
##  3) A faster, but slightly less accurate result (about 1 mas), can be
##     obtained by using instead the iauPnm00b function.
##
##  Called:
##     iauPn00a     bias/precession/nutation, IAU 2000A
##
##  Reference:
##
##     IAU: Trans. International Astronomical Union, Vol. XXIVB;  Proc.
##     24th General Assembly, Manchester, UK.  Resolutions B1.3, B1.6.
##     (2000)
##
##  This revision:  2013 June 18
##
##  SOFA release 2019-07-22
##
##  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
####################################
    tmp <- Nut00a(tt)/DAS2R
    Pn00a(tt)
}

#Nut06a
sofa_Nut06a <- function(tt){
####################################
##  IAU 2000A nutation with adjustments to match the IAU 2006
##  precession.
##
##  Given:
##     date1,date2   double   TT as a 2-part Julian Date (Note 1)
##
##  Returned:
##     dpsi,deps     double   nutation, luni-solar + planetary (Note 2)
##
##  Status:  canonical model.
##
##  Notes:
##
##  1) The TT date date1+date2 is a Julian Date, apportioned in any
##     convenient way between the two arguments.  For example,
##     JD(TT)=2450123.7 could be expressed in any of these ways,
##     among others:
##
##            date1          date2
##
##         2450123.7           0.0       (JD method)
##         2451545.0       -1421.3       (J2000 method)
##         2400000.5       50123.2       (MJD method)
##         2450123.5           0.2       (date & time method)
##
##     The JD method is the most natural and convenient to use in
##     cases where the loss of several decimal digits of resolution
##     is acceptable.  The J2000 method is best matched to the way
##     the argument is handled internally and will deliver the
##     optimum resolution.  The MJD method and the date & time methods
##     are both good compromises between resolution and convenience.
##
##  2) The nutation components in longitude and obliquity are in radians
##     and with respect to the mean equinox and ecliptic of date,
##     IAU 2006 precession model (Hilton et al. 2006, Capitaine et al.
##     2005).
##
##  3) The function first computes the IAU 2000A nutation, then applies
##     adjustments for (i) the consequences of the change in obliquity
##     from the IAU 1980 ecliptic to the IAU 2006 ecliptic and (ii) the
##     secular variation in the Earth's dynamical form factor J2.
##
##  4) The present function provides classical nutation, complementing
##     the IAU 2000 frame bias and IAU 2006 precession.  It delivers a
##     pole which is at current epochs accurate to a few tens of
##     microarcseconds, apart from the free core nutation.
##
##  Called:
##     iauNut00a    nutation, IAU 2000A
##
##  References:
##
##     Chapront, J., Chapront-Touze, M. & Francou, G. 2002,
##     Astron.Astrophys. 387, 700
##
##     Lieske, J.H., Lederle, T., Fricke, W. & Morando, B. 1977,
##     Astron.Astrophys. 58, 1-16
##
##     Mathews, P.M., Herring, T.A., Buffet, B.A. 2002, J.Geophys.Res.
##     107, B4.  The MHB_2000 code itself was obtained on 9th September
##     2002 from ftp//maia.usno.navy.mil/conv2000/chapter5/IAU2000A.
##
##     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
##     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
##
##     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
##     Astron.Astrophys.Supp.Ser. 135, 111
##
##     Wallace, P.T., "Software for Implementing the IAU 2000
##     Resolutions", in IERS Workshop 5.1 (2002)
##
##  This revision:  2013 June 18
##
##  SOFA release 2019-07-22
##
##  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
####################################
    date1 <- tt[,1]
    date2 <- tt[,2]
    ##Interval between fundamental date J2000.0 and given date (JC).
    t = ((date1 - DJ00) + date2) / DJC;

    ##Factor correcting for secular variation of J2.
    fj2 = -2.7774e-6 * t;

    ##Obtain IAU 2000A nutation.
    dpe <- iauNut00a(tt)

    ##Apply P03 adjustments (Wallace & Capitaine, 2006, Eqs.5).
    dpsi <- dp + dp * (0.4697e-6 + fj2)
    deps <- de + de * fj2
    return(c(epsi,deps))
}

sofa_Gc2gd <- function(n,xyz){
####################################
##  Transform geocentric coordinates to geodetic using the specified
##  reference ellipsoid.
##
##  This function is part of the International Astronomical Union's
##  SOFA (Standards of Fundamental Astronomy) software collection.
##
##  Status:  canonical transformation.
##
##  Given:
##     n       int        ellipsoid identifier (Note 1)
##     xyz     double[3]  geocentric vector (Note 2)
##
##  Returned:
##     elong   double     longitude (radians, east +ve, Note 3)
##     phi     double     latitude (geodetic, radians, Note 3)
##     height  double     height above ellipsoid (geodetic, Notes 2,3)
##
##  Returned (function value):
##            int         status:  0 = OK
##                                -1 = illegal identifier (Note 3)
##                                -2 = internal error (Note 3)
##
##  Notes:
##
##  1) The identifier n is a number that specifies the choice of
##     reference ellipsoid.  The following are supported:
##
##        n    ellipsoid
##
##        1     WGS84
##        2     GRS80
##        3     WGS72
##
##     The n value has no significance outside the SOFA software.  For
##     convenience, symbols WGS84 etc. are defined in sofam.h.
##
##  2) The geocentric vector (xyz, given) and height (height, returned)
##     are in meters.
##
##  3) An error status -1 means that the identifier n is illegal.  An
##     error status -2 is theoretically impossible.  In all error cases,
##     all three results are set to -1e9.
##
##  4) The inverse transformation is performed in the function iauGd2gc.
##
##  Called:
##     iauEform     Earth reference ellipsoids
##     iauGc2gde    geocentric to geodetic transformation, general
##
##  This revision:  2013 September 1
##
##  SOFA release 2019-07-22
##
##  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
####################################

    af <- sofa_Eform( n )
    sofa_Gc2gde( af[1], af[2], xyz)
}


sofa_Gc2gde <- function(a,f,xyz){
####################################
##  Transform geocentric coordinates to geodetic for a reference
##  ellipsoid of specified form.
##
##  This function is part of the International Astronomical Union's
##  SOFA (Standards of Fundamental Astronomy) software collection.
##
##  Status:  support function.
##
##  Given:
##     a       double     equatorial radius (Notes 2,4)
##     f       double     flattening (Note 3)
##     xyz     double[3]  geocentric vector (Note 4)
##
##  Returned:
##     elong   double     longitude (radians, east +ve)
##     phi     double     latitude (geodetic, radians)
##     height  double     height above ellipsoid (geodetic, Note 4)
##
##  Returned (function value):
##             int        status:  0 = OK
##                                -1 = illegal f
##                                -2 = illegal a
##
##  Notes:
##
##  1) This function is based on the GCONV2H Fortran subroutine by
##     Toshio Fukushima (see reference).
##
##  2) The equatorial radius, a, can be in any units, but meters is
##     the conventional choice.
##
##  3) The flattening, f, is (for the Earth) a value around 0.00335,
##     i.e. around 1/298.
##
##  4) The equatorial radius, a, and the geocentric vector, xyz,
##     must be given in the same units, and determine the units of
##     the returned height, height.
##
##  5) If an error occurs (status < 0), elong, phi and height are
##     unchanged.
##
##  6) The inverse transformation is performed in the function
##     iauGd2gce.
##
##  7) The transformation for a standard ellipsoid (such as WGS84) can
##     more conveniently be performed by calling iauGc2gd, which uses a
##     numerical code to identify the required A and F values.
##
##  Reference:
##
##     Fukushima, T., "Transformation from Cartesian to geodetic
##     coordinates accelerated by Halley's method", J.Geodesy (2006)
##     79: 689-693
##
##  This revision:  2014 November 7
##
##  SOFA release 2019-07-22
##
##  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
####################################

##Functions of ellipsoid parameters (with further validation of f).
   aeps2 = a*a * 1e-32
   e2 = (2.0 - f) * f
   e4t = e2*e2 * 1.5
   ec2 = 1.0 - e2
   if ( ec2 <= 0.0 ) return(-1)
   ec = sqrt(ec2)
   b = a * ec

##Cartesian components.
   x = xyz[1]
   y = xyz[2]
   z = xyz[3]

##Distance from polar axis squared.
   p2 <- x*x + y*y

##Longitude.
   elong <- atan2(y, x)

##Unsigned z-coordinate.
   absz <- abs(z)
##Proceed unless polar case.
   if ( p2 > aeps2 ){
##Distance from polar axis.
      p <- sqrt(p2)
##Normalization.
      s0 = absz/a
      pn = p/a
      zc = ec*s0

##Prepare Newton correction factors.
      c0 = ec * pn
      c02 = c0 * c0
      c03 = c02 * c0
      s02 = s0 * s0
      s03 = s02 * s0
      a02 = c02 + s02
      a0 = sqrt(a02)
      a03 = a02 * a0
      d0 = zc*a03 + e2*s03
      f0 = pn*a03 - e2*c03

##Prepare Halley correction factor.
      b0 <- e4t * s02 * c02 * pn * (a0 - ec)
      s1 <- d0*f0 - b0*s0
      cc <- ec * (f0*f0 - b0*c0)

##Evaluate latitude and height.
      phi <- atan(s1/cc)
      s12 <- s1 * s1
      cc2 <- cc * cc
      height <- (p*cc + absz*s1 - a * sqrt(ec2*s12 + cc2)) /sqrt(s12 + cc2)
   }else{
##Exception: pole
      phi <- DPI / 2.0
      height <- absz - b
   }

##Restore sign of latitude.
   if ( z < 0 ) phi <- -phi
   return(c(elong=elong,phi=phi,height=height))
}

sofa_C2ixy <- function(tt,x,y){
####################################
##  Form the celestial to intermediate-frame-of-date matrix for a given
##  date when the CIP X,Y coordinates are known.  IAU 2000.
##
##  This function is part of the International Astronomical Union's
##  SOFA (Standards Of Fundamental Astronomy) software collection.
##
##  Status:  support function.
##
##  Given:
##     date1,date2 double       TT as a 2-part Julian Date (Note 1)
##     x,y         double       Celestial Intermediate Pole (Note 2)
##
##  Returned:
##     rc2i        double[3][3] celestial-to-intermediate matrix (Note 3)
##
##  Notes:
##
##  1) The TT date date1+date2 is a Julian Date, apportioned in any
##     convenient way between the two arguments.  For example,
##     JD(TT)=2450123.7 could be expressed in any of these ways,
##     among others:
##
##            date1          date2
##
##         2450123.7           0.0       (JD method)
##         2451545.0       -1421.3       (J2000 method)
##         2400000.5       50123.2       (MJD method)
##         2450123.5           0.2       (date & time method)
##
##     The JD method is the most natural and convenient to use in
##     cases where the loss of several decimal digits of resolution
##     is acceptable.  The J2000 method is best matched to the way
##     the argument is handled internally and will deliver the
##     optimum resolution.  The MJD method and the date & time methods
##     are both good compromises between resolution and convenience.
##
##  2) The Celestial Intermediate Pole coordinates are the x,y components
##     of the unit vector in the Geocentric Celestial Reference System.
##
##  3) The matrix rc2i is the first stage in the transformation from
##     celestial to terrestrial coordinates:
##
##        [TRS] = RPOM * R_3(ERA) * rc2i * [CRS]
##
##              = RC2T * [CRS]
##
##     where [CRS] is a vector in the Geocentric Celestial Reference
##     System and [TRS] is a vector in the International Terrestrial
##     Reference System (see IERS Conventions 2003), ERA is the Earth
##     Rotation Angle and RPOM is the polar motion matrix.
##
##  4) Although its name does not include "00", This function is in fact
##     specific to the IAU 2000 models.
##
##  Called:
##     iauC2ixys    celestial-to-intermediate matrix, given X,Y and s
##     iauS00       the CIO locator s, given X,Y, IAU 2000A
##
##  Reference:
##
##     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
##     IERS Technical Note No. 32, BKG (2004)
##
##  This revision:  2013 June 18
##
##  SOFA release 2019-07-22
##
##  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
####################################
    Nt <- nrow(tt)
    rc2i <- array(NA,c(3,3,Nt))
    s <- S06(tt, x, y)
    ## Obtain the spherical angles E and d. */
    r2 <- x*x + y*y
    e <- atan2(y, x)
    d <- atan(sqrt(r2/(1.0 - r2)))

    ## Form the matrix.
    for(j in 1:Nt){
        rc2i[,,j] <- R3(-(e[j]+s[j]))%*%R2(d[j])%*%R3(e[j])
    }
    return(rc2i)
}

sofa_Nut00a <- function(tt){
####################################
##  Nutation, IAU 2000A model (MHB2000 luni-solar and planetary nutation
##  with free core nutation omitted).
##
##  This function is part of the International Astronomical Union's
##  SOFA (Standards Of Fundamental Astronomy) software collection.
##
##  Status:  canonical model.
##
##  Given:
##     date1,date2   double   TT as a 2-part Julian Date (Note 1)
##
##  Returned:
##     dpsi,deps     double   nutation, luni-solar + planetary (Note 2)
##
##  Notes:
##
##  1) The TT date date1+date2 is a Julian Date, apportioned in any
##     convenient way between the two arguments.  For example,
##     JD(TT)=2450123.7 could be expressed in any of these ways,
##     among others:
##
##            date1          date2
##
##         2450123.7           0.0       (JD method)
##         2451545.0       -1421.3       (J2000 method)
##         2400000.5       50123.2       (MJD method)
##         2450123.5           0.2       (date & time method)
##
##     The JD method is the most natural and convenient to use in
##     cases where the loss of several decimal digits of resolution
##     is acceptable.  The J2000 method is best matched to the way
##     the argument is handled internally and will deliver the
##     optimum resolution.  The MJD method and the date & time methods
##     are both good compromises between resolution and convenience.
##
##  2) The nutation components in longitude and obliquity are in radians
##     and with respect to the equinox and ecliptic of date.  The
##     obliquity at J2000.0 is assumed to be the Lieske et al. (1977)
##     value of 84381.448 arcsec.
##
##     Both the luni-solar and planetary nutations are included.  The
##     latter are due to direct planetary nutations and the
##     perturbations of the lunar and terrestrial orbits.
##
##  3) The function computes the MHB2000 nutation series with the
##     associated corrections for planetary nutations.  It is an
##     implementation of the nutation part of the IAU 2000A precession-
##     nutation model, formally adopted by the IAU General Assembly in
##     2000, namely MHB2000 (Mathews et al. 2002), but with the free
##     core nutation (FCN - see Note 4) omitted.
##
##  4) The full MHB2000 model also contains contributions to the
##     nutations in longitude and obliquity due to the free-excitation
##     of the free-core-nutation during the period 1979-2000.  These FCN
##     terms, which are time-dependent and unpredictable, are NOT
##     included in the present function and, if required, must be
##     independently computed.  With the FCN corrections included, the
##     present function delivers a pole which is at current epochs
##     accurate to a few hundred microarcseconds.  The omission of FCN
##     introduces further errors of about that size.
##
##  5) The present function provides classical nutation.  The MHB2000
##     algorithm, from which it is adapted, deals also with (i) the
##     offsets between the GCRS and mean poles and (ii) the adjustments
##     in longitude and obliquity due to the changed precession rates.
##     These additional functions, namely frame bias and precession
##     adjustments, are supported by the SOFA functions iauBi00  and
##     iauPr00.
##
##  6) The MHB2000 algorithm also provides "total" nutations, comprising
##     the arithmetic sum of the frame bias, precession adjustments,
##     luni-solar nutation and planetary nutation.  These total
##     nutations can be used in combination with an existing IAU 1976
##     precession implementation, such as iauPmat76,  to deliver GCRS-
##     to-true predictions of sub-mas accuracy at current dates.
##     However, there are three shortcomings in the MHB2000 model that
##     must be taken into account if more accurate or definitive results
##     are required (see Wallace 2002):
##
##       (i) The MHB2000 total nutations are simply arithmetic sums,
##           yet in reality the various components are successive Euler
##           rotations.  This slight lack of rigor leads to cross terms
##           that exceed 1 mas after a century.  The rigorous procedure
##           is to form the GCRS-to-true rotation matrix by applying the
##           bias, precession and nutation in that order.
##
##      (ii) Although the precession adjustments are stated to be with
##           respect to Lieske et al. (1977), the MHB2000 model does
##           not specify which set of Euler angles are to be used and
##           how the adjustments are to be applied.  The most literal
##           and straightforward procedure is to adopt the 4-rotation
##           epsilon_0, psi_A, omega_A, xi_A option, and to add DPSIPR
##           to psi_A and DEPSPR to both omega_A and eps_A.
##
##     (iii) The MHB2000 model predates the determination by Chapront
##           et al. (2002) of a 14.6 mas displacement between the
##           J2000.0 mean equinox and the origin of the ICRS frame.  It
##           should, however, be noted that neglecting this displacement
##           when calculating star coordinates does not lead to a
##           14.6 mas change in right ascension, only a small second-
##           order distortion in the pattern of the precession-nutation
##           effect.
##
##     For these reasons, the SOFA functions do not generate the "total
##     nutations" directly, though they can of course easily be
##     generated by calling iauBi00, iauPr00 and the present function
##     and adding the results.
##
##  7) The MHB2000 model contains 41 instances where the same frequency
##     appears multiple times, of which 38 are duplicates and three are
##     triplicates.  To keep the present code close to the original MHB
##     algorithm, this small inefficiency has not been corrected.
##
##  Called:
##     iauFal03     mean anomaly of the Moon
##     iauFaf03     mean argument of the latitude of the Moon
##     iauFaom03    mean longitude of the Moon's ascending node
##     iauFame03    mean longitude of Mercury
##     iauFave03    mean longitude of Venus
##     iauFae03     mean longitude of Earth
##     iauFama03    mean longitude of Mars
##     iauFaju03    mean longitude of Jupiter
##     iauFasa03    mean longitude of Saturn
##     iauFaur03    mean longitude of Uranus
##     iauFapa03    general accumulated precession in longitude
##
##  References:
##
##     Chapront, J., Chapront-Touze, M. & Francou, G. 2002,
##     Astron.Astrophys. 387, 700
##
##     Lieske, J.H., Lederle, T., Fricke, W. & Morando, B. 1977,
##     Astron.Astrophys. 58, 1-16
##
##     Mathews, P.M., Herring, T.A., Buffet, B.A. 2002, J.Geophys.Res.
##     107, B4.  The MHB_2000 code itself was obtained on 9th September
##     2002 from ftp//maia.usno.navy.mil/conv2000/chapter5/IAU2000A.
##
##     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
##     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
##
##     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
##     Astron.Astrophys.Supp.Ser. 135, 111
##
##     Wallace, P.T., "Software for Implementing the IAU 2000
##     Resolutions", in IERS Workshop 5.1 (2002)
##
##  This revision:  2019 June 23
##
##  SOFA release 2019-07-22
##
##  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
####################################

    date1 <- tt[,1]
    date2 <- tt[,2]
    U2R <- DAS2R/1e7
#Part I: Luni-Solar nutation model
###The units for the sine and cosine coefficients are
###0.1 microarcsecond and the same per Julian century
### coefficients of l,l',F,D,Om
    nl <- xls[,'nl']
    nlp <- xls[,'nlp']
    nf <- xls[,'nf']
    nd <- xls[,'nd']
    nom <- xls[,'nom']
                                        #longitude sin, t*sin, cos coefficients
    sp <- xls[,'sp']
    spt <- xls[,'spt']
    cp <- xls[,'cp']
### obliquity cos, t*cos, sin coefficients
    ce <- xls[,'ce']
    cet <- xls[,'cet']
    se <- xls[,'se']
                                        #Number of terms in the luni-solar nutation model
    NLS <- nrow(xls)

                                        # Interval between fundamental date J2000.0 and given date (JC). */
    ##Days per Julian century */
    DJC <- 36525.0
    DJ00 <- 2451545.0
    t <- ((date1 - DJ00) + date2) / DJC

    ##LUNI-SOLAR NUTATION */
    ##Fundamental (Delaunay) arguments */
    ##Mean anomaly of the Moon (IERS 2003). */
                                        #Arcseconds in a full circle
                                        #   el <- iauFal03(t);
    TURNAS <- 1296000.0
    el <- 485868.249036+t * ( 1717915923.2178 +t * (31.8792 +t * (0.051635 +t * (-0.00024470 ) ) ) )
    el <- (el%%TURNAS)*DAS2R

                                        #Mean anomaly of the Sun (MHB2000).
    elp <- 1287104.79305 +t * (129596581.0481  +t * (-0.5532  +t * (0.000136  +t * (-0.00001149))))
    elp <- (elp%%TURNAS)*DAS2R

    ##Mean longitude of the Moon minus that of the ascending node(IERS 2003)
    ##f = iauFaf03(t);
    f <- 335779.526232 +t * ( 1739527262.8478 +t * (-12.7512 +t * (-0.001037 +t * (0.00000417 ) ) ) )
    f <- (f%%TURNAS)*DAS2R

    ##Mean elongation of the Moon from the Sun (MHB2000).
    d <- 1072260.70369  +t * (1602961601.2090  +t * (-6.3706  +t * (0.006593  +t * (-0.00003169))))
    d <- (d%%TURNAS)*DAS2R

    ##Mean longitude of the ascending node of the Moon (IERS 2003).
                                        #om = iauFaom03(t);
    om <- 450160.398036 +t*(-6962890.5431+t*(7.4722 +t*(0.007702 +t*(-0.00005939 ) ) ) )
    om <- (om%%TURNAS)*DAS2R

    ##Initialize the nutation values.
    dp <- 0.0
    de <- 0.0

    ##Summation of luni-solar nutation series (in reverse order).
    for (i in NLS:1){
        ##Argument and functions.
        arg <- xls[i,'nl']  * el +nlp[i] * elp +nf[i]  * f +nd[i]  * d +nom[i] * om
        arg <- arg%%D2PI
        sarg <- sin(arg)
        carg <- cos(arg)

###Term
        dp <- dp+(sp[i] + spt[i] * t) * sarg + cp[i] * carg
        de <- de+(ce[i] + cet[i] * t) * carg + se[i] * sarg
    }

    ##Convert from 0.1 microarcsec units to radians.
    dpsils <- dp * U2R
    depsls <- de * U2R

###
###PLANETARY NUTATION
    ## n.b.  The MHB2000 code computes the luni-solar and planetary nutation */
    ## in different functions, using slightly different Delaunay */
    ## arguments in the two cases.  This behaviour is faithfully */
    ## reproduced here.  Use of the IERS 2003 expressions for both */
    ## cases leads to negligible changes, well below */
    ## 0.1 microarcsecond. */

    ##Mean anomaly of the Moon (MHB2000). */
    al <- 2.35555598 + 8328.6914269554 * t
    al <- al%%D2PI

    ## Mean longitude of the Moon minus that of the ascending node (MHB2000).
    af <- 1.627905234 + 8433.466158131 * t
    af <- af%%D2PI

    ## Mean elongation of the Moon from the Sun (MHB2000).
    ad <- 5.198466741 + 7771.3771468121 * t
    ad <- ad%%D2PI

    ## Mean longitude of the ascending node of the Moon (MHB2000). */
    aom <- (2.18243920 - 33.757045 * t)%%D2PI

    ## General accumulated precession in longitude (IERS 2003). */
    apa <- (0.024381750 + 0.00000538691 * t) * t

    ## Planetary longitudes, Mercury through Uranus (IERS 2003). */
    alme <- (4.402608842 + 2608.7903141574 * t)%%D2PI
    alve <- (3.176146697 + 1021.3285546211 * t)%%D2PI
    alea <- (1.753470314 + 628.3075849991 * t)%%D2PI
    alma <- (6.203480913 + 334.0612426700 * t)%%D2PI
    alju <- (0.599546497 + 52.9690962641 * t)%%D2PI
    alsa <- (0.874016757 + 21.3299104960 * t)%%D2PI
    alur <- (5.481293872 + 7.4781598567 * t)%%D2PI

    ## Neptune longitude (MHB2000).
    alne <- (5.321159000 + 3.8127774000 * t)%%D2PI


######Part II: Planetary nutation model
######The units for the sine and cosine coefficients are
######0.1 microarcsecond
###xpl: coefficients of l,l',F,D,Om
### coefficients of longitude sin, t*sin, cos coefficients
### obliquity cos, t*cos, sin coefficients
    nl <- xpl[,'nl']
    nf <- xpl[,'nf']
    nd <- xpl[,'nd']
    nom <- xpl[,'nom']
    ## /* coefficients of planetary longitudes */
    nme <- xpl[,'nme']
    nve <- xpl[,'nve']
    nea <- xpl[,'nea']
    nma <- xpl[,'nma']
    nju <- xpl[,'nju']
    nsa <- xpl[,'nsa']
    nur <- xpl[,'nur']
    nne <- xpl[,'nne']
    ## /* coefficient of general precession */
    npa <- xpl[,'npa']
    ## /* longitude sin, cos coefficients */
    sp <- xpl[,'sp']
    cp <- xpl[,'cp']
    ##/* obliquity sin, cos coefficients */
    ce <- xpl[,'ce']
    se <- xpl[,'se']
    NPL <- nrow(xpl)

    ## Initialize the nutation values. */
    dp <- 0.0
    de <- 0.0
    ## Summation of planetary nutation series (in reverse order). */
    for (i in NPL:1) {
        ##Argument and functions. */
        arg <- (nl[i]  * al   +
                   nf[i]  * af   +
                       nd[i]  * ad   +
                           nom[i] * aom  +
                               nme[i] * alme +
                                   nve[i] * alve +
                                       nea[i] * alea +
                                           nma[i] * alma +
                                               nju[i] * alju +
                                                   nsa[i] * alsa +
                                                       nur[i] * alur +
                                                           nne[i] * alne +
                                                               npa[i] * apa)%%D2PI
        sarg <- sin(arg)
        carg <- cos(arg)

        ##term
        dp <- dp+sp[i] * sarg + cp[i] * carg
        de <- de+se[i] * sarg + ce[i] * carg
    }

    ## Convert from 0.1 microarcsec units to radians. */
    dpsipl <- dp * U2R
    depspl <- de * U2R

    ##results
    ##Add luni-solar and planetary components.
    dpsi <- dpsils + dpsipl
    deps <- depsls + depspl
    return(list(dpsi=dpsi,deps=deps))
}

sofa_Nut00b <- function(tt){
####################################
##  Nutation, IAU 2000B model.
##
##  This function is part of the International Astronomical Union's
##  SOFA (Standards Of Fundamental Astronomy) software collection.
##
##  Status:  canonical model.
##
##  Given:
##     date1,date2   double    TT as a 2-part Julian Date (Note 1)
##
##  Returned:
##     dpsi,deps     double    nutation, luni-solar + planetary (Note 2)
##
##  Notes:
##
##  1) The TT date date1+date2 is a Julian Date, apportioned in any
##     convenient way between the two arguments.  For example,
##     JD(TT)=2450123.7 could be expressed in any of these ways,
##     among others:
##
##            date1          date2
##
##         2450123.7           0.0       (JD method)
##         2451545.0       -1421.3       (J2000 method)
##         2400000.5       50123.2       (MJD method)
##         2450123.5           0.2       (date & time method)
##
##     The JD method is the most natural and convenient to use in
##     cases where the loss of several decimal digits of resolution
##     is acceptable.  The J2000 method is best matched to the way
##     the argument is handled internally and will deliver the
##     optimum resolution.  The MJD method and the date & time methods
##     are both good compromises between resolution and convenience.
##
##  2) The nutation components in longitude and obliquity are in radians
##     and with respect to the equinox and ecliptic of date.  The
##     obliquity at J2000.0 is assumed to be the Lieske et al. (1977)
##     value of 84381.448 arcsec.  (The errors that result from using
##     this function with the IAU 2006 value of 84381.406 arcsec can be
##     neglected.)
##
##     The nutation model consists only of luni-solar terms, but
##     includes also a fixed offset which compensates for certain long-
##     period planetary terms (Note 7).
##
##  3) This function is an implementation of the IAU 2000B abridged
##     nutation model formally adopted by the IAU General Assembly in
##     2000.  The function computes the MHB_2000_SHORT luni-solar
##     nutation series (Luzum 2001), but without the associated
##     corrections for the precession rate adjustments and the offset
##     between the GCRS and J2000.0 mean poles.
##
##  4) The full IAU 2000A (MHB2000) nutation model contains nearly 1400
##     terms.  The IAU 2000B model (McCarthy & Luzum 2003) contains only
##     77 terms, plus additional simplifications, yet still delivers
##     results of 1 mas accuracy at present epochs.  This combination of
##     accuracy and size makes the IAU 2000B abridged nutation model
##     suitable for most practical applications.
##
##     The function delivers a pole accurate to 1 mas from 1900 to 2100
##     (usually better than 1 mas, very occasionally just outside
##     1 mas).  The full IAU 2000A model, which is implemented in the
##     function iauNut00a (q.v.), delivers considerably greater accuracy
##     at current dates;  however, to realize this improved accuracy,
##     corrections for the essentially unpredictable free-core-nutation
##     (FCN) must also be included.
##
##  5) The present function provides classical nutation.  The
##     MHB_2000_SHORT algorithm, from which it is adapted, deals also
##     with (i) the offsets between the GCRS and mean poles and (ii) the
##     adjustments in longitude and obliquity due to the changed
##     precession rates.  These additional functions, namely frame bias
##     and precession adjustments, are supported by the SOFA functions
##     iauBi00  and iauPr00.
##
##  6) The MHB_2000_SHORT algorithm also provides "total" nutations,
##     comprising the arithmetic sum of the frame bias, precession
##     adjustments, and nutation (luni-solar + planetary).  These total
##     nutations can be used in combination with an existing IAU 1976
##     precession implementation, such as iauPmat76,  to deliver GCRS-
##     to-true predictions of mas accuracy at current epochs.  However,
##     for symmetry with the iauNut00a  function (q.v. for the reasons),
##     the SOFA functions do not generate the "total nutations"
##     directly.  Should they be required, they could of course easily
##     be generated by calling iauBi00, iauPr00 and the present function
##     and adding the results.
##
##  7) The IAU 2000B model includes "planetary bias" terms that are
##     fixed in size but compensate for long-period nutations.  The
##     amplitudes quoted in McCarthy & Luzum (2003), namely
##     Dpsi = -1.5835 mas and Depsilon = +1.6339 mas, are optimized for
##     the "total nutations" method described in Note 6.  The Luzum
##     (2001) values used in this SOFA implementation, namely -0.135 mas
##     and +0.388 mas, are optimized for the "rigorous" method, where
##     frame bias, precession and nutation are applied separately and in
##     that order.  During the interval 1995-2050, the SOFA
##     implementation delivers a maximum error of 1.001 mas (not
##     including FCN).
##
##  References:
##
##     Lieske, J.H., Lederle, T., Fricke, W., Morando, B., "Expressions
##     for the precession quantities based upon the IAU /1976/ system of
##     astronomical constants", Astron.Astrophys. 58, 1-2, 1-16. (1977)
##
##     Luzum, B., private communication, 2001 (Fortran code
##     MHB_2000_SHORT)
##
##     McCarthy, D.D. & Luzum, B.J., "An abridged model of the
##     precession-nutation of the celestial pole", Cel.Mech.Dyn.Astron.
##     85, 37-49 (2003)
##
##     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
##     Francou, G., Laskar, J., Astron.Astrophys. 282, 663-683 (1994)
##
##  This revision:  2019 June 23
##
##  SOFA release 2019-07-22
##
##  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
####################################
    date1 <- tt[,1]
    date2 <- tt[,2]
    U2R <- DAS2R/1e7
    DPPLAN = -0.135 * DMAS2R
    DEPLAN =  0.388 * DMAS2R
       #nl,nlp,nf,nd,nom; /* coefficients of l,l',F,D,Om */
      #ps,pst,pc;     /* longitude sin, t*sin, cos coefficients */
      # ec,ect,es;     /* obliquity cos, t*cos, sin coefficients */
x <- rbind(
      c( 0, 0, 0, 0,1,-172064161.0, -174666.0, 33386.0, 92052331.0, 9086.0, 15377.0),
      c( 0, 0, 2,-2,2, -13170906.0, -1675.0, -13696.0, 5730336.0, -3015.0, -4587.0),
      c( 0, 0, 2, 0,2,-2276413.0,-234.0, 2796.0, 978459.0,-485.0,1374.0),
      c( 0, 0, 0, 0,2,2074554.0,  207.0, -698.0,-897492.0, 470.0,-291.0),
      c( 0, 1, 0, 0,0,1475877.0,-3633.0,11817.0, 73871.0,-184.0,-1924.0),
      c( 0, 1, 2,-2,2,-516821.0, 1226.0, -524.0, 224386.0,-677.0,-174.0),
      c( 1, 0, 0, 0,0, 711159.0,   73.0, -872.0,  -6750.0,   0.0, 358.0),
      c( 0, 0, 2, 0,1,-387298.0, -367.0,  380.0, 200728.0,  18.0, 318.0),
      c( 1, 0, 2, 0,2,-301461.0,  -36.0,  816.0, 129025.0, -63.0, 367.0),
      c( 0,-1, 2,-2,2, 215829.0, -494.0,  111.0, -95929.0, 299.0, 132.0),
      c( 0, 0, 2,-2,1, 128227.0,  137.0,  181.0, -68982.0,  -9.0,  39.0),
      c(-1, 0, 2, 0,2, 123457.0,   11.0,   19.0, -53311.0,  32.0,  -4.0),
      c(-1, 0, 0, 2,0, 156994.0,   10.0, -168.0,  -1235.0,   0.0,  82.0),
      c( 1, 0, 0, 0,1,  63110.0,   63.0,   27.0, -33228.0,   0.0,  -9.0),
      c(-1, 0, 0, 0,1, -57976.0,  -63.0, -189.0,  31429.0,   0.0, -75.0),
      c(-1, 0, 2, 2,2, -59641.0,  -11.0,  149.0,  25543.0, -11.0,  66.0),
      c( 1, 0, 2, 0,1, -51613.0,  -42.0,  129.0,  26366.0,   0.0,  78.0),
      c(-2, 0, 2, 0,1,  45893.0,   50.0,   31.0, -24236.0, -10.0,  20.0),
      c( 0, 0, 0, 2,0,  63384.0,   11.0, -150.0,  -1220.0,   0.0,  29.0),
      c( 0, 0, 2, 2,2, -38571.0,   -1.0,  158.0,  16452.0, -11.0,  68.0),
      c( 0,-2, 2,-2,2,  32481.0,    0.0,    0.0, -13870.0,   0.0,   0.0),
      c(-2, 0, 0, 2,0, -47722.0,    0.0,  -18.0,    477.0,   0.0, -25.0),
      c( 2, 0, 2, 0,2, -31046.0,   -1.0,  131.0,  13238.0, -11.0,  59.0),
      c( 1, 0, 2,-2,2,  28593.0,    0.0,   -1.0, -12338.0,  10.0,  -3.0),
      c(-1, 0, 2, 0,1,  20441.0,   21.0,   10.0, -10758.0,   0.0,  -3.0),
      c( 2, 0, 0, 0,0,  29243.0,    0.0,  -74.0,   -609.0,   0.0,  13.0),
      c( 0, 0, 2, 0,0,  25887.0,    0.0,  -66.0,   -550.0,   0.0,  11.0),
      c( 0, 1, 0, 0,1, -14053.0,  -25.0,   79.0,   8551.0,  -2.0, -45.0),
      c(-1, 0, 0, 2,1,  15164.0,   10.0,   11.0,  -8001.0,   0.0,  -1.0),
      c( 0, 2, 2,-2,2, -15794.0,   72.0,  -16.0,   6850.0, -42.0,  -5.0),
      c( 0, 0,-2, 2,0,  21783.0,    0.0,   13.0,   -167.0,   0.0,  13.0),
      c( 1, 0, 0,-2,1, -12873.0,  -10.0,  -37.0,   6953.0,   0.0, -14.0),
      c( 0,-1, 0, 0,1, -12654.0,   11.0,   63.0,   6415.0,   0.0,  26.0),
      c(-1, 0, 2, 2,1, -10204.0,    0.0,   25.0,   5222.0,   0.0,  15.0),
      c( 0, 2, 0, 0,0,  16707.0,  -85.0,  -10.0,    168.0,  -1.0,  10.0),
      c( 1, 0, 2, 2,2,  -7691.0,    0.0,   44.0,   3268.0,   0.0,  19.0),
      c(-2, 0, 2, 0,0, -11024.0,    0.0,  -14.0,    104.0,   0.0,   2.0),
      c( 0, 1, 2, 0,2,   7566.0,  -21.0,  -11.0,  -3250.0,   0.0,  -5.0),
      c( 0, 0, 2, 2,1,  -6637.0,  -11.0,   25.0,   3353.0,   0.0,  14.0),
      c( 0,-1, 2, 0,2,  -7141.0,   21.0,    8.0,   3070.0,   0.0,   4.0),
      c( 0, 0, 0, 2,1,  -6302.0,  -11.0,    2.0,   3272.0,   0.0,   4.0),
      c( 1, 0, 2,-2,1,   5800.0,   10.0,    2.0,  -3045.0,   0.0,  -1.0),
      c( 2, 0, 2,-2,2,   6443.0,    0.0,   -7.0,  -2768.0,   0.0,  -4.0),
      c(-2, 0, 0, 2,1,  -5774.0,  -11.0,  -15.0,   3041.0,   0.0,  -5.0),
      c( 2, 0, 2, 0,1,  -5350.0,    0.0,   21.0,   2695.0,   0.0,  12.0),
      c( 0,-1, 2,-2,1,  -4752.0,  -11.0,   -3.0,   2719.0,   0.0,  -3.0),
      c( 0, 0, 0,-2,1,  -4940.0,  -11.0,  -21.0,   2720.0,   0.0,  -9.0),
      c(-1,-1, 0, 2,0,   7350.0,    0.0,   -8.0,    -51.0,   0.0,   4.0),
      c( 2, 0, 0,-2,1,   4065.0,    0.0,    6.0,  -2206.0,   0.0,   1.0),
      c( 1, 0, 0, 2,0,   6579.0,    0.0,  -24.0,   -199.0,   0.0,   2.0),
      c( 0, 1, 2,-2,1,   3579.0,    0.0,    5.0,  -1900.0,   0.0,   1.0),
      c( 1,-1, 0, 0,0,   4725.0,    0.0,   -6.0,    -41.0,   0.0,   3.0),
      c(-2, 0, 2, 0,2,  -3075.0,    0.0,   -2.0,   1313.0,   0.0,  -1.0),
      c( 3, 0, 2, 0,2,  -2904.0,    0.0,   15.0,   1233.0,   0.0,   7.0),
      c( 0,-1, 0, 2,0,   4348.0,    0.0,  -10.0,    -81.0,   0.0,   2.0),
      c( 1,-1, 2, 0,2,  -2878.0,    0.0,    8.0,   1232.0,   0.0,   4.0),
      c( 0, 0, 0, 1,0,  -4230.0,    0.0,    5.0,    -20.0,   0.0,  -2.0),
      c(-1,-1, 2, 2,2,  -2819.0,    0.0,    7.0,   1207.0,   0.0,   3.0),
      c(-1, 0, 2, 0,0,  -4056.0,    0.0,    5.0,     40.0,   0.0,  -2.0),
      c( 0,-1, 2, 2,2,  -2647.0,    0.0,   11.0,   1129.0,   0.0,   5.0),
      c(-2, 0, 0, 0,1,  -2294.0,    0.0,  -10.0,   1266.0,   0.0,  -4.0),
      c( 1, 1, 2, 0,2,   2481.0,    0.0,   -7.0,  -1062.0,   0.0,  -3.0),
      c( 2, 0, 0, 0,1,   2179.0,    0.0,   -2.0,  -1129.0,   0.0,  -2.0),
      c(-1, 1, 0, 1,0,   3276.0,    0.0,    1.0,     -9.0,   0.0,   0.0),
      c( 1, 1, 0, 0,0,  -3389.0,    0.0,    5.0,     35.0,   0.0,  -2.0),
      c( 1, 0, 2, 0,0,   3339.0,    0.0,  -13.0,   -107.0,   0.0,   1.0),
      c(-1, 0, 2,-2,1,  -1987.0,    0.0,   -6.0,   1073.0,   0.0,  -2.0),
      c( 1, 0, 0, 0,2,  -1981.0,    0.0,    0.0,    854.0,   0.0,   0.0),
      c(-1, 0, 0, 1,0,   4026.0,    0.0, -353.0,   -553.0,   0.0,-139.0),
      c( 0, 0, 2, 1,2,   1660.0,    0.0,   -5.0,   -710.0,   0.0,  -2.0),
      c(-1, 0, 2, 4,2,  -1521.0,    0.0,    9.0,    647.0,   0.0,   4.0),
      c(-1, 1, 0, 1,1,   1314.0,    0.0,    0.0,   -700.0,   0.0,   0.0),
      c( 0,-2, 2,-2,1,  -1283.0,    0.0,    0.0,    672.0,   0.0,   0.0),
      c( 1, 0, 2, 2,1,  -1331.0,    0.0,    8.0,    663.0,   0.0,   4.0),
      c(-2, 0, 2, 2,2,   1383.0,    0.0,   -2.0,   -594.0,   0.0,  -2.0),
      c(-1, 0, 0, 0,2,   1405.0,    0.0,    4.0,   -610.0,   0.0,   2.0),
      c( 1, 1, 2,-2,2,   1290.0,    0.0,    0.0,   -556.0,   0.0,   0.0)
   )
colnames(x) <- c('nl','nlp','nf','nd','nom','ps','pst','pc','ec','ect','es')

##Number of terms in the series */
   NLS = nrow(x)

## Interval between fundamental epoch J2000.0 and given date (JC). */
   t <- ((date1 - DJ00) + date2) / DJC

########################
##LUNI-SOLAR NUTATION */
########################
##Fundamental (Delaunay) arguments from Simon et al. (1994)
##Mean anomaly of the Moon.
  el <- ((485868.249036 + (1717915923.2178) * t)%%TURNAS)*DAS2R

##Mean anomaly of the Sun. */
   elp <- ((1287104.79305 + (129596581.0481) * t)%%TURNAS)*DAS2R

##Mean argument of the latitude of the Moon. */
   f <- ((335779.526232 + (1739527262.8478) * t)%%TURNAS)*DAS2R

##Mean elongation of the Moon from the Sun. */
   d <- ((1072260.70369 + (1602961601.2090) * t)%%TURNAS)*DAS2R

##Mean longitude of the ascending node of the Moon. */
   om <- ((450160.398036 + (-6962890.5431) * t)%%TURNAS)*DAS2R

##Initialize the nutation values.
   dp = 0.0
   de = 0.0

##Summation of luni-solar nutation series (smallest terms first).
   for (i in NLS:1) {
##Argument and functions. */
      arg = ( x[i,'nl']  * el  +
                  x[i,'nlp']*elp +
                  x[i,'nf']*f   +
                  x[i,'nd']*d   +
                  x[i,'nom']*om)%%D2PI
      sarg = sin(arg)
      carg = cos(arg)

##Term.
      dp <- dp+(x[i,'ps'] + x[i,'pst'] * t) * sarg + x[i,'pc'] * carg
      de <- de+(x[i,'ec'] + x[i,'ect'] * t) * carg + x[i,'es'] * sarg
   }

##Convert from 0.1 microarcsec units to radians.
   dpsils = dp * U2R
   depsls = de * U2R

###IN LIEU OF PLANETARY NUTATION
##Fixed offset to correct for missing terms in truncated series. */
   dpsipl = DPPLAN
   depspl = DEPLAN

##Add luni-solar and planetary components.
   dpsi <- as.numeric(dpsils + dpsipl)
   deps <- as.numeric(depsls + depspl)
    return(list(dpsi=dpsi,deps=deps))
}
sofa_Dtdb <- function(date1, date2, ut, elong, u, v){
####################################
##  An approximation to TDB-TT, the difference between barycentric
##  dynamical time and terrestrial time, for an observer on the Earth.
##
##  The different time scales - proper, coordinate and realized - are
##  related to each other:
##
##            TAI             <-  physically realized
##             :
##          offset            <-  observed (nominally +32.184s)
##             :
##            TT              <-  terrestrial time
##             :
##    rate adjustment (L_G)   <-  definition of TT
##             :
##            TCG             <-  time scale for GCRS
##             :
##      "periodic" terms      <-  iauDtdb  is an implementation
##             :
##    rate adjustment (L_C)   <-  function of solar-system ephemeris
##             :
##            TCB             <-  time scale for BCRS
##             :
##    rate adjustment (-L_B)  <-  definition of TDB
##             :
##            TDB             <-  TCB scaled to track TT
##             :
##      "periodic" terms      <-  -iauDtdb is an approximation
##             :
##            TT              <-  terrestrial time
##
##  Adopted values for the various constants can be found in the IERS
##  Conventions (McCarthy & Petit 2003).
##
##  This function is part of the International Astronomical Union's
##  SOFA (Standards Of Fundamental Astronomy) software collection.
##
##  Status:  support function.
##
##  Given:
##     date1,date2   double  date, TDB (Notes 1-3)
##     ut            double  universal time (UT1, fraction of one day)
##     elong         double  longitude (east positive, radians)
##     u             double  distance from Earth spin axis (km)
##     v             double  distance north of equatorial plane (km)
##
##  Returned (function value):
##                   double  TDB-TT (seconds)
##
##  Notes:
##
##  1) The date date1+date2 is a Julian Date, apportioned in any
##     convenient way between the two arguments.  For example,
##     JD(TT)=2450123.7 could be expressed in any of these ways,
##     among others:
##
##            date1          date2
##
##         2450123.7           0.0       (JD method)
##         2451545.0       -1421.3       (J2000 method)
##         2400000.5       50123.2       (MJD method)
##         2450123.5           0.2       (date & time method)
##
##     The JD method is the most natural and convenient to use in
##     cases where the loss of several decimal digits of resolution
##     is acceptable.  The J2000 method is best matched to the way
##     the argument is handled internally and will deliver the
##     optimum resolution.  The MJD method and the date & time methods
##     are both good compromises between resolution and convenience.
##
##     Although the date is, formally, barycentric dynamical time (TDB),
##     the terrestrial dynamical time (TT) can be used with no practical
##     effect on the accuracy of the prediction.
##
##  2) TT can be regarded as a coordinate time that is realized as an
##     offset of 32.184s from International Atomic Time, TAI.  TT is a
##     specific linear transformation of geocentric coordinate time TCG,
##     which is the time scale for the Geocentric Celestial Reference
##     System, GCRS.
##
##  3) TDB is a coordinate time, and is a specific linear transformation
##     of barycentric coordinate time TCB, which is the time scale for
##     the Barycentric Celestial Reference System, BCRS.
##
##  4) The difference TCG-TCB depends on the masses and positions of the
##     bodies of the solar system and the velocity of the Earth.  It is
##     dominated by a rate difference, the residual being of a periodic
##     character.  The latter, which is modeled by the present function,
##     comprises a main (annual) sinusoidal term of amplitude
##     approximately 0.00166 seconds, plus planetary terms up to about
##     20 microseconds, and lunar and diurnal terms up to 2 microseconds.
##     These effects come from the changing transverse Doppler effect
##     and gravitational red-shift as the observer (on the Earth's
##     surface) experiences variations in speed (with respect to the
##     BCRS) and gravitational potential.
##
##  5) TDB can be regarded as the same as TCB but with a rate adjustment
##     to keep it close to TT, which is convenient for many applications.
##     The history of successive attempts to define TDB is set out in
##     Resolution 3 adopted by the IAU General Assembly in 2006, which
##     defines a fixed TDB(TCB) transformation that is consistent with
##     contemporary solar-system ephemerides.  Future ephemerides will
##     imply slightly changed transformations between TCG and TCB, which
##     could introduce a linear drift between TDB and TT;  however, any
##     such drift is unlikely to exceed 1 nanosecond per century.
##
##  6) The geocentric TDB-TT model used in the present function is that of
##     Fairhead & Bretagnon (1990), in its full form.  It was originally
##     supplied by Fairhead (private communications with P.T.Wallace,
##     1990) as a Fortran subroutine.  The present C function contains an
##     adaptation of the Fairhead code.  The numerical results are
##     essentially unaffected by the changes, the differences with
##     respect to the Fairhead & Bretagnon original being at the 1e-20 s
##     level.
##
##     The topocentric part of the model is from Moyer (1981) and
##     Murray (1983), with fundamental arguments adapted from
##     Simon et al. 1994.  It is an approximation to the expression
##     ( v / c ) . ( r / c ), where v is the barycentric velocity of
##     the Earth, r is the geocentric position of the observer and
##     c is the speed of light.
##
##     By supplying zeroes for u and v, the topocentric part of the
##     model can be nullified, and the function will return the Fairhead
##     & Bretagnon result alone.
##
##  7) During the interval 1950-2050, the absolute accuracy is better
##     than +/- 3 nanoseconds relative to time ephemerides obtained by
##     direct numerical integrations based on the JPL DE405 solar system
##     ephemeris.
##
##  8) It must be stressed that the present function is merely a model,
##     and that numerical integration of solar-system ephemerides is the
##     definitive method for predicting the relationship between TCG and
##     TCB and hence between TT and TDB.
##
##  References:
##
##     Fairhead, L., & Bretagnon, P., Astron.Astrophys., 229, 240-247
##     (1990).
##
##     IAU 2006 Resolution 3.
##
##     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
##     IERS Technical Note No. 32, BKG (2004)
##
##     Moyer, T.D., Cel.Mech., 23, 33 (1981).
##
##     Murray, C.A., Vectorial Astrometry, Adam Hilger (1983).
##
##     Seidelmann, P.K. et al., Explanatory Supplement to the
##     Astronomical Almanac, Chapter 2, University Science Books (1992).
##
##     Simon, J.L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
##     Francou, G. & Laskar, J., Astron.Astrophys., 282, 663-683 (1994).
##
##  This revision:  2018 January 2
##
##  SOFA release 2019-07-22
##
##  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
####################################
    fairhd <- rbind(
        ## 1, 10 ##
        c( 1656.674564e-6,     6283.075849991,  6.240054195 ),
        c(   22.417471e-6,     5753.384884897,  4.296977442 ),
        c(   13.839792e-6,    12566.151699983,  6.196904410 ),
        c(    4.770086e-6,      529.690965095,  0.444401603 ),
        c(    4.676740e-6,     6069.776754553,  4.021195093 ),
        c(    2.256707e-6,      213.299095438,  5.543113262 ),
        c(    1.694205e-6,      -3.523118349,   5.025132748 ),
        c(    1.554905e-6,    77713.771467920,  5.198467090 ),
        c(    1.276839e-6,     7860.419392439,  5.988822341 ),
        c(    1.193379e-6,     5223.693919802,  3.649823730 ),
        ## 11, 20 ##
        c(    1.115322e-6,     3930.209696220,  1.422745069 ),
        c(    0.794185e-6,    11506.769769794,  2.322313077 ),
        c(    0.447061e-6,       26.298319800,  3.615796498 ),
        c(    0.435206e-6,     -398.149003408,  4.349338347 ),
        c(    0.600309e-6,     1577.343542448,  2.678271909 ),
        c(    0.496817e-6,     6208.294251424,  5.696701824 ),
        c(    0.486306e-6,     5884.926846583,  0.520007179 ),
        c(    0.432392e-6,       74.781598567,  2.435898309 ),
        c(    0.468597e-6,     6244.942814354,  5.866398759 ),
        c(    0.375510e-6,     5507.553238667,  4.103476804 ),
        ## 21, 30 ##
        c(    0.243085e-6,     -775.522611324,  3.651837925 ),
        c(    0.173435e-6,    18849.227549974,  6.153743485 ),
        c(    0.230685e-6,     5856.477659115,  4.773852582 ),
        c(    0.203747e-6,    12036.460734888,  4.333987818 ),
        c(    0.143935e-6,     -796.298006816,  5.957517795 ),
        c(    0.159080e-6,    10977.078804699,  1.890075226 ),
        c(    0.119979e-6,       38.133035638,  4.551585768 ),
        c(    0.118971e-6,     5486.777843175,  1.914547226 ),
        c(    0.116120e-6,     1059.381930189,  0.873504123 ),
        c(    0.137927e-6,    11790.629088659,  1.135934669 ),
        ## 31, 40 ##
        c(    0.098358e-6,     2544.314419883,  0.092793886 ),
        c(    0.101868e-6,    -5573.142801634,  5.984503847 ),
        c(    0.080164e-6,      206.185548437,  2.095377709 ),
        c(    0.079645e-6,     4694.002954708,  2.949233637 ),
        c(    0.062617e-6,       20.775395492,  2.654394814 ),
        c(    0.075019e-6,     2942.463423292,  4.980931759 ),
        c(    0.064397e-6,     5746.271337896,  1.280308748 ),
        c(    0.063814e-6,     5760.498431898,  4.167901731 ),
        c(    0.048042e-6,     2146.165416475,  1.495846011 ),
        c(    0.048373e-6,      155.420399434,  2.251573730 ),
        ## 41, 50 ##
        c(    0.058844e-6,      426.598190876,  4.839650148 ),
        c(    0.046551e-6,       -0.980321068,  0.921573539 ),
        c(    0.054139e-6,    17260.154654690,  3.411091093 ),
        c(    0.042411e-6,     6275.962302991,  2.869567043 ),
        c(    0.040184e-6,       -7.113547001,  3.565975565 ),
        c(    0.036564e-6,     5088.628839767,  3.324679049 ),
        c(    0.040759e-6,    12352.852604545,  3.981496998 ),
        c(    0.036507e-6,      801.820931124,  6.248866009 ),
        c(    0.036955e-6,     3154.687084896,  5.071801441 ),
        c(    0.042732e-6,      632.783739313,  5.720622217 ),
        ## 51, 60 ##
        c(    0.042560e-6,   161000.685737473,  1.270837679 ),
        c(    0.040480e-6,    15720.838784878,  2.546610123 ),
        c(    0.028244e-6,    -6286.598968340,  5.069663519 ),
        c(    0.033477e-6,     6062.663207553,  4.144987272 ),
        c(    0.034867e-6,      522.577418094,  5.210064075 ),
        c(    0.032438e-6,     6076.890301554,  0.749317412 ),
        c(    0.030215e-6,     7084.896781115,  3.389610345 ),
        c(    0.029247e-6,   -71430.695617928,  4.183178762 ),
        c(    0.033529e-6,     9437.762934887,  2.404714239 ),
        c(    0.032423e-6,     8827.390269875,  5.541473556 ),
        ## 61, 70 ##
        c(    0.027567e-6,     6279.552731642,  5.040846034 ),
        c(    0.029862e-6,    12139.553509107,  1.770181024 ),
        c(    0.022509e-6,    10447.387839604,  1.460726241 ),
        c(    0.020937e-6,     8429.241266467,  0.652303414 ),
        c(    0.020322e-6,      419.484643875,  3.735430632 ),
        c(    0.024816e-6,    -1194.447010225,  1.087136918 ),
        c(    0.025196e-6,     1748.016413067,  2.901883301 ),
        c(    0.021691e-6,    14143.495242431,  5.952658009 ),
        c(    0.017673e-6,     6812.766815086,  3.186129845 ),
        c(    0.022567e-6,     6133.512652857,  3.307984806 ),
        ## 71, 80 ##
        c(    0.016155e-6,    10213.285546211,  1.331103168 ),
        c(    0.014751e-6,     1349.867409659,  4.308933301 ),
        c(    0.015949e-6,     -220.412642439,  4.005298270 ),
        c(    0.015974e-6,    -2352.866153772,  6.145309371 ),
        c(    0.014223e-6,    17789.845619785,  2.104551349 ),
        c(    0.017806e-6,       73.297125859,  3.475975097 ),
        c(    0.013671e-6,     -536.804512095,  5.971672571 ),
        c(    0.011942e-6,     8031.092263058,  2.053414715 ),
        c(    0.014318e-6,    16730.463689596,  3.016058075 ),
        c(    0.012462e-6,      103.092774219,  1.737438797 ),
        ## 81, 90 ##
        c(    0.010962e-6,        3.590428652,  2.196567739 ),
        c(    0.015078e-6,    19651.048481098,  3.969480770 ),
        c(    0.010396e-6,      951.718406251,  5.717799605 ),
        c(    0.011707e-6,    -4705.732307544,  2.654125618 ),
        c(    0.010453e-6,     5863.591206116,  1.913704550 ),
        c(    0.012420e-6,     4690.479836359,  4.734090399 ),
        c(    0.011847e-6,     5643.178563677,  5.489005403 ),
        c(    0.008610e-6,     3340.612426700,  3.661698944 ),
        c(    0.011622e-6,     5120.601145584,  4.863931876 ),
        c(    0.010825e-6,      553.569402842,  0.842715011 ),
        ## 91, 100 ##
        c(    0.008666e-6,     -135.065080035,  3.293406547 ),
        c(    0.009963e-6,      149.563197135,  4.870690598 ),
        c(    0.009858e-6,     6309.374169791,  1.061816410 ),
        c(    0.007959e-6,      316.391869657,  2.465042647 ),
        c(    0.010099e-6,      283.859318865,  1.942176992 ),
        c(    0.007147e-6,     -242.728603974,  3.661486981 ),
        c(    0.007505e-6,     5230.807466803,  4.920937029 ),
        c(    0.008323e-6,    11769.853693166,  1.229392026 ),
        c(    0.007490e-6,    -6256.777530192,  3.658444681 ),
        c(    0.009370e-6,   149854.400134205,  0.673880395 ),
        ## 101, 110 ##
        c(    0.007117e-6,       38.027672636,  5.294249518 ),
        c(    0.007857e-6,    12168.002696575,  0.525733528 ),
        c(    0.007019e-6,     6206.809778716,  0.837688810 ),
        c(    0.006056e-6,      955.599741609,  4.194535082 ),
        c(    0.008107e-6,    13367.972631107,  3.793235253 ),
        c(    0.006731e-6,     5650.292110678,  5.639906583 ),
        c(    0.007332e-6,       36.648562930,  0.114858677 ),
        c(    0.006366e-6,     4164.311989613,  2.262081818 ),
        c(    0.006858e-6,     5216.580372801,  0.642063318 ),
        c(    0.006919e-6,     6681.224853400,  6.018501522 ),
        ## 111, 120 ##
        c(    0.006826e-6,     7632.943259650,  3.458654112 ),
        c(    0.005308e-6,    -1592.596013633,  2.500382359 ),
        c(    0.005096e-6,    11371.704689758,  2.547107806 ),
        c(    0.004841e-6,     5333.900241022,  0.437078094 ),
        c(    0.005582e-6,     5966.683980335,  2.246174308 ),
        c(    0.006304e-6,    11926.254413669,  2.512929171 ),
        c(    0.006603e-6,    23581.258177318,  5.393136889 ),
        c(    0.005123e-6,       -1.484472708,  2.999641028 ),
        c(    0.004648e-6,     1589.072895284,  1.275847090 ),
        c(    0.005119e-6,     6438.496249426,  1.486539246 ),
        ## 121, 130 ##
        c(    0.004521e-6,     4292.330832950,  6.140635794 ),
        c(    0.005680e-6,    23013.539539587,  4.557814849 ),
        c(    0.005488e-6,       -3.455808046,  0.090675389 ),
        c(    0.004193e-6,     7234.794256242,  4.869091389 ),
        c(    0.003742e-6,     7238.675591600,  4.691976180 ),
        c(    0.004148e-6,     -110.206321219,  3.016173439 ),
        c(    0.004553e-6,    11499.656222793,  5.554998314 ),
        c(    0.004892e-6,     5436.993015240,  1.475415597 ),
        c(    0.004044e-6,     4732.030627343,  1.398784824 ),
        c(    0.004164e-6,    12491.370101415,  5.650931916 ),
        ## 131, 140 ##
        c(    0.004349e-6,    11513.883316794,  2.181745369 ),
        c(    0.003919e-6,    12528.018664345,  5.823319737 ),
        c(    0.003129e-6,     6836.645252834,  0.003844094 ),
        c(    0.004080e-6,    -7058.598461315,  3.690360123 ),
        c(    0.003270e-6,       76.266071276,  1.517189902 ),
        c(    0.002954e-6,     6283.143160294,  4.447203799 ),
        c(    0.002872e-6,       28.449187468,  1.158692983 ),
        c(    0.002881e-6,      735.876513532,  0.349250250 ),
        c(    0.003279e-6,     5849.364112115,  4.893384368 ),
        c(    0.003625e-6,     6209.778724132,  1.473760578 ),
        ## 141, 150 ##
        c(    0.003074e-6,      949.175608970,  5.185878737 ),
        c(    0.002775e-6,     9917.696874510,  1.030026325 ),
        c(    0.002646e-6,    10973.555686350,  3.918259169 ),
        c(    0.002575e-6,    25132.303399966,  6.109659023 ),
        c(    0.003500e-6,      263.083923373,  1.892100742 ),
        c(    0.002740e-6,    18319.536584880,  4.320519510 ),
        c(    0.002464e-6,      202.253395174,  4.698203059 ),
        c(    0.002409e-6,        2.542797281,  5.325009315 ),
        c(    0.003354e-6,   -90955.551694697,  1.942656623 ),
        c(    0.002296e-6,     6496.374945429,  5.061810696 ),
        ## 151, 160 ##
        c(    0.003002e-6,     6172.869528772,  2.797822767 ),
        c(    0.003202e-6,    27511.467873537,  0.531673101 ),
        c(    0.002954e-6,    -6283.008539689,  4.533471191 ),
        c(    0.002353e-6,      639.897286314,  3.734548088 ),
        c(    0.002401e-6,    16200.772724501,  2.605547070 ),
        c(    0.003053e-6,   233141.314403759,  3.029030662 ),
        c(    0.003024e-6,    83286.914269554,  2.355556099 ),
        c(    0.002863e-6,    17298.182327326,  5.240963796 ),
        c(    0.002103e-6,    -7079.373856808,  5.756641637 ),
        c(    0.002303e-6,    83996.847317911,  2.013686814 ),
        ## 161, 170 ##
        c(    0.002303e-6,    18073.704938650,  1.089100410 ),
        c(    0.002381e-6,       63.735898303,  0.759188178 ),
        c(    0.002493e-6,     6386.168624210,  0.645026535 ),
        c(    0.002366e-6,        3.932153263,  6.215885448 ),
        c(    0.002169e-6,    11015.106477335,  4.845297676 ),
        c(    0.002397e-6,     6243.458341645,  3.809290043 ),
        c(    0.002183e-6,     1162.474704408,  6.179611691 ),
        c(    0.002353e-6,     6246.427287062,  4.781719760 ),
        c(    0.002199e-6,     -245.831646229,  5.956152284 ),
        c(    0.001729e-6,     3894.181829542,  1.264976635 ),
        ## 171, 180 ##
        c(    0.001896e-6,    -3128.388765096,  4.914231596 ),
        c(    0.002085e-6,       35.164090221,  1.405158503 ),
        c(    0.002024e-6,    14712.317116458,  2.752035928 ),
        c(    0.001737e-6,     6290.189396992,  5.280820144 ),
        c(    0.002229e-6,      491.557929457,  1.571007057 ),
        c(    0.001602e-6,    14314.168113050,  4.203664806 ),
        c(    0.002186e-6,      454.909366527,  1.402101526 ),
        c(    0.001897e-6,    22483.848574493,  4.167932508 ),
        c(    0.001825e-6,    -3738.761430108,  0.545828785 ),
        c(    0.001894e-6,     1052.268383188,  5.817167450 ),
        ## 181, 190 ##
        c(    0.001421e-6,       20.355319399,  2.419886601 ),
        c(    0.001408e-6,    10984.192351700,  2.732084787 ),
        c(    0.001847e-6,    10873.986030480,  2.903477885 ),
        c(    0.001391e-6,    -8635.942003763,  0.593891500 ),
        c(    0.001388e-6,       -7.046236698,  1.166145902 ),
        c(    0.001810e-6,   -88860.057071188,  0.487355242 ),
        c(    0.001288e-6,    -1990.745017041,  3.913022880 ),
        c(    0.001297e-6,    23543.230504682,  3.063805171 ),
        c(    0.001335e-6,     -266.607041722,  3.995764039 ),
        c(    0.001376e-6,    10969.965257698,  5.152914309 ),
        ## 191, 200 ##
        c(    0.001745e-6,   244287.600007027,  3.626395673 ),
        c(    0.001649e-6,    31441.677569757,  1.952049260 ),
        c(    0.001416e-6,     9225.539273283,  4.996408389 ),
        c(    0.001238e-6,     4804.209275927,  5.503379738 ),
        c(    0.001472e-6,     4590.910180489,  4.164913291 ),
        c(    0.001169e-6,     6040.347246017,  5.841719038 ),
        c(    0.001039e-6,     5540.085789459,  2.769753519 ),
        c(    0.001004e-6,     -170.672870619,  0.755008103 ),
        c(    0.001284e-6,    10575.406682942,  5.306538209 ),
        c(    0.001278e-6,       71.812653151,  4.713486491 ),
        ## 201, 210 ##
        c(    0.001321e-6,    18209.330263660,  2.624866359 ),
        c(    0.001297e-6,    21228.392023546,  0.382603541 ),
        c(    0.000954e-6,     6282.095528923,  0.882213514 ),
        c(    0.001145e-6,     6058.731054289,  1.169483931 ),
        c(    0.000979e-6,     5547.199336460,  5.448375984 ),
        c(    0.000987e-6,    -6262.300454499,  2.656486959 ),
        c(    0.001070e-6,  -154717.609887482,  1.827624012 ),
        c(    0.000991e-6,     4701.116501708,  4.387001801 ),
        c(    0.001155e-6,      -14.227094002,  3.042700750 ),
        c(    0.001176e-6,      277.034993741,  3.335519004 ),
        ## 211, 220 ##
        c(    0.000890e-6,    13916.019109642,  5.601498297 ),
        c(    0.000884e-6,    -1551.045222648,  1.088831705 ),
        c(    0.000876e-6,     5017.508371365,  3.969902609 ),
        c(    0.000806e-6,    15110.466119866,  5.142876744 ),
        c(    0.000773e-6,    -4136.910433516,  0.022067765 ),
        c(    0.001077e-6,      175.166059800,  1.844913056 ),
        c(    0.000954e-6,    -6284.056171060,  0.968480906 ),
        c(    0.000737e-6,     5326.786694021,  4.923831588 ),
        c(    0.000845e-6,     -433.711737877,  4.749245231 ),
        c(    0.000819e-6,     8662.240323563,  5.991247817 ),
        ## 221, 230 ##
        c(    0.000852e-6,      199.072001436,  2.189604979 ),
        c(    0.000723e-6,    17256.631536341,  6.068719637 ),
        c(    0.000940e-6,     6037.244203762,  6.197428148 ),
        c(    0.000885e-6,    11712.955318231,  3.280414875 ),
        c(    0.000706e-6,    12559.038152982,  2.824848947 ),
        c(    0.000732e-6,     2379.164473572,  2.501813417 ),
        c(    0.000764e-6,    -6127.655450557,  2.236346329 ),
        c(    0.000908e-6,      131.541961686,  2.521257490 ),
        c(    0.000907e-6,    35371.887265976,  3.370195967 ),
        c(    0.000673e-6,     1066.495477190,  3.876512374 ),
        ## 231, 240 ##
        c(    0.000814e-6,    17654.780539750,  4.627122566 ),
        c(    0.000630e-6,       36.027866677,  0.156368499 ),
        c(    0.000798e-6,      515.463871093,  5.151962502 ),
        c(    0.000798e-6,      148.078724426,  5.909225055 ),
        c(    0.000806e-6,      309.278322656,  6.054064447 ),
        c(    0.000607e-6,      -39.617508346,  2.839021623 ),
        c(    0.000601e-6,      412.371096874,  3.984225404 ),
        c(    0.000646e-6,    11403.676995575,  3.852959484 ),
        c(    0.000704e-6,    13521.751441591,  2.300991267 ),
        c(    0.000603e-6,   -65147.619767937,  4.140083146 ),
        ## 241, 250 ##
        c(    0.000609e-6,    10177.257679534,  0.437122327 ),
        c(    0.000631e-6,     5767.611978898,  4.026532329 ),
        c(    0.000576e-6,    11087.285125918,  4.760293101 ),
        c(    0.000674e-6,    14945.316173554,  6.270510511 ),
        c(    0.000726e-6,     5429.879468239,  6.039606892 ),
        c(    0.000710e-6,    28766.924424484,  5.672617711 ),
        c(    0.000647e-6,    11856.218651625,  3.397132627 ),
        c(    0.000678e-6,    -5481.254918868,  6.249666675 ),
        c(    0.000618e-6,    22003.914634870,  2.466427018 ),
        c(    0.000738e-6,     6134.997125565,  2.242668890 ),
        ## 251, 260 ##
        c(    0.000660e-6,      625.670192312,  5.864091907 ),
        c(    0.000694e-6,     3496.032826134,  2.668309141 ),
        c(    0.000531e-6,     6489.261398429,  1.681888780 ),
        c(    0.000611e-6,  -143571.324284214,  2.424978312 ),
        c(    0.000575e-6,    12043.574281889,  4.216492400 ),
        c(    0.000553e-6,    12416.588502848,  4.772158039 ),
        c(    0.000689e-6,     4686.889407707,  6.224271088 ),
        c(    0.000495e-6,     7342.457780181,  3.817285811 ),
        c(    0.000567e-6,     3634.621024518,  1.649264690 ),
        c(    0.000515e-6,    18635.928454536,  3.945345892 ),
        ## 261, 270 ##
        c(    0.000486e-6,     -323.505416657,  4.061673868 ),
        c(    0.000662e-6,    25158.601719765,  1.794058369 ),
        c(    0.000509e-6,      846.082834751,  3.053874588 ),
        c(    0.000472e-6,   -12569.674818332,  5.112133338 ),
        c(    0.000461e-6,     6179.983075773,  0.513669325 ),
        c(    0.000641e-6,    83467.156352816,  3.210727723 ),
        c(    0.000520e-6,    10344.295065386,  2.445597761 ),
        c(    0.000493e-6,    18422.629359098,  1.676939306 ),
        c(    0.000478e-6,     1265.567478626,  5.487314569 ),
        c(    0.000472e-6,      -18.159247265,  1.999707589 ),
        ## 271, 280 ##
        c(    0.000559e-6,    11190.377900137,  5.783236356 ),
        c(    0.000494e-6,     9623.688276691,  3.022645053 ),
        c(    0.000463e-6,     5739.157790895,  1.411223013 ),
        c(    0.000432e-6,    16858.482532933,  1.179256434 ),
        c(    0.000574e-6,    72140.628666286,  1.758191830 ),
        c(    0.000484e-6,    17267.268201691,  3.290589143 ),
        c(    0.000550e-6,     4907.302050146,  0.864024298 ),
        c(    0.000399e-6,       14.977853527,  2.094441910 ),
        c(    0.000491e-6,      224.344795702,  0.878372791 ),
        c(    0.000432e-6,    20426.571092422,  6.003829241 ),
        ## 281, 290 ##
        c(    0.000481e-6,     5749.452731634,  4.309591964 ),
        c(    0.000480e-6,     5757.317038160,  1.142348571 ),
        c(    0.000485e-6,     6702.560493867,  0.210580917 ),
        c(    0.000426e-6,     6055.549660552,  4.274476529 ),
        c(    0.000480e-6,     5959.570433334,  5.031351030 ),
        c(    0.000466e-6,    12562.628581634,  4.959581597 ),
        c(    0.000520e-6,    39302.096962196,  4.788002889 ),
        c(    0.000458e-6,    12132.439962106,  1.880103788 ),
        c(    0.000470e-6,    12029.347187887,  1.405611197 ),
        c(    0.000416e-6,    -7477.522860216,  1.082356330 ),
        ## 291, 300 ##
        c(    0.000449e-6,    11609.862544012,  4.179989585 ),
        c(    0.000465e-6,    17253.041107690,  0.353496295 ),
        c(    0.000362e-6,    -4535.059436924,  1.583849576 ),
        c(    0.000383e-6,    21954.157609398,  3.747376371 ),
        c(    0.000389e-6,       17.252277143,  1.395753179 ),
        c(    0.000331e-6,    18052.929543158,  0.566790582 ),
        c(    0.000430e-6,    13517.870106233,  0.685827538 ),
        c(    0.000368e-6,    -5756.908003246,  0.731374317 ),
        c(    0.000330e-6,    10557.594160824,  3.710043680 ),
        c(    0.000332e-6,    20199.094959633,  1.652901407 ),
        ## 301, 310 ##
        c(    0.000384e-6,    11933.367960670,  5.827781531 ),
        c(    0.000387e-6,    10454.501386605,  2.541182564 ),
        c(    0.000325e-6,    15671.081759407,  2.178850542 ),
        c(    0.000318e-6,      138.517496871,  2.253253037 ),
        c(    0.000305e-6,     9388.005909415,  0.578340206 ),
        c(    0.000352e-6,     5749.861766548,  3.000297967 ),
        c(    0.000311e-6,     6915.859589305,  1.693574249 ),
        c(    0.000297e-6,    24072.921469776,  1.997249392 ),
        c(    0.000363e-6,     -640.877607382,  5.071820966 ),
        c(    0.000323e-6,    12592.450019783,  1.072262823 ),
        ## 311, 320 ##
        c(    0.000341e-6,    12146.667056108,  4.700657997 ),
        c(    0.000290e-6,     9779.108676125,  1.812320441 ),
        c(    0.000342e-6,     6132.028180148,  4.322238614 ),
        c(    0.000329e-6,     6268.848755990,  3.033827743 ),
        c(    0.000374e-6,    17996.031168222,  3.388716544 ),
        c(    0.000285e-6,     -533.214083444,  4.687313233 ),
        c(    0.000338e-6,     6065.844601290,  0.877776108 ),
        c(    0.000276e-6,       24.298513841,  0.770299429 ),
        c(    0.000336e-6,    -2388.894020449,  5.353796034 ),
        c(    0.000290e-6,     3097.883822726,  4.075291557 ),
        ## 321, 330 ##
        c(    0.000318e-6,      709.933048357,  5.941207518 ),
        c(    0.000271e-6,    13095.842665077,  3.208912203 ),
        c(    0.000331e-6,     6073.708907816,  4.007881169 ),
        c(    0.000292e-6,      742.990060533,  2.714333592 ),
        c(    0.000362e-6,    29088.811415985,  3.215977013 ),
        c(    0.000280e-6,    12359.966151546,  0.710872502 ),
        c(    0.000267e-6,    10440.274292604,  4.730108488 ),
        c(    0.000262e-6,      838.969287750,  1.327720272 ),
        c(    0.000250e-6,    16496.361396202,  0.898769761 ),
        c(    0.000325e-6,    20597.243963041,  0.180044365 ),
        ## 331, 340 ##
        c(    0.000268e-6,     6148.010769956,  5.152666276 ),
        c(    0.000284e-6,     5636.065016677,  5.655385808 ),
        c(    0.000301e-6,     6080.822454817,  2.135396205 ),
        c(    0.000294e-6,     -377.373607916,  3.708784168 ),
        c(    0.000236e-6,     2118.763860378,  1.733578756 ),
        c(    0.000234e-6,     5867.523359379,  5.575209112 ),
        c(    0.000268e-6,  -226858.238553767,  0.069432392 ),
        c(    0.000265e-6,   167283.761587465,  4.369302826 ),
        c(    0.000280e-6,    28237.233459389,  5.304829118 ),
        c(    0.000292e-6,    12345.739057544,  4.096094132 ),
        ## 341, 350 ##
        c(    0.000223e-6,    19800.945956225,  3.069327406 ),
        c(    0.000301e-6,    43232.306658416,  6.205311188 ),
        c(    0.000264e-6,    18875.525869774,  1.417263408 ),
        c(    0.000304e-6,    -1823.175188677,  3.409035232 ),
        c(    0.000301e-6,      109.945688789,  0.510922054 ),
        c(    0.000260e-6,      813.550283960,  2.389438934 ),
        c(    0.000299e-6,   316428.228673312,  5.384595078 ),
        c(    0.000211e-6,     5756.566278634,  3.789392838 ),
        c(    0.000209e-6,     5750.203491159,  1.661943545 ),
        c(    0.000240e-6,    12489.885628707,  5.684549045 ),
        ## 351, 360 ##
        c(    0.000216e-6,     6303.851245484,  3.862942261 ),
        c(    0.000203e-6,     1581.959348283,  5.549853589 ),
        c(    0.000200e-6,     5642.198242609,  1.016115785 ),
        c(    0.000197e-6,      -70.849445304,  4.690702525 ),
        c(    0.000227e-6,     6287.008003254,  2.911891613 ),
        c(    0.000197e-6,      533.623118358,  1.048982898 ),
        c(    0.000205e-6,    -6279.485421340,  1.829362730 ),
        c(    0.000209e-6,   -10988.808157535,  2.636140084 ),
        c(    0.000208e-6,     -227.526189440,  4.127883842 ),
        c(    0.000191e-6,      415.552490612,  4.401165650 ),
        ## 361, 370 ##
        c(    0.000190e-6,    29296.615389579,  4.175658539 ),
        c(    0.000264e-6,    66567.485864652,  4.601102551 ),
        c(    0.000256e-6,    -3646.350377354,  0.506364778 ),
        c(    0.000188e-6,    13119.721102825,  2.032195842 ),
        c(    0.000185e-6,     -209.366942175,  4.694756586 ),
        c(    0.000198e-6,    25934.124331089,  3.832703118 ),
        c(    0.000195e-6,     4061.219215394,  3.308463427 ),
        c(    0.000234e-6,     5113.487598583,  1.716090661 ),
        c(    0.000188e-6,     1478.866574064,  5.686865780 ),
        c(    0.000222e-6,    11823.161639450,  1.942386641 ),
        ## 371, 380 ##
        c(    0.000181e-6,    10770.893256262,  1.999482059 ),
        c(    0.000171e-6,     6546.159773364,  1.182807992 ),
        c(    0.000206e-6,       70.328180442,  5.934076062 ),
        c(    0.000169e-6,    20995.392966449,  2.169080622 ),
        c(    0.000191e-6,    10660.686935042,  5.405515999 ),
        c(    0.000228e-6,    33019.021112205,  4.656985514 ),
        c(    0.000184e-6,    -4933.208440333,  3.327476868 ),
        c(    0.000220e-6,     -135.625325010,  1.765430262 ),
        c(    0.000166e-6,    23141.558382925,  3.454132746 ),
        c(    0.000191e-6,     6144.558353121,  5.020393445 ),
        ## 381, 390 ##
        c(    0.000180e-6,     6084.003848555,  0.602182191 ),
        c(    0.000163e-6,    17782.732072784,  4.960593133 ),
        c(    0.000225e-6,    16460.333529525,  2.596451817 ),
        c(    0.000222e-6,     5905.702242076,  3.731990323 ),
        c(    0.000204e-6,      227.476132789,  5.636192701 ),
        c(    0.000159e-6,    16737.577236597,  3.600691544 ),
        c(    0.000200e-6,     6805.653268085,  0.868220961 ),
        c(    0.000187e-6,    11919.140866668,  2.629456641 ),
        c(    0.000161e-6,      127.471796607,  2.862574720 ),
        c(    0.000205e-6,     6286.666278643,  1.742882331 ),
        ## 391, 400 ##
        c(    0.000189e-6,      153.778810485,  4.812372643 ),
        c(    0.000168e-6,    16723.350142595,  0.027860588 ),
        c(    0.000149e-6,    11720.068865232,  0.659721876 ),
        c(    0.000189e-6,     5237.921013804,  5.245313000 ),
        c(    0.000143e-6,     6709.674040867,  4.317625647 ),
        c(    0.000146e-6,     4487.817406270,  4.815297007 ),
        c(    0.000144e-6,     -664.756045130,  5.381366880 ),
        c(    0.000175e-6,     5127.714692584,  4.728443327 ),
        c(    0.000162e-6,     6254.626662524,  1.435132069 ),
        c(    0.000187e-6,    47162.516354635,  1.354371923 ),
        ## 401, 410 ##
        c(    0.000146e-6,    11080.171578918,  3.369695406 ),
        c(    0.000180e-6,     -348.924420448,  2.490902145 ),
        c(    0.000148e-6,      151.047669843,  3.799109588 ),
        c(    0.000157e-6,     6197.248551160,  1.284375887 ),
        c(    0.000167e-6,      146.594251718,  0.759969109 ),
        c(    0.000133e-6,    -5331.357443741,  5.409701889 ),
        c(    0.000154e-6,       95.979227218,  3.366890614 ),
        c(    0.000148e-6,    -6418.140930027,  3.384104996 ),
        c(    0.000128e-6,    -6525.804453965,  3.803419985 ),
        c(    0.000130e-6,    11293.470674356,  0.939039445 ),
        ## 411, 420 ##
        c(    0.000152e-6,    -5729.506447149,  0.734117523 ),
        c(    0.000138e-6,      210.117701700,  2.564216078 ),
        c(    0.000123e-6,     6066.595360816,  4.517099537 ),
        c(    0.000140e-6,    18451.078546566,  0.642049130 ),
        c(    0.000126e-6,    11300.584221356,  3.485280663 ),
        c(    0.000119e-6,    10027.903195729,  3.217431161 ),
        c(    0.000151e-6,     4274.518310832,  4.404359108 ),
        c(    0.000117e-6,     6072.958148291,  0.366324650 ),
        c(    0.000165e-6,    -7668.637425143,  4.298212528 ),
        c(    0.000117e-6,    -6245.048177356,  5.379518958 ),
        ## 421, 430 ##
        c(    0.000130e-6,    -5888.449964932,  4.527681115 ),
        c(    0.000121e-6,     -543.918059096,  6.109429504 ),
        c(    0.000162e-6,     9683.594581116,  5.720092446 ),
        c(    0.000141e-6,     6219.339951688,  0.679068671 ),
        c(    0.000118e-6,    22743.409379516,  4.881123092 ),
        c(    0.000129e-6,     1692.165669502,  0.351407289 ),
        c(    0.000126e-6,     5657.405657679,  5.146592349 ),
        c(    0.000114e-6,      728.762966531,  0.520791814 ),
        c(    0.000120e-6,       52.596639600,  0.948516300 ),
        c(    0.000115e-6,       65.220371012,  3.504914846 ),
        ## 431, 440 ##
        c(    0.000126e-6,     5881.403728234,  5.577502482 ),
        c(    0.000158e-6,   163096.180360983,  2.957128968 ),
        c(    0.000134e-6,    12341.806904281,  2.598576764 ),
        c(    0.000151e-6,    16627.370915377,  3.985702050 ),
        c(    0.000109e-6,     1368.660252845,  0.014730471 ),
        c(    0.000131e-6,     6211.263196841,  0.085077024 ),
        c(    0.000146e-6,     5792.741760812,  0.708426604 ),
        c(    0.000146e-6,      -77.750543984,  3.121576600 ),
        c(    0.000107e-6,     5341.013788022,  0.288231904 ),
        c(    0.000138e-6,     6281.591377283,  2.797450317 ),
        ## 441, 450 ##
        c(    0.000113e-6,    -6277.552925684,  2.788904128 ),
        c(    0.000115e-6,     -525.758811831,  5.895222200 ),
        c(    0.000138e-6,     6016.468808270,  6.096188999 ),
        c(    0.000139e-6,    23539.707386333,  2.028195445 ),
        c(    0.000146e-6,    -4176.041342449,  4.660008502 ),
        c(    0.000107e-6,    16062.184526117,  4.066520001 ),
        c(    0.000142e-6,    83783.548222473,  2.936315115 ),
        c(    0.000128e-6,     9380.959672717,  3.223844306 ),
        c(    0.000135e-6,     6205.325306007,  1.638054048 ),
        c(    0.000101e-6,     2699.734819318,  5.481603249 ),
        ## 451, 460 ##
        c(    0.000104e-6,     -568.821874027,  2.205734493 ),
        c(    0.000103e-6,     6321.103522627,  2.440421099 ),
        c(    0.000119e-6,     6321.208885629,  2.547496264 ),
        c(    0.000138e-6,     1975.492545856,  2.314608466 ),
        c(    0.000121e-6,      137.033024162,  4.539108237 ),
        c(    0.000123e-6,    19402.796952817,  4.538074405 ),
        c(    0.000119e-6,    22805.735565994,  2.869040566 ),
        c(    0.000133e-6,    64471.991241142,  6.056405489 ),
        c(    0.000129e-6,      -85.827298831,  2.540635083 ),
        c(    0.000131e-6,    13613.804277336,  4.005732868 ),
        ## 461, 470 ##
        c(    0.000104e-6,     9814.604100291,  1.959967212 ),
        c(    0.000112e-6,    16097.679950283,  3.589026260 ),
        c(    0.000123e-6,     2107.034507542,  1.728627253 ),
        c(    0.000121e-6,    36949.230808424,  6.072332087 ),
        c(    0.000108e-6,   -12539.853380183,  3.716133846 ),
        c(    0.000113e-6,    -7875.671863624,  2.725771122 ),
        c(    0.000109e-6,     4171.425536614,  4.033338079 ),
        c(    0.000101e-6,     6247.911759770,  3.441347021 ),
        c(    0.000113e-6,     7330.728427345,  0.656372122 ),
        c(    0.000113e-6,    51092.726050855,  2.791483066 ),
        ## 471, 480 ##
        c(    0.000106e-6,     5621.842923210,  1.815323326 ),
        c(    0.000101e-6,      111.430161497,  5.711033677 ),
        c(    0.000103e-6,      909.818733055,  2.812745443 ),
        c(    0.000101e-6,     1790.642637886,  1.965746028 ),

        ## T ##
        c(  102.156724e-6,     6283.075849991,  4.249032005 ),
        c(    1.706807e-6,    12566.151699983,  4.205904248 ),
        c(    0.269668e-6,      213.299095438,  3.400290479 ),
        c(    0.265919e-6,      529.690965095,  5.836047367 ),
        c(    0.210568e-6,       -3.523118349,  6.262738348 ),
        c(    0.077996e-6,     5223.693919802,  4.670344204 ),
        ## 481, 490 ##
        c(    0.054764e-6,     1577.343542448,  4.534800170 ),
        c(    0.059146e-6,       26.298319800,  1.083044735 ),
        c(    0.034420e-6,     -398.149003408,  5.980077351 ),
        c(    0.032088e-6,    18849.227549974,  4.162913471 ),
        c(    0.033595e-6,     5507.553238667,  5.980162321 ),
        c(    0.029198e-6,     5856.477659115,  0.623811863 ),
        c(    0.027764e-6,      155.420399434,  3.745318113 ),
        c(    0.025190e-6,     5746.271337896,  2.980330535 ),
        c(    0.022997e-6,     -796.298006816,  1.174411803 ),
        c(    0.024976e-6,     5760.498431898,  2.467913690 ),
        ## 491, 500 ##
        c(    0.021774e-6,      206.185548437,  3.854787540 ),
        c(    0.017925e-6,     -775.522611324,  1.092065955 ),
        c(    0.013794e-6,      426.598190876,  2.699831988 ),
        c(    0.013276e-6,     6062.663207553,  5.845801920 ),
        c(    0.011774e-6,    12036.460734888,  2.292832062 ),
        c(    0.012869e-6,     6076.890301554,  5.333425680 ),
        c(    0.012152e-6,     1059.381930189,  6.222874454 ),
        c(    0.011081e-6,       -7.113547001,  5.154724984 ),
        c(    0.010143e-6,     4694.002954708,  4.044013795 ),
        c(    0.009357e-6,     5486.777843175,  3.416081409 ),
        ## 501, 510 ##
        c(    0.010084e-6,      522.577418094,  0.749320262 ),
        c(    0.008587e-6,    10977.078804699,  2.777152598 ),
        c(    0.008628e-6,     6275.962302991,  4.562060226 ),
        c(    0.008158e-6,     -220.412642439,  5.806891533 ),
        c(    0.007746e-6,     2544.314419883,  1.603197066 ),
        c(    0.007670e-6,     2146.165416475,  3.000200440 ),
        c(    0.007098e-6,       74.781598567,  0.443725817 ),
        c(    0.006180e-6,     -536.804512095,  1.302642751 ),
        c(    0.005818e-6,     5088.628839767,  4.827723531 ),
        c(    0.004945e-6,    -6286.598968340,  0.268305170 ),
        ## 511, 520 ##
        c(    0.004774e-6,     1349.867409659,  5.808636673 ),
        c(    0.004687e-6,     -242.728603974,  5.154890570 ),
        c(    0.006089e-6,     1748.016413067,  4.403765209 ),
        c(    0.005975e-6,    -1194.447010225,  2.583472591 ),
        c(    0.004229e-6,      951.718406251,  0.931172179 ),
        c(    0.005264e-6,      553.569402842,  2.336107252 ),
        c(    0.003049e-6,     5643.178563677,  1.362634430 ),
        c(    0.002974e-6,     6812.766815086,  1.583012668 ),
        c(    0.003403e-6,    -2352.866153772,  2.552189886 ),
        c(    0.003030e-6,      419.484643875,  5.286473844 ),
        ## 521, 530 ##
        c(    0.003210e-6,       -7.046236698,  1.863796539 ),
        c(    0.003058e-6,     9437.762934887,  4.226420633 ),
        c(    0.002589e-6,    12352.852604545,  1.991935820 ),
        c(    0.002927e-6,     5216.580372801,  2.319951253 ),
        c(    0.002425e-6,     5230.807466803,  3.084752833 ),
        c(    0.002656e-6,     3154.687084896,  2.487447866 ),
        c(    0.002445e-6,    10447.387839604,  2.347139160 ),
        c(    0.002990e-6,     4690.479836359,  6.235872050 ),
        c(    0.002890e-6,     5863.591206116,  0.095197563 ),
        c(    0.002498e-6,     6438.496249426,  2.994779800 ),
        ## 531, 540 ##
        c(    0.001889e-6,     8031.092263058,  3.569003717 ),
        c(    0.002567e-6,      801.820931124,  3.425611498 ),
        c(    0.001803e-6,   -71430.695617928,  2.192295512 ),
        c(    0.001782e-6,        3.932153263,  5.180433689 ),
        c(    0.001694e-6,    -4705.732307544,  4.641779174 ),
        c(    0.001704e-6,    -1592.596013633,  3.997097652 ),
        c(    0.001735e-6,     5849.364112115,  0.417558428 ),
        c(    0.001643e-6,     8429.241266467,  2.180619584 ),
        c(    0.001680e-6,       38.133035638,  4.164529426 ),
        c(    0.002045e-6,     7084.896781115,  0.526323854 ),
        ## 541, 550 ##
        c(    0.001458e-6,     4292.330832950,  1.356098141 ),
        c(    0.001437e-6,       20.355319399,  3.895439360 ),
        c(    0.001738e-6,     6279.552731642,  0.087484036 ),
        c(    0.001367e-6,    14143.495242431,  3.987576591 ),
        c(    0.001344e-6,     7234.794256242,  0.090454338 ),
        c(    0.001438e-6,    11499.656222793,  0.974387904 ),
        c(    0.001257e-6,     6836.645252834,  1.509069366 ),
        c(    0.001358e-6,    11513.883316794,  0.495572260 ),
        c(    0.001628e-6,     7632.943259650,  4.968445721 ),
        c(    0.001169e-6,      103.092774219,  2.838496795 ),
        ## 551, 560 ##
        c(    0.001162e-6,     4164.311989613,  3.408387778 ),
        c(    0.001092e-6,     6069.776754553,  3.617942651 ),
        c(    0.001008e-6,    17789.845619785,  0.286350174 ),
        c(    0.001008e-6,      639.897286314,  1.610762073 ),
        c(    0.000918e-6,    10213.285546211,  5.532798067 ),
        c(    0.001011e-6,    -6256.777530192,  0.661826484 ),
        c(    0.000753e-6,    16730.463689596,  3.905030235 ),
        c(    0.000737e-6,    11926.254413669,  4.641956361 ),
        c(    0.000694e-6,     3340.612426700,  2.111120332 ),
        c(    0.000701e-6,     3894.181829542,  2.760823491 ),
        ## 561, 570 ##
        c(    0.000689e-6,     -135.065080035,  4.768800780 ),
        c(    0.000700e-6,    13367.972631107,  5.760439898 ),
        c(    0.000664e-6,     6040.347246017,  1.051215840 ),
        c(    0.000654e-6,     5650.292110678,  4.911332503 ),
        c(    0.000788e-6,     6681.224853400,  4.699648011 ),
        c(    0.000628e-6,     5333.900241022,  5.024608847 ),
        c(    0.000755e-6,     -110.206321219,  4.370971253 ),
        c(    0.000628e-6,     6290.189396992,  3.660478857 ),
        c(    0.000635e-6,    25132.303399966,  4.121051532 ),
        c(    0.000534e-6,     5966.683980335,  1.173284524 ),
        ## 571, 580 ##
        c(    0.000543e-6,     -433.711737877,  0.345585464 ),
        c(    0.000517e-6,    -1990.745017041,  5.414571768 ),
        c(    0.000504e-6,     5767.611978898,  2.328281115 ),
        c(    0.000485e-6,     5753.384884897,  1.685874771 ),
        c(    0.000463e-6,     7860.419392439,  5.297703006 ),
        c(    0.000604e-6,      515.463871093,  0.591998446 ),
        c(    0.000443e-6,    12168.002696575,  4.830881244 ),
        c(    0.000570e-6,      199.072001436,  3.899190272 ),
        c(    0.000465e-6,    10969.965257698,  0.476681802 ),
        c(    0.000424e-6,    -7079.373856808,  1.112242763 ),
        ## 581, 590 ##
        c(    0.000427e-6,      735.876513532,  1.994214480 ),
        c(    0.000478e-6,    -6127.655450557,  3.778025483 ),
        c(    0.000414e-6,    10973.555686350,  5.441088327 ),
        c(    0.000512e-6,     1589.072895284,  0.107123853 ),
        c(    0.000378e-6,    10984.192351700,  0.915087231 ),
        c(    0.000402e-6,    11371.704689758,  4.107281715 ),
        c(    0.000453e-6,     9917.696874510,  1.917490952 ),
        c(    0.000395e-6,      149.563197135,  2.763124165 ),
        c(    0.000371e-6,     5739.157790895,  3.112111866 ),
        c(    0.000350e-6,    11790.629088659,  0.440639857 ),
        ## 591, 600 ##
        c(    0.000356e-6,     6133.512652857,  5.444568842 ),
        c(    0.000344e-6,      412.371096874,  5.676832684 ),
        c(    0.000383e-6,      955.599741609,  5.559734846 ),
        c(    0.000333e-6,     6496.374945429,  0.261537984 ),
        c(    0.000340e-6,     6055.549660552,  5.975534987 ),
        c(    0.000334e-6,     1066.495477190,  2.335063907 ),
        c(    0.000399e-6,    11506.769769794,  5.321230910 ),
        c(    0.000314e-6,    18319.536584880,  2.313312404 ),
        c(    0.000424e-6,     1052.268383188,  1.211961766 ),
        c(    0.000307e-6,       63.735898303,  3.169551388 ),
        ## 601, 610 ##
        c(    0.000329e-6,       29.821438149,  6.106912080 ),
        c(    0.000357e-6,     6309.374169791,  4.223760346 ),
        c(    0.000312e-6,    -3738.761430108,  2.180556645 ),
        c(    0.000301e-6,      309.278322656,  1.499984572 ),
        c(    0.000268e-6,    12043.574281889,  2.447520648 ),
        c(    0.000257e-6,    12491.370101415,  3.662331761 ),
        c(    0.000290e-6,      625.670192312,  1.272834584 ),
        c(    0.000256e-6,     5429.879468239,  1.913426912 ),
        c(    0.000339e-6,     3496.032826134,  4.165930011 ),
        c(    0.000283e-6,     3930.209696220,  4.325565754 ),
        ## 611, 620 ##
        c(    0.000241e-6,    12528.018664345,  3.832324536 ),
        c(    0.000304e-6,     4686.889407707,  1.612348468 ),
        c(    0.000259e-6,    16200.772724501,  3.470173146 ),
        c(    0.000238e-6,    12139.553509107,  1.147977842 ),
        c(    0.000236e-6,     6172.869528772,  3.776271728 ),
        c(    0.000296e-6,    -7058.598461315,  0.460368852 ),
        c(    0.000306e-6,    10575.406682942,  0.554749016 ),
        c(    0.000251e-6,    17298.182327326,  0.834332510 ),
        c(    0.000290e-6,     4732.030627343,  4.759564091 ),
        c(    0.000261e-6,     5884.926846583,  0.298259862 ),
        ## 621, 630 ##
        c(    0.000249e-6,     5547.199336460,  3.749366406 ),
        c(    0.000213e-6,    11712.955318231,  5.415666119 ),
        c(    0.000223e-6,     4701.116501708,  2.703203558 ),
        c(    0.000268e-6,     -640.877607382,  0.283670793 ),
        c(    0.000209e-6,     5636.065016677,  1.238477199 ),
        c(    0.000193e-6,    10177.257679534,  1.943251340 ),
        c(    0.000182e-6,     6283.143160294,  2.456157599 ),
        c(    0.000184e-6,     -227.526189440,  5.888038582 ),
        c(    0.000182e-6,    -6283.008539689,  0.241332086 ),
        c(    0.000228e-6,    -6284.056171060,  2.657323816 ),
        ## 631, 640 ##
        c(    0.000166e-6,     7238.675591600,  5.930629110 ),
        c(    0.000167e-6,     3097.883822726,  5.570955333 ),
        c(    0.000159e-6,     -323.505416657,  5.786670700 ),
        c(    0.000154e-6,    -4136.910433516,  1.517805532 ),
        c(    0.000176e-6,    12029.347187887,  3.139266834 ),
        c(    0.000167e-6,    12132.439962106,  3.556352289 ),
        c(    0.000153e-6,      202.253395174,  1.463313961 ),
        c(    0.000157e-6,    17267.268201691,  1.586837396 ),
        c(    0.000142e-6,    83996.847317911,  0.022670115 ),
        c(    0.000152e-6,    17260.154654690,  0.708528947 ),
        ## 641, 650 ##
        c(    0.000144e-6,     6084.003848555,  5.187075177 ),
        c(    0.000135e-6,     5756.566278634,  1.993229262 ),
        c(    0.000134e-6,     5750.203491159,  3.457197134 ),
        c(    0.000144e-6,     5326.786694021,  6.066193291 ),
        c(    0.000160e-6,    11015.106477335,  1.710431974 ),
        c(    0.000133e-6,     3634.621024518,  2.836451652 ),
        c(    0.000134e-6,    18073.704938650,  5.453106665 ),
        c(    0.000134e-6,     1162.474704408,  5.326898811 ),
        c(    0.000128e-6,     5642.198242609,  2.511652591 ),
        c(    0.000160e-6,      632.783739313,  5.628785365 ),
        ## 651, 660 ##
        c(    0.000132e-6,    13916.019109642,  0.819294053 ),
        c(    0.000122e-6,    14314.168113050,  5.677408071 ),
        c(    0.000125e-6,    12359.966151546,  5.251984735 ),
        c(    0.000121e-6,     5749.452731634,  2.210924603 ),
        c(    0.000136e-6,     -245.831646229,  1.646502367 ),
        c(    0.000120e-6,     5757.317038160,  3.240883049 ),
        c(    0.000134e-6,    12146.667056108,  3.059480037 ),
        c(    0.000137e-6,     6206.809778716,  1.867105418 ),
        c(    0.000141e-6,    17253.041107690,  2.069217456 ),
        c(    0.000129e-6,    -7477.522860216,  2.781469314 ),
        ## 661, 670 ##
        c(    0.000116e-6,     5540.085789459,  4.281176991 ),
        c(    0.000116e-6,     9779.108676125,  3.320925381 ),
        c(    0.000129e-6,     5237.921013804,  3.497704076 ),
        c(    0.000113e-6,     5959.570433334,  0.983210840 ),
        c(    0.000122e-6,     6282.095528923,  2.674938860 ),
        c(    0.000140e-6,      -11.045700264,  4.957936982 ),
        c(    0.000108e-6,    23543.230504682,  1.390113589 ),
        c(    0.000106e-6,   -12569.674818332,  0.429631317 ),
        c(    0.000110e-6,     -266.607041722,  5.501340197 ),
        c(    0.000115e-6,    12559.038152982,  4.691456618 ),
        ## 671, 680 ##
        c(    0.000134e-6,    -2388.894020449,  0.577313584 ),
        c(    0.000109e-6,    10440.274292604,  6.218148717 ),
        c(    0.000102e-6,     -543.918059096,  1.477842615 ),
        c(    0.000108e-6,    21228.392023546,  2.237753948 ),
        c(    0.000101e-6,    -4535.059436924,  3.100492232 ),
        c(    0.000103e-6,       76.266071276,  5.594294322 ),
        c(    0.000104e-6,      949.175608970,  5.674287810 ),
        c(    0.000101e-6,    13517.870106233,  2.196632348 ),
        c(    0.000100e-6,    11933.367960670,  4.056084160 ),

        ## T^2 ##
        c(    4.322990e-6,     6283.075849991,  2.642893748 ),
        ## 681, 690 ##
        c(    0.406495e-6,        0.000000000,  4.712388980 ),
        c(    0.122605e-6,    12566.151699983,  2.438140634 ),
        c(    0.019476e-6,      213.299095438,  1.642186981 ),
        c(    0.016916e-6,      529.690965095,  4.510959344 ),
        c(    0.013374e-6,       -3.523118349,  1.502210314 ),
        c(    0.008042e-6,       26.298319800,  0.478549024 ),
        c(    0.007824e-6,      155.420399434,  5.254710405 ),
        c(    0.004894e-6,     5746.271337896,  4.683210850 ),
        c(    0.004875e-6,     5760.498431898,  0.759507698 ),
        c(    0.004416e-6,     5223.693919802,  6.028853166 ),
        ## 691, 700 ##
        c(    0.004088e-6,       -7.113547001,  0.060926389 ),
        c(    0.004433e-6,    77713.771467920,  3.627734103 ),
        c(    0.003277e-6,    18849.227549974,  2.327912542 ),
        c(    0.002703e-6,     6062.663207553,  1.271941729 ),
        c(    0.003435e-6,     -775.522611324,  0.747446224 ),
        c(    0.002618e-6,     6076.890301554,  3.633715689 ),
        c(    0.003146e-6,      206.185548437,  5.647874613 ),
        c(    0.002544e-6,     1577.343542448,  6.232904270 ),
        c(    0.002218e-6,     -220.412642439,  1.309509946 ),
        c(    0.002197e-6,     5856.477659115,  2.407212349 ),
        ## 701, 710 ##
        c(    0.002897e-6,     5753.384884897,  5.863842246 ),
        c(    0.001766e-6,      426.598190876,  0.754113147 ),
        c(    0.001738e-6,     -796.298006816,  2.714942671 ),
        c(    0.001695e-6,      522.577418094,  2.629369842 ),
        c(    0.001584e-6,     5507.553238667,  1.341138229 ),
        c(    0.001503e-6,     -242.728603974,  0.377699736 ),
        c(    0.001552e-6,     -536.804512095,  2.904684667 ),
        c(    0.001370e-6,     -398.149003408,  1.265599125 ),
        c(    0.001889e-6,    -5573.142801634,  4.413514859 ),
        c(    0.001722e-6,     6069.776754553,  2.445966339 ),
        ## 711, 720 ##
        c(    0.001124e-6,     1059.381930189,  5.041799657 ),
        c(    0.001258e-6,      553.569402842,  3.849557278 ),
        c(    0.000831e-6,      951.718406251,  2.471094709 ),
        c(    0.000767e-6,     4694.002954708,  5.363125422 ),
        c(    0.000756e-6,     1349.867409659,  1.046195744 ),
        c(    0.000775e-6,      -11.045700264,  0.245548001 ),
        c(    0.000597e-6,     2146.165416475,  4.543268798 ),
        c(    0.000568e-6,     5216.580372801,  4.178853144 ),
        c(    0.000711e-6,     1748.016413067,  5.934271972 ),
        c(    0.000499e-6,    12036.460734888,  0.624434410 ),
        ## 721, 730 ##
        c(    0.000671e-6,    -1194.447010225,  4.136047594 ),
        c(    0.000488e-6,     5849.364112115,  2.209679987 ),
        c(    0.000621e-6,     6438.496249426,  4.518860804 ),
        c(    0.000495e-6,    -6286.598968340,  1.868201275 ),
        c(    0.000456e-6,     5230.807466803,  1.271231591 ),
        c(    0.000451e-6,     5088.628839767,  0.084060889 ),
        c(    0.000435e-6,     5643.178563677,  3.324456609 ),
        c(    0.000387e-6,    10977.078804699,  4.052488477 ),
        c(    0.000547e-6,   161000.685737473,  2.841633844 ),
        c(    0.000522e-6,     3154.687084896,  2.171979966 ),
        ## 731, 740 ##
        c(    0.000375e-6,     5486.777843175,  4.983027306 ),
        c(    0.000421e-6,     5863.591206116,  4.546432249 ),
        c(    0.000439e-6,     7084.896781115,  0.522967921 ),
        c(    0.000309e-6,     2544.314419883,  3.172606705 ),
        c(    0.000347e-6,     4690.479836359,  1.479586566 ),
        c(    0.000317e-6,      801.820931124,  3.553088096 ),
        c(    0.000262e-6,      419.484643875,  0.606635550 ),
        c(    0.000248e-6,     6836.645252834,  3.014082064 ),
        c(    0.000245e-6,    -1592.596013633,  5.519526220 ),
        c(    0.000225e-6,     4292.330832950,  2.877956536 ),
        ## 741, 750 ##
        c(    0.000214e-6,     7234.794256242,  1.605227587 ),
        c(    0.000205e-6,     5767.611978898,  0.625804796 ),
        c(    0.000180e-6,    10447.387839604,  3.499954526 ),
        c(    0.000229e-6,      199.072001436,  5.632304604 ),
        c(    0.000214e-6,      639.897286314,  5.960227667 ),
        c(    0.000175e-6,     -433.711737877,  2.162417992 ),
        c(    0.000209e-6,      515.463871093,  2.322150893 ),
        c(    0.000173e-6,     6040.347246017,  2.556183691 ),
        c(    0.000184e-6,     6309.374169791,  4.732296790 ),
        c(    0.000227e-6,   149854.400134205,  5.385812217 ),
        ## 751, 760 ##
        c(    0.000154e-6,     8031.092263058,  5.120720920 ),
        c(    0.000151e-6,     5739.157790895,  4.815000443 ),
        c(    0.000197e-6,     7632.943259650,  0.222827271 ),
        c(    0.000197e-6,       74.781598567,  3.910456770 ),
        c(    0.000138e-6,     6055.549660552,  1.397484253 ),
        c(    0.000149e-6,    -6127.655450557,  5.333727496 ),
        c(    0.000137e-6,     3894.181829542,  4.281749907 ),
        c(    0.000135e-6,     9437.762934887,  5.979971885 ),
        c(    0.000139e-6,    -2352.866153772,  4.715630782 ),
        c(    0.000142e-6,     6812.766815086,  0.513330157 ),
        ## 761, 770 ##
        c(    0.000120e-6,    -4705.732307544,  0.194160689 ),
        c(    0.000131e-6,   -71430.695617928,  0.000379226 ),
        c(    0.000124e-6,     6279.552731642,  2.122264908 ),
        c(    0.000108e-6,    -6256.777530192,  0.883445696 ),

        ## T^3 ##
        c(    0.143388e-6,     6283.075849991,  1.131453581 ),
        c(    0.006671e-6,    12566.151699983,  0.775148887 ),
        c(    0.001480e-6,      155.420399434,  0.480016880 ),
        c(    0.000934e-6,      213.299095438,  6.144453084 ),
        c(    0.000795e-6,      529.690965095,  2.941595619 ),
        c(    0.000673e-6,     5746.271337896,  0.120415406 ),
        ## 771, 780 ##
        c(    0.000672e-6,     5760.498431898,  5.317009738 ),
        c(    0.000389e-6,     -220.412642439,  3.090323467 ),
        c(    0.000373e-6,     6062.663207553,  3.003551964 ),
        c(    0.000360e-6,     6076.890301554,  1.918913041 ),
        c(    0.000316e-6,      -21.340641002,  5.545798121 ),
        c(    0.000315e-6,     -242.728603974,  1.884932563 ),
        c(    0.000278e-6,      206.185548437,  1.266254859 ),
        c(    0.000238e-6,     -536.804512095,  4.532664830 ),
        c(    0.000185e-6,      522.577418094,  4.578313856 ),
        c(    0.000245e-6,    18849.227549974,  0.587467082 ),
        ## 781, 787 ##
        c(    0.000180e-6,      426.598190876,  5.151178553 ),
        c(    0.000200e-6,      553.569402842,  5.355983739 ),
        c(    0.000141e-6,     5223.693919802,  1.336556009 ),
        c(    0.000104e-6,     5856.477659115,  4.239842759 ),

        ## T^4 ##
        c(    0.003826e-6,     6283.075849991,  5.705257275 ),
        c(    0.000303e-6,    12566.151699983,  5.407132842 ),
        c(    0.000209e-6,      155.420399434,  1.989815753 ))


    ## Time since J2000.0 in Julian millennia. ##
    t  <-  ((date1 - DJ00) + date2) / DJM
    t[t==0] <- 1e-20

    DD2R <- pi/180
    ## ================= ##
    ## Topocentric terms ##
    ## ================= ##

    ## Convert UT to local solar time in radians. ##
    tsol = (ut%%1)*pi + elong

    ## FUNDAMENTAL ARGUMENTS:  Simon et al. 1994. ##

    ## Combine time argument (millennia) with deg/arcsec factor. ##
    w = t / 3600.0

    ## Sun Mean Longitude. ##
    elsun = ((280.46645683 + 1296027711.03429 * w)%%360.0) * DD2R

    ## Sun Mean Anomaly. ##
    emsun = ((357.52910918 + 1295965810.481 * w)%%360.0)*DD2R

    ## Mean Elongation of Moon from Sun. ##
    d = ((297.85019547 + 16029616012.090 * w)%%360.0) * DD2R

    ## Mean Longitude of Jupiter. ##
    elj = ((34.35151874 + 109306899.89453 * w)%%360.0) * DD2R

    ## Mean Longitude of Saturn. ##
    els = ((50.07744430 + 44046398.47038 * w)%%360.0) * DD2R

    ## TOPOCENTRIC TERMS:  Moyer 1981 and Murray 1983. ##
    wt  <-  0.00029e-10 * u * sin(tsol + elsun - els)
    +  0.00100e-10 * u * sin(tsol - 2.0 * emsun)
    +  0.00133e-10 * u * sin(tsol - d)
    +  0.00133e-10 * u * sin(tsol + elsun - elj)
    -  0.00229e-10 * u * sin(tsol + 2.0 * elsun + emsun)
    -  0.02200e-10 * v * cos(elsun + emsun)
    +  0.05312e-10 * u * sin(tsol - emsun)
    -  0.13677e-10 * u * sin(tsol + 2.0 * elsun)
    -  1.31840e-10 * v * cos(elsun)
    +  3.17679e-10 * u * sin(tsol)
    ## ===================== ##
    ## Fairhead et al. model ##
    ## ===================== ##
    ## T**0 ##
    w0 = 0
    dw0 <- 0
    for (j in 474:1){
        w0  <- w0+ fairhd[j,1] * sin(fairhd[j,2] * t + fairhd[j,3])
        dw0 <- dw0+fairhd[j,1] * cos(fairhd[j,2] * t + fairhd[j,3])*fairhd[j,2]
    }

    ## T**1 ##
    w1 = 0
    dw1 <- 0
    for(j in 679:475){
        w1  <- w1+ fairhd[j,1] * sin(fairhd[j,2] * t + fairhd[j,3])
        dw1 <- dw1+fairhd[j,1] * cos(fairhd[j,2] * t + fairhd[j,3])*fairhd[j,2]
    }

    ## T**2 ##
    w2 = 0
    dw2 <- 0
    for(j in 764:680){
        w2  <- w2+ fairhd[j,1] * sin(fairhd[j,2] * t + fairhd[j,3])
        dw2 <- dw2+fairhd[j,1] * cos(fairhd[j,2] * t + fairhd[j,3])*fairhd[j,2]
    }

    ## T**3 ##
    w3 = 0
    dw3 <- 0
    for(j in 784:765){
        w3  <- w3+ fairhd[j,1] * sin(fairhd[j,2] * t + fairhd[j,3])
        dw3 <- dw3+fairhd[j,1] * cos(fairhd[j,2] * t + fairhd[j,3])*fairhd[j,2]
    }

    ## T**4 ##
    w4 = 0
    dw4 <- 0
    for(j in 787:785){
        w4  <- w4+ fairhd[j,1] * sin(fairhd[j,2] * t + fairhd[j,3])
        dw4 <- dw4+fairhd[j,1] * cos(fairhd[j,2] * t + fairhd[j,3])*fairhd[j,2]
    }

    ## Multiply by powers of T and combine. ##
    wf = t * (t * (t * (t * w4 + w3) + w2) + w1) + w0
    dwf = t * (t * (t * (t * dw4 + dw3) + dw2) + dw1)+w1+2*w2*t+3*w3*t^2+4*w4*t^3

    ## Adjustments to use JPL planetary masses instead of IAU. ##
    wj =   0.00065e-6 * sin(6069.776754 * t + 4.021194) +
        0.00033e-6 * sin( 213.299095 * t + 5.543132) +
        (-0.00196e-6 * sin(6208.294251 * t + 5.696701)) +
        (-0.00173e-6 * sin(  74.781599 * t + 2.435900)) +
        0.03638e-6 * t * t
    dwj <-  0.00065e-6 * 6069.776754*cos(6069.776754 * t + 4.021194) +
        0.00033e-6 * 213.299095*cos( 213.299095 * t + 5.543132) +
        (-0.00196e-6 * 6208.294251*cos(6208.294251 * t + 5.696701)) +
        (-0.00173e-6 * 74.781599*cos(  74.781599 * t + 2.435900)) +
        2*0.03638e-6 * t

    ## ============ ##
    ## Final result ##
    ## ============ ##

    ## TDB-TT in seconds. ##
    w  <-  wt + wf + wj#s
#    dw <- (dwf+dwj)/DAYSEC/DJM#s/s; ignore dwt which is difficult to derive analytic form
    dw <- dwf/DAYSEC/DJM#s/s; ignore dwt which is difficult to derive analytic form
    return(cbind(w,dw))
}

sofa_Dtdb2 <- function(jd){
####################################
## NAME:
##   TDB2TDT
##
## AUTHOR:
##   Craig B. Markwardt, NASA/GSFC Code 662, Greenbelt, MD 20770
##   craigm@lheamail.gsfc.nasa.gov
##   UPDATED VERSIONs can be found on my WEB PAGE:
##      http://cow.physics.wisc.edu/~craigm/idl/idl.html
##
## PURPOSE:
##   Relativistic clock corrections due to Earth motion in solar system
##
## MAJOR TOPICS:
##   Planetary Orbits
##
## CALLING SEQUENCE:
##   corr = TDB2TDT(JD, TBASE=, DERIV=deriv)
##
## DESCRIPTION:
##
##   The function TDB2TDT computes relativistic corrections that must
##   be applied when performing high precision absolute timing in the
##   solar system.
##
##   According to general relativity, moving clocks, and clocks at
##   different gravitational potentials, will run at different rates
##   with respect to each other.  A clock placed on the earth will run
##   at a time-variable rate because of the non-constant influence of
##   the sun and other planets.  Thus, for the most demanding
##   astrophysical timing applications -- high precision pulsar timing
##   -- times in the accelerating earth observer's frame must be
##   corrected to an inertial frame, such as the solar system
##   barycenter (SSB).  This correction is also convenient because the
##   coordinate time at the SSB is the ephemeris time of the JPL
##   Planetary Ephemeris.
##
##   In general, the difference in the rate of Ti, the time kept by an
##   arbitrary clock, and the rate of T, the ephemeris time, is given
##   by the expression (Standish 1998):
##
##      dTi/dT = 1 - (Ui + vi^2/2) / c^2
##
##   where Ui is the potential of clock i, and vi is the velocity of
##   clock i.  However, when integrated, this expression depends on the
##   position of an individual clock.  A more convenient approximate
##   expression is:
##
##     T = Ti + (robs(Ti) . vearth(T))/c^2 + dtgeo(Ti) + TDB2TDT(Ti)
##
##   where robs is the vector from the geocenter to the observer##
##   vearth is the vector velocity of the earth## and dtgeo is a
##   correction to convert from the observer's clock to geocentric TT
##   time.  TDB2TDT is the value computed by this function, the
##   correction to convert from the geocenter to the solar system
##   barycenter.
##
##   As the above equation shows, while this function provides an
##   important component of the correction, the user must also be
##   responsible for (a) correcting their times to the geocenter (ie,
##   by maintaining atomic clock corrections)## (b) estimating the
##   observatory position vector## and and (c) estimating earth's
##   velocity vector (using JPLEPHINTERP).
##
##   Users may note a circularity to the above equation, since
##   vearth(T) is expressed in terms of the SSB coordinate time.  This
##   appears to be a chicken and egg problem since in order to get the
##   earth's velocity, the ephemeris time is needed to begin with.
##   However, to the precision of the above equation, < 25 ns, it is
##   acceptable to replace vearth(T) with vearth(TT).
##
##   The method of computation of TDB2TDT in this function is based on
##   the analytical formulation by Fairhead, Bretagnon & Lestrade, 1988
##   (so-called FBL model) and Fairhead & Bretagnon 1990, in terms of
##   sinusoids of various amplitudes.  TDB2TDT has a dominant periodic
##   component of period 1 year and amplitude 1.7 ms.  The set of 791
##   coefficients used here were drawn from the Princeton pulsar timing
##   program TEMPO version 11.005 (Taylor & Weisberg 1989).
##
##   Because the TDB2TDT quantity is rather expensive to compute but
##   slowly varying, users may wish to also retrieve the time
##   derivative using the DERIV keyword, if they have many times to
##   convert over a short baseline.
##
## Verification
##
##   This implementation has been compared against a set of FBL test
##   data found in the 1996 IERS Conventions, Chapter 11, provided by
##   T. Fukushima.  It has been verified that this routine reproduces
##   the Fukushima numbers to the accuracy of the table, within
##   10^{-14} seconds.
##
##   Fukushima (1995) has found that the 791-term Fairhead & Bretagnon
##   analytical approximation use here has a maximum error of 23
##   nanoseconds in the time range 1980-2000, compared to a numerical
##   integration.  In comparison the truncated 127-term approximation
##   has an error of ~130 nanoseconds.
##
##
## PARAMETERS:
##
##   JD - Geocentric time TT, scalar or vector, expressed in Julian
##        days.  The actual time used is (JD + TBASE).  For maximum
##        precision, TBASE should be used to express a fixed epoch in
##        whole day numbers, and JD should express fractional offset
##        days from that epoch.
##
##
## KEYWORD PARAMETERS:
##
##   TBASE - scalar Julian day of a fixed epoch, which provides the
##           origin for times passed in JD.
##          Default: 0
##
##   DERIV - upon return, contains the derivative of TDB2TDT in units
##           of seconds per day.  As many derivatives are returned as
##           values passed in JD.
##
##
## RETURNS:
##   The correction offset(s) in units of seconds, to be applied as
##   noted above.
##
##
## EXAMPLE:
##
##   Find the correction at ephemeris time 2451544.5 (JD):
##     IDL> print, tdb2tdt(2451544.5d)
##       -0.00011376314
##   or 0.11 ms.
##
##
## REFERENCES:
##
##   Princeton TEMPO Program
##      http://tempo.sourceforge.net/tempo_idx.html
##
##   FBL Test Data Set
##      ftp://maia.usno.navy.mil/conventions/chapter11/fbl.results
##
##   Fairhead, L. & Bretagnon, P. 1990, A&A, 229, 240
##     (basis of this routine)
##
##   Fairhead, L. Bretagnon, P. & Lestrade, J.-F. 1988, in *The Earth's
##     Rotation and Reference Frames for Geodesy and Geodynamics*,
##     ed. A. K. Babcock and G. A. Wilkins, (Dordrecht: Kluwer), p. 419
##     (original "FBL" paper)
##
##   Fukushima, T. 1995, A&A, 294, 895  (error analysis)
##
##   Irwin, A. W. & Fukushima, T. 1999, A&A, 348, 642  (error analysis)
##
##   Standish, E. M. 1998, A&A, 336, 381 (description of time scales)
##
##   Taylor, J. H. & Weisberg, J. M. 1989, ApJ, 345, 434 (pulsar timing)
##
##
## SEE ALSO
##   JPLEPHREAD, JPLEPHINTERP, JPLEPHTEST
##
## MODIFICATION HISTORY:
##   Original logic from Fairhead & Bretagnon, 1990
##   Drawn from TEMPO v. 11.005, copied 20 Jun 2001
##   Documented and vectorized, 30 Jun 2001
##
##
##  $Id: tdb2tdt.pro,v 1.4 2001/07/01 07:37:40 craigm Exp $
##
##-
## Copyright (C) 2001, Craig Markwardt
## This software is provided as is without any warranty whatsoever.
## Permission to use, copy and distribute unmodified copies for
## non-commercial purposes, and to modify and use for personal or
## internal use, is granted.  All other rights are reserved.
####################################
fbldata = rbind(c(
1656.674564,   6283.075849991, 6.240054195),c(
  22.417471,   5753.384884897, 4.296977442),c(
  13.839792,  12566.151699983, 6.196904410),c(
   4.770086,    529.690965095, 0.444401603),c(
   4.676740,   6069.776754553, 4.021195093),c(
   2.256707,    213.299095438, 5.543113262),c(
   1.694205,     -3.523118349, 5.025132748),c(
   1.554905,  77713.771467920, 5.198467090),c(
   1.276839,   7860.419392439, 5.988822341),c(
   1.193379,   5223.693919802, 3.649823730),c(
   1.115322,   3930.209696220, 1.422745069),c(
   0.794185,  11506.769769794, 2.322313077),c(
   0.447061,     26.298319800, 3.615796498),c(
   0.435206,   -398.149003408, 4.349338347),c(
   0.600309,   1577.343542448, 2.678271909),c(
   0.496817,   6208.294251424, 5.696701824),c(
   0.486306,   5884.926846583, 0.520007179),c(
   0.432392,     74.781598567, 2.435898309),c(
   0.468597,   6244.942814354, 5.866398759),c(
   0.375510,   5507.553238667, 4.103476804),c(
   0.243085,   -775.522611324, 3.651837925),c(
   0.173435,  18849.227549974, 6.153743485),c(
   0.230685,   5856.477659115, 4.773852582),c(
   0.203747,  12036.460734888, 4.333987818),c(
   0.143935,   -796.298006816, 5.957517795  ))
fbldata = rbind( fbldata,c(
   0.159080,  10977.078804699, 1.890075226),c(
   0.119979,     38.133035638, 4.551585768),c(
   0.118971,   5486.777843175, 1.914547226),c(
   0.116120,   1059.381930189, 0.873504123),c(
   0.137927,  11790.629088659, 1.135934669),c(
   0.098358,   2544.314419883, 0.092793886),c(
   0.101868,  -5573.142801634, 5.984503847),c(
   0.080164,    206.185548437, 2.095377709),c(
   0.079645,   4694.002954708, 2.949233637),c(
   0.062617,     20.775395492, 2.654394814),c(
   0.075019,   2942.463423292, 4.980931759),c(
   0.064397,   5746.271337896, 1.280308748),c(
   0.063814,   5760.498431898, 4.167901731),c(
   0.048042,   2146.165416475, 1.495846011),c(
   0.048373,    155.420399434, 2.251573730),c(
   0.058844,    426.598190876, 4.839650148),c(
   0.046551,     -0.980321068, 0.921573539),c(
   0.054139,  17260.154654690, 3.411091093),c(
   0.042411,   6275.962302991, 2.869567043),c(
   0.040184,     -7.113547001, 3.565975565),c(
   0.036564,   5088.628839767, 3.324679049),c(
   0.040759,  12352.852604545, 3.981496998),c(
   0.036507,    801.820931124, 6.248866009),c(
   0.036955,   3154.687084896, 5.071801441),c(
   0.042732,    632.783739313, 5.720622217  ))
fbldata = rbind( fbldata,c(
   0.042560, 161000.685737473, 1.270837679),c(
   0.040480,  15720.838784878, 2.546610123),c(
   0.028244,  -6286.598968340, 5.069663519),c(
   0.033477,   6062.663207553, 4.144987272),c(
   0.034867,    522.577418094, 5.210064075),c(
   0.032438,   6076.890301554, 0.749317412),c(
   0.030215,   7084.896781115, 3.389610345),c(
   0.029247, -71430.695617928, 4.183178762),c(
   0.033529,   9437.762934887, 2.404714239),c(
   0.032423,   8827.390269875, 5.541473556),c(
   0.027567,   6279.552731642, 5.040846034),c(
   0.029862,  12139.553509107, 1.770181024),c(
   0.022509,  10447.387839604, 1.460726241),c(
   0.020937,   8429.241266467, 0.652303414),c(
   0.020322,    419.484643875, 3.735430632),c(
   0.024816,  -1194.447010225, 1.087136918),c(
   0.025196,   1748.016413067, 2.901883301),c(
   0.021691,  14143.495242431, 5.952658009),c(
   0.017673,   6812.766815086, 3.186129845),c(
   0.022567,   6133.512652857, 3.307984806),c(
   0.016155,  10213.285546211, 1.331103168),c(
   0.014751,   1349.867409659, 4.308933301),c(
   0.015949,   -220.412642439, 4.005298270),c(
   0.015974,  -2352.866153772, 6.145309371),c(
   0.014223,  17789.845619785, 2.104551349  ))
fbldata = rbind( fbldata,c(
   0.017806,     73.297125859, 3.475975097),c(
   0.013671,   -536.804512095, 5.971672571),c(
   0.011942,   8031.092263058, 2.053414715),c(
   0.014318,  16730.463689596, 3.016058075),c(
   0.012462,    103.092774219, 1.737438797),c(
   0.010962,      3.590428652, 2.196567739),c(
   0.015078,  19651.048481098, 3.969480770),c(
   0.010396,    951.718406251, 5.717799605),c(
   0.011707,  -4705.732307544, 2.654125618),c(
   0.010453,   5863.591206116, 1.913704550),c(
   0.012420,   4690.479836359, 4.734090399),c(
   0.011847,   5643.178563677, 5.489005403),c(
   0.008610,   3340.612426700, 3.661698944),c(
   0.011622,   5120.601145584, 4.863931876),c(
   0.010825,    553.569402842, 0.842715011),c(
   0.008666,   -135.065080035, 3.293406547),c(
   0.009963,    149.563197135, 4.870690598),c(
   0.009858,   6309.374169791, 1.061816410),c(
   0.007959,    316.391869657, 2.465042647),c(
   0.010099,    283.859318865, 1.942176992),c(
   0.007147,   -242.728603974, 3.661486981),c(
   0.007505,   5230.807466803, 4.920937029),c(
   0.008323,  11769.853693166, 1.229392026),c(
   0.007490,  -6256.777530192, 3.658444681),c(
   0.009370, 149854.400134205, 0.673880395  ))
fbldata = rbind( fbldata,c(
   0.007117,     38.027672636, 5.294249518),c(
   0.007857,  12168.002696575, 0.525733528),c(
   0.007019,   6206.809778716, 0.837688810),c(
   0.006056,    955.599741609, 4.194535082),c(
   0.008107,  13367.972631107, 3.793235253),c(
   0.006731,   5650.292110678, 5.639906583),c(
   0.007332,     36.648562930, 0.114858677),c(
   0.006366,   4164.311989613, 2.262081818),c(
   0.006858,   5216.580372801, 0.642063318),c(
   0.006919,   6681.224853400, 6.018501522),c(
   0.006826,   7632.943259650, 3.458654112),c(
   0.005308,  -1592.596013633, 2.500382359),c(
   0.005096,  11371.704689758, 2.547107806),c(
   0.004841,   5333.900241022, 0.437078094),c(
   0.005582,   5966.683980335, 2.246174308),c(
   0.006304,  11926.254413669, 2.512929171),c(
   0.006603,  23581.258177318, 5.393136889),c(
   0.005123,     -1.484472708, 2.999641028),c(
   0.004648,   1589.072895284, 1.275847090),c(
   0.005119,   6438.496249426, 1.486539246),c(
   0.004521,   4292.330832950, 6.140635794),c(
   0.005680,  23013.539539587, 4.557814849),c(
   0.005488,     -3.455808046, 0.090675389),c(
   0.004193,   7234.794256242, 4.869091389),c(
   0.003742,   7238.675591600, 4.691976180  ))
fbldata = rbind( fbldata,c(
   0.004148,   -110.206321219, 3.016173439),c(
   0.004553,  11499.656222793, 5.554998314),c(
   0.004892,   5436.993015240, 1.475415597),c(
   0.004044,   4732.030627343, 1.398784824),c(
   0.004164,  12491.370101415, 5.650931916),c(
   0.004349,  11513.883316794, 2.181745369),c(
   0.003919,  12528.018664345, 5.823319737),c(
   0.003129,   6836.645252834, 0.003844094),c(
   0.004080,  -7058.598461315, 3.690360123),c(
   0.003270,     76.266071276, 1.517189902),c(
   0.002954,   6283.143160294, 4.447203799),c(
   0.002872,     28.449187468, 1.158692983),c(
   0.002881,    735.876513532, 0.349250250),c(
   0.003279,   5849.364112115, 4.893384368),c(
   0.003625,   6209.778724132, 1.473760578),c(
   0.003074,    949.175608970, 5.185878737),c(
   0.002775,   9917.696874510, 1.030026325),c(
   0.002646,  10973.555686350, 3.918259169),c(
   0.002575,  25132.303399966, 6.109659023),c(
   0.003500,    263.083923373, 1.892100742),c(
   0.002740,  18319.536584880, 4.320519510),c(
   0.002464,    202.253395174, 4.698203059),c(
   0.002409,      2.542797281, 5.325009315),c(
   0.003354, -90955.551694697, 1.942656623),c(
   0.002296,   6496.374945429, 5.061810696  ))
fbldata = rbind( fbldata,c(
   0.003002,   6172.869528772, 2.797822767),c(
   0.003202,  27511.467873537, 0.531673101),c(
   0.002954,  -6283.008539689, 4.533471191),c(
   0.002353,    639.897286314, 3.734548088),c(
   0.002401,  16200.772724501, 2.605547070),c(
   0.003053, 233141.314403759, 3.029030662),c(
   0.003024,  83286.914269554, 2.355556099),c(
   0.002863,  17298.182327326, 5.240963796),c(
   0.002103,  -7079.373856808, 5.756641637),c(
   0.002303,  83996.847317911, 2.013686814),c(
   0.002303,  18073.704938650, 1.089100410),c(
   0.002381,     63.735898303, 0.759188178),c(
   0.002493,   6386.168624210, 0.645026535),c(
   0.002366,      3.932153263, 6.215885448),c(
   0.002169,  11015.106477335, 4.845297676),c(
   0.002397,   6243.458341645, 3.809290043),c(
   0.002183,   1162.474704408, 6.179611691),c(
   0.002353,   6246.427287062, 4.781719760),c(
   0.002199,   -245.831646229, 5.956152284),c(
   0.001729,   3894.181829542, 1.264976635),c(
   0.001896,  -3128.388765096, 4.914231596),c(
   0.002085,     35.164090221, 1.405158503),c(
   0.002024,  14712.317116458, 2.752035928),c(
   0.001737,   6290.189396992, 5.280820144),c(
   0.002229,    491.557929457, 1.571007057 ))
fbldata = rbind( fbldata,c(
   0.001602,  14314.168113050, 4.203664806),c(
   0.002186,    454.909366527, 1.402101526),c(
   0.001897,  22483.848574493, 4.167932508),c(
   0.001825,  -3738.761430108, 0.545828785),c(
   0.001894,   1052.268383188, 5.817167450),c(
   0.001421,     20.355319399, 2.419886601),c(
   0.001408,  10984.192351700, 2.732084787),c(
   0.001847,  10873.986030480, 2.903477885),c(
   0.001391,  -8635.942003763, 0.593891500),c(
   0.001388,     -7.046236698, 1.166145902),c(
   0.001810, -88860.057071188, 0.487355242),c(
   0.001288,  -1990.745017041, 3.913022880),c(
   0.001297,  23543.230504682, 3.063805171),c(
   0.001335,   -266.607041722, 3.995764039),c(
   0.001376,  10969.965257698, 5.152914309),c(
   0.001745, 244287.600007027, 3.626395673),c(
   0.001649,  31441.677569757, 1.952049260),c(
   0.001416,   9225.539273283, 4.996408389),c(
   0.001238,   4804.209275927, 5.503379738),c(
   0.001472,   4590.910180489, 4.164913291),c(
   0.001169,   6040.347246017, 5.841719038),c(
   0.001039,   5540.085789459, 2.769753519),c(
   0.001004,   -170.672870619, 0.755008103),c(
   0.001284,  10575.406682942, 5.306538209),c(
   0.001278,     71.812653151, 4.713486491 ))
fbldata = rbind( fbldata,c(
   0.001321,  18209.330263660, 2.624866359),c(
   0.001297,  21228.392023546, 0.382603541),c(
   0.000954,   6282.095528923, 0.882213514),c(
   0.001145,   6058.731054289, 1.169483931),c(
   0.000979,   5547.199336460, 5.448375984),c(
   0.000987,  -6262.300454499, 2.656486959),c(
   0.001070,-154717.609887482, 1.827624012),c(
   0.000991,   4701.116501708, 4.387001801),c(
   0.001155,    -14.227094002, 3.042700750),c(
   0.001176,    277.034993741, 3.335519004),c(
   0.000890,  13916.019109642, 5.601498297),c(
   0.000884,  -1551.045222648, 1.088831705),c(
   0.000876,   5017.508371365, 3.969902609),c(
   0.000806,  15110.466119866, 5.142876744),c(
   0.000773,  -4136.910433516, 0.022067765),c(
   0.001077,    175.166059800, 1.844913056),c(
   0.000954,  -6284.056171060, 0.968480906),c(
   0.000737,   5326.786694021, 4.923831588),c(
   0.000845,   -433.711737877, 4.749245231),c(
   0.000819,   8662.240323563, 5.991247817),c(
   0.000852,    199.072001436, 2.189604979),c(
   0.000723,  17256.631536341, 6.068719637),c(
   0.000940,   6037.244203762, 6.197428148),c(
   0.000885,  11712.955318231, 3.280414875),c(
   0.000706,  12559.038152982, 2.824848947 ))
fbldata = rbind( fbldata,c(
   0.000732,   2379.164473572, 2.501813417),c(
   0.000764,  -6127.655450557, 2.236346329),c(
   0.000908,    131.541961686, 2.521257490),c(
   0.000907,  35371.887265976, 3.370195967),c(
   0.000673,   1066.495477190, 3.876512374),c(
   0.000814,  17654.780539750, 4.627122566),c(
   0.000630,     36.027866677, 0.156368499),c(
   0.000798,    515.463871093, 5.151962502),c(
   0.000798,    148.078724426, 5.909225055),c(
   0.000806,    309.278322656, 6.054064447),c(
   0.000607,    -39.617508346, 2.839021623),c(
   0.000601,    412.371096874, 3.984225404),c(
   0.000646,  11403.676995575, 3.852959484),c(
   0.000704,  13521.751441591, 2.300991267),c(
   0.000603, -65147.619767937, 4.140083146),c(
   0.000609,  10177.257679534, 0.437122327),c(
   0.000631,   5767.611978898, 4.026532329),c(
   0.000576,  11087.285125918, 4.760293101),c(
   0.000674,  14945.316173554, 6.270510511),c(
   0.000726,   5429.879468239, 6.039606892),c(
   0.000710,  28766.924424484, 5.672617711),c(
   0.000647,  11856.218651625, 3.397132627),c(
   0.000678,  -5481.254918868, 6.249666675),c(
   0.000618,  22003.914634870, 2.466427018),c(
   0.000738,   6134.997125565, 2.242668890 ))
fbldata = rbind( fbldata,c(
   0.000660,    625.670192312, 5.864091907),c(
   0.000694,   3496.032826134, 2.668309141),c(
   0.000531,   6489.261398429, 1.681888780),c(
   0.000611,-143571.324284214, 2.424978312),c(
   0.000575,  12043.574281889, 4.216492400),c(
   0.000553,  12416.588502848, 4.772158039),c(
   0.000689,   4686.889407707, 6.224271088),c(
   0.000495,   7342.457780181, 3.817285811),c(
   0.000567,   3634.621024518, 1.649264690),c(
   0.000515,  18635.928454536, 3.945345892),c(
   0.000486,   -323.505416657, 4.061673868),c(
   0.000662,  25158.601719765, 1.794058369),c(
   0.000509,    846.082834751, 3.053874588),c(
   0.000472, -12569.674818332, 5.112133338),c(
   0.000461,   6179.983075773, 0.513669325),c(
   0.000641,  83467.156352816, 3.210727723),c(
   0.000520,  10344.295065386, 2.445597761),c(
   0.000493,  18422.629359098, 1.676939306),c(
   0.000478,   1265.567478626, 5.487314569),c(
   0.000472,    -18.159247265, 1.999707589),c(
   0.000559,  11190.377900137, 5.783236356),c(
   0.000494,   9623.688276691, 3.022645053),c(
   0.000463,   5739.157790895, 1.411223013),c(
   0.000432,  16858.482532933, 1.179256434),c(
   0.000574,  72140.628666286, 1.758191830 ))
fbldata = rbind( fbldata,c(
   0.000484,  17267.268201691, 3.290589143),c(
   0.000550,   4907.302050146, 0.864024298),c(
   0.000399,     14.977853527, 2.094441910),c(
   0.000491,    224.344795702, 0.878372791),c(
   0.000432,  20426.571092422, 6.003829241),c(
   0.000481,   5749.452731634, 4.309591964),c(
   0.000480,   5757.317038160, 1.142348571),c(
   0.000485,   6702.560493867, 0.210580917),c(
   0.000426,   6055.549660552, 4.274476529),c(
   0.000480,   5959.570433334, 5.031351030),c(
   0.000466,  12562.628581634, 4.959581597),c(
   0.000520,  39302.096962196, 4.788002889),c(
   0.000458,  12132.439962106, 1.880103788),c(
   0.000470,  12029.347187887, 1.405611197),c(
   0.000416,  -7477.522860216, 1.082356330),c(
   0.000449,  11609.862544012, 4.179989585),c(
   0.000465,  17253.041107690, 0.353496295),c(
   0.000362,  -4535.059436924, 1.583849576),c(
   0.000383,  21954.157609398, 3.747376371),c(
   0.000389,     17.252277143, 1.395753179),c(
   0.000331,  18052.929543158, 0.566790582),c(
   0.000430,  13517.870106233, 0.685827538),c(
   0.000368,  -5756.908003246, 0.731374317),c(
   0.000330,  10557.594160824, 3.710043680),c(
   0.000332,  20199.094959633, 1.652901407 ))
fbldata = rbind( fbldata,c(
   0.000384,  11933.367960670, 5.827781531),c(
   0.000387,  10454.501386605, 2.541182564),c(
   0.000325,  15671.081759407, 2.178850542),c(
   0.000318,    138.517496871, 2.253253037),c(
   0.000305,   9388.005909415, 0.578340206),c(
   0.000352,   5749.861766548, 3.000297967),c(
   0.000311,   6915.859589305, 1.693574249),c(
   0.000297,  24072.921469776, 1.997249392),c(
   0.000363,   -640.877607382, 5.071820966),c(
   0.000323,  12592.450019783, 1.072262823),c(
   0.000341,  12146.667056108, 4.700657997),c(
   0.000290,   9779.108676125, 1.812320441),c(
   0.000342,   6132.028180148, 4.322238614),c(
   0.000329,   6268.848755990, 3.033827743),c(
   0.000374,  17996.031168222, 3.388716544),c(
   0.000285,   -533.214083444, 4.687313233),c(
   0.000338,   6065.844601290, 0.877776108),c(
   0.000276,     24.298513841, 0.770299429),c(
   0.000336,  -2388.894020449, 5.353796034),c(
   0.000290,   3097.883822726, 4.075291557),c(
   0.000318,    709.933048357, 5.941207518),c(
   0.000271,  13095.842665077, 3.208912203),c(
   0.000331,   6073.708907816, 4.007881169),c(
   0.000292,    742.990060533, 2.714333592),c(
   0.000362,  29088.811415985, 3.215977013 ))
fbldata = rbind( fbldata,c(
   0.000280,  12359.966151546, 0.710872502),c(
   0.000267,  10440.274292604, 4.730108488),c(
   0.000262,    838.969287750, 1.327720272),c(
   0.000250,  16496.361396202, 0.898769761),c(
   0.000325,  20597.243963041, 0.180044365),c(
   0.000268,   6148.010769956, 5.152666276),c(
   0.000284,   5636.065016677, 5.655385808),c(
   0.000301,   6080.822454817, 2.135396205),c(
   0.000294,   -377.373607916, 3.708784168),c(
   0.000236,   2118.763860378, 1.733578756),c(
   0.000234,   5867.523359379, 5.575209112),c(
   0.000268,-226858.238553767, 0.069432392),c(
   0.000265, 167283.761587465, 4.369302826),c(
   0.000280,  28237.233459389, 5.304829118),c(
   0.000292,  12345.739057544, 4.096094132),c(
   0.000223,  19800.945956225, 3.069327406),c(
   0.000301,  43232.306658416, 6.205311188),c(
   0.000264,  18875.525869774, 1.417263408),c(
   0.000304,  -1823.175188677, 3.409035232),c(
   0.000301,    109.945688789, 0.510922054),c(
   0.000260,    813.550283960, 2.389438934),c(
   0.000299, 316428.228673312, 5.384595078),c(
   0.000211,   5756.566278634, 3.789392838),c(
   0.000209,   5750.203491159, 1.661943545),c(
   0.000240,  12489.885628707, 5.684549045 ))
fbldata = rbind( fbldata,c(
   0.000216,   6303.851245484, 3.862942261),c(
   0.000203,   1581.959348283, 5.549853589),c(
   0.000200,   5642.198242609, 1.016115785),c(
   0.000197,    -70.849445304, 4.690702525),c(
   0.000227,   6287.008003254, 2.911891613),c(
   0.000197,    533.623118358, 1.048982898),c(
   0.000205,  -6279.485421340, 1.829362730),c(
   0.000209, -10988.808157535, 2.636140084),c(
   0.000208,   -227.526189440, 4.127883842),c(
   0.000191,    415.552490612, 4.401165650),c(
   0.000190,  29296.615389579, 4.175658539),c(
   0.000264,  66567.485864652, 4.601102551),c(
   0.000256,  -3646.350377354, 0.506364778),c(
   0.000188,  13119.721102825, 2.032195842),c(
   0.000185,   -209.366942175, 4.694756586),c(
   0.000198,  25934.124331089, 3.832703118),c(
   0.000195,   4061.219215394, 3.308463427),c(
   0.000234,   5113.487598583, 1.716090661),c(
   0.000188,   1478.866574064, 5.686865780),c(
   0.000222,  11823.161639450, 1.942386641),c(
   0.000181,  10770.893256262, 1.999482059),c(
   0.000171,   6546.159773364, 1.182807992),c(
   0.000206,     70.328180442, 5.934076062),c(
   0.000169,  20995.392966449, 2.169080622),c(
   0.000191,  10660.686935042, 5.405515999 ))
fbldata = rbind( fbldata,c(
   0.000228,  33019.021112205, 4.656985514),c(
   0.000184,  -4933.208440333, 3.327476868),c(
   0.000220,   -135.625325010, 1.765430262),c(
   0.000166,  23141.558382925, 3.454132746),c(
   0.000191,   6144.558353121, 5.020393445),c(
   0.000180,   6084.003848555, 0.602182191),c(
   0.000163,  17782.732072784, 4.960593133),c(
   0.000225,  16460.333529525, 2.596451817),c(
   0.000222,   5905.702242076, 3.731990323),c(
   0.000204,    227.476132789, 5.636192701),c(
   0.000159,  16737.577236597, 3.600691544),c(
   0.000200,   6805.653268085, 0.868220961),c(
   0.000187,  11919.140866668, 2.629456641),c(
   0.000161,    127.471796607, 2.862574720),c(
   0.000205,   6286.666278643, 1.742882331),c(
   0.000189,    153.778810485, 4.812372643),c(
   0.000168,  16723.350142595, 0.027860588),c(
   0.000149,  11720.068865232, 0.659721876),c(
   0.000189,   5237.921013804, 5.245313000),c(
   0.000143,   6709.674040867, 4.317625647),c(
   0.000146,   4487.817406270, 4.815297007),c(
   0.000144,   -664.756045130, 5.381366880),c(
   0.000175,   5127.714692584, 4.728443327),c(
   0.000162,   6254.626662524, 1.435132069),c(
   0.000187,  47162.516354635, 1.354371923 ))
fbldata = rbind( fbldata,c(
   0.000146,  11080.171578918, 3.369695406),c(
   0.000180,   -348.924420448, 2.490902145),c(
   0.000148,    151.047669843, 3.799109588),c(
   0.000157,   6197.248551160, 1.284375887),c(
   0.000167,    146.594251718, 0.759969109),c(
   0.000133,  -5331.357443741, 5.409701889),c(
   0.000154,     95.979227218, 3.366890614),c(
   0.000148,  -6418.140930027, 3.384104996),c(
   0.000128,  -6525.804453965, 3.803419985),c(
   0.000130,  11293.470674356, 0.939039445),c(
   0.000152,  -5729.506447149, 0.734117523),c(
   0.000138,    210.117701700, 2.564216078),c(
   0.000123,   6066.595360816, 4.517099537),c(
   0.000140,  18451.078546566, 0.642049130),c(
   0.000126,  11300.584221356, 3.485280663),c(
   0.000119,  10027.903195729, 3.217431161),c(
   0.000151,   4274.518310832, 4.404359108),c(
   0.000117,   6072.958148291, 0.366324650),c(
   0.000165,  -7668.637425143, 4.298212528),c(
   0.000117,  -6245.048177356, 5.379518958),c(
   0.000130,  -5888.449964932, 4.527681115),c(
   0.000121,   -543.918059096, 6.109429504),c(
   0.000162,   9683.594581116, 5.720092446),c(
   0.000141,   6219.339951688, 0.679068671),c(
   0.000118,  22743.409379516, 4.881123092 ))
fbldata = rbind( fbldata,c(
   0.000129,   1692.165669502, 0.351407289),c(
   0.000126,   5657.405657679, 5.146592349),c(
   0.000114,    728.762966531, 0.520791814),c(
   0.000120,     52.596639600, 0.948516300),c(
   0.000115,     65.220371012, 3.504914846),c(
   0.000126,   5881.403728234, 5.577502482),c(
   0.000158, 163096.180360983, 2.957128968),c(
   0.000134,  12341.806904281, 2.598576764),c(
   0.000151,  16627.370915377, 3.985702050),c(
   0.000109,   1368.660252845, 0.014730471),c(
   0.000131,   6211.263196841, 0.085077024),c(
   0.000146,   5792.741760812, 0.708426604),c(
   0.000146,    -77.750543984, 3.121576600),c(
   0.000107,   5341.013788022, 0.288231904),c(
   0.000138,   6281.591377283, 2.797450317),c(
   0.000113,  -6277.552925684, 2.788904128),c(
   0.000115,   -525.758811831, 5.895222200),c(
   0.000138,   6016.468808270, 6.096188999),c(
   0.000139,  23539.707386333, 2.028195445),c(
   0.000146,  -4176.041342449, 4.660008502),c(
   0.000107,  16062.184526117, 4.066520001),c(
   0.000142,  83783.548222473, 2.936315115),c(
   0.000128,   9380.959672717, 3.223844306),c(
   0.000135,   6205.325306007, 1.638054048),c(
   0.000101,   2699.734819318, 5.481603249 ))
fbldata = rbind( fbldata,c(
   0.000104,   -568.821874027, 2.205734493),c(
   0.000103,   6321.103522627, 2.440421099),c(
   0.000119,   6321.208885629, 2.547496264),c(
   0.000138,   1975.492545856, 2.314608466),c(
   0.000121,    137.033024162, 4.539108237),c(
   0.000123,  19402.796952817, 4.538074405),c(
   0.000119,  22805.735565994, 2.869040566),c(
   0.000133,  64471.991241142, 6.056405489),c(
   0.000129,    -85.827298831, 2.540635083),c(
   0.000131,  13613.804277336, 4.005732868),c(
   0.000104,   9814.604100291, 1.959967212),c(
   0.000112,  16097.679950283, 3.589026260),c(
   0.000123,   2107.034507542, 1.728627253),c(
   0.000121,  36949.230808424, 6.072332087),c(
   0.000108, -12539.853380183, 3.716133846),c(
   0.000113,  -7875.671863624, 2.725771122),c(
   0.000109,   4171.425536614, 4.033338079),c(
   0.000101,   6247.911759770, 3.441347021),c(
   0.000113,   7330.728427345, 0.656372122),c(
   0.000113,  51092.726050855, 2.791483066),c(
   0.000106,   5621.842923210, 1.815323326),c(
   0.000101,    111.430161497, 5.711033677),c(
   0.000103,    909.818733055, 2.812745443),c(
   0.000101,   1790.642637886, 1.965746028 ))
fbldata = rbind( fbldata,c(  #### From enof TDB1NS.F
   0.00065,    6069.776754,    4.021194),c(
   0.00033,     213.299095,    5.543132),c(
  -0.00196,    6208.294251,    5.696701),c(
  -0.00173,      74.781599,    2.435900 ))

i1terms = nrow(fbldata)
## T**1
fbldata = rbind( fbldata,c(
 102.156724,   6283.075849991, 4.249032005),c(
   1.706807,  12566.151699983, 4.205904248),c(
   0.269668,    213.299095438, 3.400290479),c(
   0.265919,    529.690965095, 5.836047367),c(
   0.210568,     -3.523118349, 6.262738348),c(
   0.077996,   5223.693919802, 4.670344204),c(
   0.054764,   1577.343542448, 4.534800170),c(
   0.059146,     26.298319800, 1.083044735),c(
   0.034420,   -398.149003408, 5.980077351),c(
   0.032088,  18849.227549974, 4.162913471),c(
   0.033595,   5507.553238667, 5.980162321),c(
   0.029198,   5856.477659115, 0.623811863),c(
   0.027764,    155.420399434, 3.745318113),c(
   0.025190,   5746.271337896, 2.980330535),c(
   0.022997,   -796.298006816, 1.174411803),c(
   0.024976,   5760.498431898, 2.467913690),c(
   0.021774,    206.185548437, 3.854787540),c(
   0.017925,   -775.522611324, 1.092065955),c(
   0.013794,    426.598190876, 2.699831988),c(
   0.013276,   6062.663207553, 5.845801920),c(
   0.011774,  12036.460734888, 2.292832062),c(
   0.012869,   6076.890301554, 5.333425680),c(
   0.012152,   1059.381930189, 6.222874454),c(
   0.011081,     -7.113547001, 5.154724984),c(
   0.010143,   4694.002954708, 4.044013795 ))
fbldata = rbind( fbldata,c(
   0.009357,   5486.777843175, 3.416081409),c(
   0.010084,    522.577418094, 0.749320262),c(
   0.008587,  10977.078804699, 2.777152598),c(
   0.008628,   6275.962302991, 4.562060226),c(
   0.008158,   -220.412642439, 5.806891533),c(
   0.007746,   2544.314419883, 1.603197066),c(
   0.007670,   2146.165416475, 3.000200440),c(
   0.007098,     74.781598567, 0.443725817),c(
   0.006180,   -536.804512095, 1.302642751),c(
   0.005818,   5088.628839767, 4.827723531),c(
   0.004945,  -6286.598968340, 0.268305170),c(
   0.004774,   1349.867409659, 5.808636673),c(
   0.004687,   -242.728603974, 5.154890570),c(
   0.006089,   1748.016413067, 4.403765209),c(
   0.005975,  -1194.447010225, 2.583472591),c(
   0.004229,    951.718406251, 0.931172179),c(
   0.005264,    553.569402842, 2.336107252),c(
   0.003049,   5643.178563677, 1.362634430),c(
   0.002974,   6812.766815086, 1.583012668),c(
   0.003403,  -2352.866153772, 2.552189886),c(
   0.003030,    419.484643875, 5.286473844),c(
   0.003210,     -7.046236698, 1.863796539),c(
   0.003058,   9437.762934887, 4.226420633),c(
   0.002589,  12352.852604545, 1.991935820),c(
   0.002927,   5216.580372801, 2.319951253 ))
fbldata = rbind( fbldata,c(
   0.002425,   5230.807466803, 3.084752833),c(
   0.002656,   3154.687084896, 2.487447866),c(
   0.002445,  10447.387839604, 2.347139160),c(
   0.002990,   4690.479836359, 6.235872050),c(
   0.002890,   5863.591206116, 0.095197563),c(
   0.002498,   6438.496249426, 2.994779800),c(
   0.001889,   8031.092263058, 3.569003717),c(
   0.002567,    801.820931124, 3.425611498),c(
   0.001803, -71430.695617928, 2.192295512),c(
   0.001782,      3.932153263, 5.180433689),c(
   0.001694,  -4705.732307544, 4.641779174),c(
   0.001704,  -1592.596013633, 3.997097652),c(
   0.001735,   5849.364112115, 0.417558428),c(
   0.001643,   8429.241266467, 2.180619584),c(
   0.001680,     38.133035638, 4.164529426),c(
   0.002045,   7084.896781115, 0.526323854),c(
   0.001458,   4292.330832950, 1.356098141),c(
   0.001437,     20.355319399, 3.895439360),c(
   0.001738,   6279.552731642, 0.087484036),c(
   0.001367,  14143.495242431, 3.987576591),c(
   0.001344,   7234.794256242, 0.090454338),c(
   0.001438,  11499.656222793, 0.974387904),c(
   0.001257,   6836.645252834, 1.509069366),c(
   0.001358,  11513.883316794, 0.495572260),c(
   0.001628,   7632.943259650, 4.968445721 ))
fbldata = rbind( fbldata,c(
   0.001169,    103.092774219, 2.838496795),c(
   0.001162,   4164.311989613, 3.408387778),c(
   0.001092,   6069.776754553, 3.617942651),c(
   0.001008,  17789.845619785, 0.286350174),c(
   0.001008,    639.897286314, 1.610762073),c(
   0.000918,  10213.285546211, 5.532798067),c(
   0.001011,  -6256.777530192, 0.661826484),c(
   0.000753,  16730.463689596, 3.905030235),c(
   0.000737,  11926.254413669, 4.641956361),c(
   0.000694,   3340.612426700, 2.111120332),c(
   0.000701,   3894.181829542, 2.760823491),c(
   0.000689,   -135.065080035, 4.768800780),c(
   0.000700,  13367.972631107, 5.760439898),c(
   0.000664,   6040.347246017, 1.051215840),c(
   0.000654,   5650.292110678, 4.911332503),c(
   0.000788,   6681.224853400, 4.699648011),c(
   0.000628,   5333.900241022, 5.024608847),c(
   0.000755,   -110.206321219, 4.370971253),c(
   0.000628,   6290.189396992, 3.660478857),c(
   0.000635,  25132.303399966, 4.121051532),c(
   0.000534,   5966.683980335, 1.173284524),c(
   0.000543,   -433.711737877, 0.345585464),c(
   0.000517,  -1990.745017041, 5.414571768),c(
   0.000504,   5767.611978898, 2.328281115),c(
   0.000485,   5753.384884897, 1.685874771 ))
fbldata = rbind( fbldata,c(
   0.000463,   7860.419392439, 5.297703006),c(
   0.000604,    515.463871093, 0.591998446),c(
   0.000443,  12168.002696575, 4.830881244),c(
   0.000570,    199.072001436, 3.899190272),c(
   0.000465,  10969.965257698, 0.476681802),c(
   0.000424,  -7079.373856808, 1.112242763),c(
   0.000427,    735.876513532, 1.994214480),c(
   0.000478,  -6127.655450557, 3.778025483),c(
   0.000414,  10973.555686350, 5.441088327),c(
   0.000512,   1589.072895284, 0.107123853),c(
   0.000378,  10984.192351700, 0.915087231),c(
   0.000402,  11371.704689758, 4.107281715),c(
   0.000453,   9917.696874510, 1.917490952),c(
   0.000395,    149.563197135, 2.763124165),c(
   0.000371,   5739.157790895, 3.112111866),c(
   0.000350,  11790.629088659, 0.440639857),c(
   0.000356,   6133.512652857, 5.444568842),c(
   0.000344,    412.371096874, 5.676832684),c(
   0.000383,    955.599741609, 5.559734846),c(
   0.000333,   6496.374945429, 0.261537984),c(
   0.000340,   6055.549660552, 5.975534987),c(
   0.000334,   1066.495477190, 2.335063907),c(
   0.000399,  11506.769769794, 5.321230910),c(
   0.000314,  18319.536584880, 2.313312404),c(
   0.000424,   1052.268383188, 1.211961766 ))
fbldata = rbind( fbldata,c(
   0.000307,     63.735898303, 3.169551388),c(
   0.000329,     29.821438149, 6.106912080),c(
   0.000357,   6309.374169791, 4.223760346),c(
   0.000312,  -3738.761430108, 2.180556645),c(
   0.000301,    309.278322656, 1.499984572),c(
   0.000268,  12043.574281889, 2.447520648),c(
   0.000257,  12491.370101415, 3.662331761),c(
   0.000290,    625.670192312, 1.272834584),c(
   0.000256,   5429.879468239, 1.913426912),c(
   0.000339,   3496.032826134, 4.165930011),c(
   0.000283,   3930.209696220, 4.325565754),c(
   0.000241,  12528.018664345, 3.832324536),c(
   0.000304,   4686.889407707, 1.612348468),c(
   0.000259,  16200.772724501, 3.470173146),c(
   0.000238,  12139.553509107, 1.147977842),c(
   0.000236,   6172.869528772, 3.776271728),c(
   0.000296,  -7058.598461315, 0.460368852),c(
   0.000306,  10575.406682942, 0.554749016),c(
   0.000251,  17298.182327326, 0.834332510),c(
   0.000290,   4732.030627343, 4.759564091),c(
   0.000261,   5884.926846583, 0.298259862),c(
   0.000249,   5547.199336460, 3.749366406),c(
   0.000213,  11712.955318231, 5.415666119),c(
   0.000223,   4701.116501708, 2.703203558),c(
   0.000268,   -640.877607382, 0.283670793 ))
fbldata = rbind( fbldata,c(
   0.000209,   5636.065016677, 1.238477199),c(
   0.000193,  10177.257679534, 1.943251340),c(
   0.000182,   6283.143160294, 2.456157599),c(
   0.000184,   -227.526189440, 5.888038582),c(
   0.000182,  -6283.008539689, 0.241332086),c(
   0.000228,  -6284.056171060, 2.657323816),c(
   0.000166,   7238.675591600, 5.930629110),c(
   0.000167,   3097.883822726, 5.570955333),c(
   0.000159,   -323.505416657, 5.786670700),c(
   0.000154,  -4136.910433516, 1.517805532),c(
   0.000176,  12029.347187887, 3.139266834),c(
   0.000167,  12132.439962106, 3.556352289),c(
   0.000153,    202.253395174, 1.463313961),c(
   0.000157,  17267.268201691, 1.586837396),c(
   0.000142,  83996.847317911, 0.022670115),c(
   0.000152,  17260.154654690, 0.708528947),c(
   0.000144,   6084.003848555, 5.187075177),c(
   0.000135,   5756.566278634, 1.993229262),c(
   0.000134,   5750.203491159, 3.457197134),c(
   0.000144,   5326.786694021, 6.066193291),c(
   0.000160,  11015.106477335, 1.710431974),c(
   0.000133,   3634.621024518, 2.836451652),c(
   0.000134,  18073.704938650, 5.453106665),c(
   0.000134,   1162.474704408, 5.326898811),c(
   0.000128,   5642.198242609, 2.511652591 ))
fbldata = rbind( fbldata,c(
   0.000160,    632.783739313, 5.628785365),c(
   0.000132,  13916.019109642, 0.819294053),c(
   0.000122,  14314.168113050, 5.677408071),c(
   0.000125,  12359.966151546, 5.251984735),c(
   0.000121,   5749.452731634, 2.210924603),c(
   0.000136,   -245.831646229, 1.646502367),c(
   0.000120,   5757.317038160, 3.240883049),c(
   0.000134,  12146.667056108, 3.059480037),c(
   0.000137,   6206.809778716, 1.867105418),c(
   0.000141,  17253.041107690, 2.069217456),c(
   0.000129,  -7477.522860216, 2.781469314),c(
   0.000116,   5540.085789459, 4.281176991),c(
   0.000116,   9779.108676125, 3.320925381),c(
   0.000129,   5237.921013804, 3.497704076),c(
   0.000113,   5959.570433334, 0.983210840),c(
   0.000122,   6282.095528923, 2.674938860),c(
   0.000140,    -11.045700264, 4.957936982),c(
   0.000108,  23543.230504682, 1.390113589),c(
   0.000106, -12569.674818332, 0.429631317),c(
   0.000110,   -266.607041722, 5.501340197),c(
   0.000115,  12559.038152982, 4.691456618),c(
   0.000134,  -2388.894020449, 0.577313584),c(
   0.000109,  10440.274292604, 6.218148717),c(
   0.000102,   -543.918059096, 1.477842615),c(
   0.000108,  21228.392023546, 2.237753948 ))
fbldata = rbind( fbldata,c(
   0.000101,  -4535.059436924, 3.100492232),c(
   0.000103,     76.266071276, 5.594294322),c(
   0.000104,    949.175608970, 5.674287810),c(
   0.000101,  13517.870106233, 2.196632348),c(
   0.000100,  11933.367960670, 4.056084160 ))

i2terms = nrow(fbldata)
## T**2
fbldata = rbind( fbldata,c(
   4.322990,   6283.075849991, 2.642893748),c(
   0.406495,      0.000000000, 4.712388980),c(
   0.122605,  12566.151699983, 2.438140634),c(
   0.019476,    213.299095438, 1.642186981),c(
   0.016916,    529.690965095, 4.510959344),c(
   0.013374,     -3.523118349, 1.502210314),c(
   0.008042,     26.298319800, 0.478549024),c(
   0.007824,    155.420399434, 5.254710405),c(
   0.004894,   5746.271337896, 4.683210850),c(
   0.004875,   5760.498431898, 0.759507698),c(
   0.004416,   5223.693919802, 6.028853166),c(
   0.004088,     -7.113547001, 0.060926389),c(
   0.004433,  77713.771467920, 3.627734103),c(
   0.003277,  18849.227549974, 2.327912542),c(
   0.002703,   6062.663207553, 1.271941729),c(
   0.003435,   -775.522611324, 0.747446224),c(
   0.002618,   6076.890301554, 3.633715689),c(
   0.003146,    206.185548437, 5.647874613),c(
   0.002544,   1577.343542448, 6.232904270),c(
   0.002218,   -220.412642439, 1.309509946),c(
   0.002197,   5856.477659115, 2.407212349),c(
   0.002897,   5753.384884897, 5.863842246),c(
   0.001766,    426.598190876, 0.754113147),c(
   0.001738,   -796.298006816, 2.714942671),c(
   0.001695,    522.577418094, 2.629369842 ))
fbldata = rbind( fbldata,c(
   0.001584,   5507.553238667, 1.341138229),c(
   0.001503,   -242.728603974, 0.377699736),c(
   0.001552,   -536.804512095, 2.904684667),c(
   0.001370,   -398.149003408, 1.265599125),c(
   0.001889,  -5573.142801634, 4.413514859),c(
   0.001722,   6069.776754553, 2.445966339),c(
   0.001124,   1059.381930189, 5.041799657),c(
   0.001258,    553.569402842, 3.849557278),c(
   0.000831,    951.718406251, 2.471094709),c(
   0.000767,   4694.002954708, 5.363125422),c(
   0.000756,   1349.867409659, 1.046195744),c(
   0.000775,    -11.045700264, 0.245548001),c(
   0.000597,   2146.165416475, 4.543268798),c(
   0.000568,   5216.580372801, 4.178853144),c(
   0.000711,   1748.016413067, 5.934271972),c(
   0.000499,  12036.460734888, 0.624434410),c(
   0.000671,  -1194.447010225, 4.136047594),c(
   0.000488,   5849.364112115, 2.209679987),c(
   0.000621,   6438.496249426, 4.518860804),c(
   0.000495,  -6286.598968340, 1.868201275),c(
   0.000456,   5230.807466803, 1.271231591),c(
   0.000451,   5088.628839767, 0.084060889),c(
   0.000435,   5643.178563677, 3.324456609),c(
   0.000387,  10977.078804699, 4.052488477),c(
   0.000547, 161000.685737473, 2.841633844 ))
fbldata = rbind( fbldata,c(
   0.000522,   3154.687084896, 2.171979966),c(
   0.000375,   5486.777843175, 4.983027306),c(
   0.000421,   5863.591206116, 4.546432249),c(
   0.000439,   7084.896781115, 0.522967921),c(
   0.000309,   2544.314419883, 3.172606705),c(
   0.000347,   4690.479836359, 1.479586566),c(
   0.000317,    801.820931124, 3.553088096),c(
   0.000262,    419.484643875, 0.606635550),c(
   0.000248,   6836.645252834, 3.014082064),c(
   0.000245,  -1592.596013633, 5.519526220),c(
   0.000225,   4292.330832950, 2.877956536),c(
   0.000214,   7234.794256242, 1.605227587),c(
   0.000205,   5767.611978898, 0.625804796),c(
   0.000180,  10447.387839604, 3.499954526),c(
   0.000229,    199.072001436, 5.632304604),c(
   0.000214,    639.897286314, 5.960227667),c(
   0.000175,   -433.711737877, 2.162417992),c(
   0.000209,    515.463871093, 2.322150893),c(
   0.000173,   6040.347246017, 2.556183691),c(
   0.000184,   6309.374169791, 4.732296790),c(
   0.000227, 149854.400134205, 5.385812217),c(
   0.000154,   8031.092263058, 5.120720920),c(
   0.000151,   5739.157790895, 4.815000443),c(
   0.000197,   7632.943259650, 0.222827271),c(
   0.000197,     74.781598567, 3.910456770 ))
fbldata = rbind( fbldata,c(
   0.000138,   6055.549660552, 1.397484253),c(
   0.000149,  -6127.655450557, 5.333727496),c(
   0.000137,   3894.181829542, 4.281749907),c(
   0.000135,   9437.762934887, 5.979971885),c(
   0.000139,  -2352.866153772, 4.715630782),c(
   0.000142,   6812.766815086, 0.513330157),c(
   0.000120,  -4705.732307544, 0.194160689),c(
   0.000131, -71430.695617928, 0.000379226),c(
   0.000124,   6279.552731642, 2.122264908),c(
   0.000108,  -6256.777530192, 0.883445696 ))

i3terms = nrow(fbldata)
## T**3
fbldata = rbind( fbldata,c(
   0.143388,   6283.075849991, 1.131453581),c(
   0.006671,  12566.151699983, 0.775148887),c(
   0.001480,    155.420399434, 0.480016880),c(
   0.000934,    213.299095438, 6.144453084),c(
   0.000795,    529.690965095, 2.941595619),c(
   0.000673,   5746.271337896, 0.120415406),c(
   0.000672,   5760.498431898, 5.317009738),c(
   0.000389,   -220.412642439, 3.090323467),c(
   0.000373,   6062.663207553, 3.003551964),c(
   0.000360,   6076.890301554, 1.918913041),c(
   0.000316,    -21.340641002, 5.545798121),c(
   0.000315,   -242.728603974, 1.884932563),c(
   0.000278,    206.185548437, 1.266254859),c(
   0.000238,   -536.804512095, 4.532664830),c(
   0.000185,    522.577418094, 4.578313856),c(
   0.000245,  18849.227549974, 0.587467082),c(
   0.000180,    426.598190876, 5.151178553),c(
   0.000200,    553.569402842, 5.355983739),c(
   0.000141,   5223.693919802, 1.336556009),c(
   0.000104,   5856.477659115, 4.239842759 ))

i4terms = nrow(fbldata)/3
## T**4
fbldata = rbind( fbldata,c(
   0.003826,   6283.075849991, 5.705257275),c(
   0.000303,  12566.151699983, 5.407132842),c(
   0.000209,    155.420399434, 1.989815753 ))

nterms = nrow(fbldata)
const0 = fbldata[,1]
freq0  = fbldata[,2]
phase0 = fbldata[,3]

texp = rep(0,nterms)
texp[(i1terms+1):i2terms] = 1
texp[(i2terms+1):i3terms] = 2
texp[(i3terms+1):i4terms] = 3
texp[(i4terms+1):length(texp)] = 4

dt <- z <- rep(NA,nrow(jd))
for(j in 1:nrow(jd)){
  t = ((jd[j,1]-2451545) + jd[j,2])/365250.0
  if(t==0) t <-  1e-50

  ph = freq0 * t + phase0
  sint = sin( ph )
  sinf = const0 * t^texp

  dt[j] = sum(sinf*sint)*1e-6
  z[j] = sum(sinf*(texp*sint/t + freq0*cos(ph)))*(1e-6/365250.0/DAYSEC)
}
return(cbind(dt,z))
}
