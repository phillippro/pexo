c*******************************************************************************
cStart of the header for a fortran source file for a subroutine
cof the Irwin pulsar epoch correction and precise radial velocity
c       barycentric code
cCopyright (C) 2003 by Alan W. Irwin
c
c$Id: fortran.header,v 2.3 1996/08/20 21:39:10 eos Exp $
c
cFor the latest version of this source code, please contact
cAlan W. Irwin
cDepartment of Physics and Astronomy
cUniversity of Victoria,
cBox 3055
cVictoria, B.C., Canada
cV8W 3P6
ce-mail: irwin@uvastro.phys.uvic.ca.
c
c    This program is free software; you can redistribute it and/or modify
c    it under the terms of the GNU General Public License as published by
c    the Free Software Foundation; either version 2 of the License, or
c    (at your option) any later version.
c
c    This program is distributed in the hope that it will be useful,
c    but WITHOUT ANY WARRANTY; without even the implied warranty of
c    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c    GNU General Public License for more details.
c
c    You should have received a copy of the GNU General Public License
c    along with this program; if not, write to the Free Software
c    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
c
cEnd of the header for a fortran source file for a subroutine
cof the Irwin epoch and barycentric correction code
c*******************************************************************************
      program fb2001
      implicit none
      real*8 ampcut, value, value_t0, deriv
      real*8 t0, teph(2)
      real*8 tstart, tstop, dt
      logical iftcb
      integer*4 iday, maxday
C      Input epoch range
      read(*,*) tstart, tstop, dt, iftcb
      if(dt.le.0.d0) stop 'dt must be positive'
      if(tstart.ge.tstop) stop 'tstart must be less than tstop'
      maxday = int((tstop-tstart)/dt)
C       ampcut in microsecs so 10 nanosecs cutoff = 0.010.
C       use full series
      ampcut = 0.d0
C       epoch where all time scales agree
      teph(1) = 2443144.5d0
      teph(2) = 0.d0
      call fbext(iftcb, ampcut, teph, value, deriv)
      teph(2) = 0.0003725d0
      call fbext(iftcb, ampcut, teph, value_t0, deriv)
      write(*,'(a,1p2d24.16)')
     &  'time ephemeris and derivative @ t0 =',
     &  0.d0, deriv
      write(*,'(a,1p2d24.16)')
     &  'integral from midnight to t0 and derivative =',
     &  value_t0-value, deriv
      write(*,'(a,1p2d24.16)') 'epoch range = ',
     &  tstart, tstop
      t0 = tstart
      teph(2) = 0.d0
      do iday = 0,maxday
        teph(1) = t0 + dble(iday)*dt
        call fbext(iftcb, ampcut, teph, value, deriv)
        value = value - value_t0
        write(*,'(f9.1,1p2d24.16)') teph(1), value, deriv
      enddo

      end
      subroutine fbext(iftcb, ampcut, teph, result, deriv)
C       routine uses series with coefficients in standard form
C       read in from stdin by serfbext_read.
C       iftcb .true. means series is in tcb so must transform teph to tcb
C       If the amplitude (microsecs)
C       is greater than or equal to the amplitude cutoff ampcut (microsecs)
C       the term is used as part of the series which is used to calculate
C       fbext.
C       teph is the julian date (teph) in days.
C       result (days) is 3rd term in eq. 3 for Teph - TT.
C       deriv (dimensionless) is the derivative of result wrt teph.
C       the equation for result appears in Standish, 1998 
C       and is calculated without the multiplication by K.
C       following jplephem, teph is split into two components.  For best
C       precision teph(1) should be half integral
      implicit none
      logical iftcb
      real*8 ampcut
      real*8 teph(2)
      integer*4 maxterms, maxorder, iorder, iterm
      parameter (maxterms = 15000, maxorder=10)
      integer*4 nterms(maxorder), norder
      real*8 series(3, maxterms, maxorder), t, result, deriv
      integer*4 iffirst
      data iffirst /1/
      real*8 t0, kmone, teph0mt0
C       t0 and best estimates of K and teph0-t0 are given in Irwin&Fukushima 1999.
      parameter (t0 = 2 443 144.500 372 5d0)
      parameter (kmone = 1.550 519 791 54d-8)
      parameter (teph0mt0 = -65.564 518d-6/86400.d0)
      save
      if(iffirst.eq.1) then
        iffirst = 0
        call fbext_read(
     &    series, nterms, norder, maxterms, maxorder)
      endif
C       t is in kilo julian years since J2000
      if(iftcb) then
C         tcb-j2000 = K (teph-j2000) + (1-K) (t0-j2000) - K (teph0-t0)
        t = ((1.d0 + kmone)*
     &    ((teph(1)-2451545.d0) + teph(2) - teph0mt0) -
     &    kmone*teph0mt0)/365250.d0
      else
        t=((teph(1)-2451545.d0) + teph(2))/365250.d0
      endif
      result = 0.d0
      deriv = 0.d0
      do iorder = norder,1, -1
C         programme differentiation checked to give correct result
        deriv = result + t*deriv
        result = t*result
        do iterm = nterms(iorder),1,-1
          if(abs(series(1,iterm, iorder)).ge.ampcut) then
            result = result + series(1, iterm, iorder)*
     &        sin(series(2, iterm, iorder)*t + series(3, iterm, iorder))
            deriv = deriv +
     &        series(1, iterm, iorder)*series(2, iterm, iorder)*
     &        cos(series(2, iterm, iorder)*t + series(3, iterm, iorder))
          endif
        enddo
      enddo
C       result in microseconds and derivative of that result wrt t
C       convert from microseconds to days and t derivative to teph derivative
      result = result/86400.d6
      deriv = deriv/(86400.d6*365250.d0)
      if(iftcb) then
        deriv = deriv*(1.d0+kmone)
      endif
      end
      subroutine fbext_read(
     &  series, nterms, norder, maxterms, maxorder)
C       read in TCB-TT series coefficients from fb2001.dat
C       n.b. the data must be in a specified (see below)
C       sorted order, otherwise this routine fails.
C       series(nterms(norder), norder), series coefficients
C       nterms(maxorder) = number of terms for a given order
C       norder = number of orders
C       maxterms = (input) maximum number of terms for a given order
C       maxorder = (input) maximum number of orders.
      implicit none
      integer*4 maxorder, nterms(maxorder), norder, norderm1,
     &  maxterms, ifeof,
     &  iterm, iterm_old, idum, norder_old, index(12)
      real*8 series(3, maxterms, maxorder), serin(3)
      save
      open(1, file='fb2001.dat', form='formatted', status='old')
      iterm = 1
      read(1,'(2x, 2i3, f13.6, f18.9, f13.9, 1x, 3(2x,4i3))',
     &  iostat=ifeof) idum, norderm1, serin, index
      norder_old = norderm1+1
C       first record must be of specific form
      if(iterm.ne.1.or.norderm1.ne.0)
     &  stop 'fbext_read: first order must be zero'
      do while(ifeof.eq.0)
        norder = norderm1+1
C         checks.
        if(norder.gt.maxorder) stop 'fbext_read: norder too large'
        if(iterm.gt.maxterms) stop 'fbext_read: iterm too large'
        if(norder.gt.norder_old+1.or.norder.lt.norder_old) then
          stop 'fbext_read: file must only ascend by 1 in order'
        elseif(norder.eq.norder_old+1) then
C           this can only happen for norder = 2 or higher
          norder_old = norder
          nterms(norder-1) = iterm_old
          iterm = 1
        else
C           same order
        endif
        series(1,iterm,norder) = serin(1)
        series(2,iterm,norder) = serin(2)
        series(3,iterm,norder) = serin(3)
        iterm_old = iterm
        iterm = iterm + 1
        read(1,'(2x, 2i3, f13.6, f18.9, f13.9, 1x, 3(2x,4i3))',
     &    iostat=ifeof) idum, norderm1, serin, index
      enddo
C       finish nterms
      nterms(norder) = iterm_old
      end
