	program xhf2002
c	test driver of 'hf2002', a routine computing approximated TCB-TCG for TE405
	real*8 t0,tjd,x
	integer iter,jyear
	real*8 hf2002
c
      t0=2451545.d0
      write(*,100)'Year','JD','TCB-TCG (sec)'
  100 format(t2,a,t12,a,t20,a)
c
      do jyear=1600,2200,100
		tjd=t0+dble(jyear-2000)*365.25
		x=hf2002(tjd)
		write(*,200) jyear,tjd,x
  200		format(i5,f10.1,f21.13)
	enddo
      end
c
      real*8 function hf2002(tjd)
c
c	Routine to provide TCB-TCG in seconds (Harada and Fukushima 2002)
c
c	The routine approximates TE405 (including the trend)
c	with errors of 0.453 ns (RMS) and 2.248 ns (Max) for the period 1600-2200
c
c	Note that an offset '0.3960 ns' remains at the epoch JD=2443144.5003725
c
      real*8 tjd
c
	integer ne,nx,i,j
	real*8 offset,a1,a2,f,t_start,pi2,cN,cNdti,dt
	real*8 t,pi2t,xi,s0,s1,s2,arg
	logical first/.true./
      parameter(n=73000,ne=463,nx=36)
c  ne: number of Fourier (=sin,cos) terms
c  nx: number of mixed secular (=t*sin,t*cos) terms
      common /hf2002c/a1(ne+nx+1),a2(ne+nx+1),f(ne+nx+1),t_start
c  a1: coefficient of sin or t*sin
c  a2: coefficient of cos or t*cos
c  f : frequency of sin/cos or t*sin/t*cos terms
c
	if(first) then
		first=.false.
		open(11,file='hf2002.dat')
		read(11,110)(a1(i),a2(i),f(i),i=1,ne+nx+1)
  110		format(3e24.16)
		close(11)
c	    offset=3.960d-10
	    offset=0.d-10
		t_start=2305450.5d0
		pi2=8.d0*datan(1.d0)
		cN1dt=3.d0*dble(n-1)
		cNdti=4.d0/cN1dt
	endif
c
      t=(tjd-t_start)-0.5d0*cN1dt
	pi2t=pi2*t
	xi=cNdti*t
c
	s0=(a1(1)-offset)+xi*(a2(1)+f(1)*xi)
c	Fourier part
	s1=0.d0
	do i=1,ne
		j=ne+2-i
		arg=pi2t*f(j)
		s1=s1+(a1(j)*dsin(arg)+a2(j)*dcos(arg))
	enddo
c	Mixed secular part
	s2=0.d0
      do i=1,nx
         j=ne+nx+2-i
	   arg=pi2t*f(j)
         s2=s2+xi*(a1(j)*dsin(arg)+a2(j)*dcos(arg))
	enddo
c
	hf2002=(s2+s1)+s0
c
	return
      end
