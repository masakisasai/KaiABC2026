  module param

  real(kind=8),parameter::pi=4.0d0*ATAN(1.0d0)
  real(kind=8),parameter::tpi=8.0d0*ATAN(1.0d0)
  
  integer,parameter::ntmax=10**7             !simulation time length
  integer,parameter::interval=10             !interval steps
  integer,parameter::iwarm=2*10**6           !warming up
  integer,parameter::Nvar=10**6
  integer,parameter::ninterval=ntmax/interval
  
  real(kind=8),parameter::dt=1.0d0/10**3     !time step  
  integer,parameter::Nc=100
  
  real(kind=8),parameter::op=tpi/12.0d0     
  real(kind=8),parameter::bt=0.50d0
  real(kind=8),parameter::bp=0.50d0
  character(len=*), parameter :: ofile = 'amplitude050.dat'
  
  real(kind=8),parameter::Difft=0.01d0
  real(kind=8),parameter::Diffp=0.01d0 

  integer,parameter::iprint=0
  
  integer,parameter ::nslide=5  
  integer,parameter ::nsamp=ninterval/nslide  !amplitude evaluation
  integer,parameter ::nwidth=ninterval/nslide

  integer,parameter ::Nh=19      !FFT parameter
  integer,parameter ::N2=2**Nh
  integer,parameter ::nmaxsqrt = 2**Nh
  integer,parameter ::nkh2=N2/2-1
  
  real(kind=8),parameter::st=sqrt(2*Difft/dt)  
  real(kind=8),parameter::sp=sqrt(2*Diffp/dt)
  real(kind=8),parameter::eop=(bt/(bt+bp))*op*dt
  
  end module param

  !#################################################################################

  program main
    use param
    
    integer::i,j,k,l,m,n,it,iit,ia,ib,ic   
    integer::seed
    real(kind=8)::r8_normal_01
    real(kind=8)::r8_uniform_01
    real(kind=8)::xa,xc 
    real(kind=8)::time,beta,betaD
    real(kind=8)::dxp,dxt,dmxp,dmxt,axp,axt,a2xp,a2xt,delp,delt   
    real(kind=8),dimension(1:Nc)::xt,xp,sxt,cxp
    real(kind=8),dimension(1:ninterval)::ensxt,ensxp


    real(kind=8)::amax,amin,ampp,ampt,xmax,xmin
    real(kind=8)::cxtp,acxtp,betaB,ampB
    integer::nmax,nmin
    real(kind=8),dimension(1:nsamp)::vsamp        
    open(unit=64, file=ofile)
    
    seed=273417            
    do i=1,100
       f1=r8_uniform_01(seed)
    enddo

 do ib=1,50
       
    do i=1,Nc
      xt(i)=2*pi*r8_uniform_01(seed)
      xp(i)=2*pi*r8_uniform_01(seed)
    enddo

    beta=0.0d0+ib*0.001d0
    betaD=beta/Diffp 
    do it=1,iwarm
       call langevinstep(xt,xp,beta,seed)
    enddo

    iit=0
    do it=1,ntmax
       call langevinstep(xt,xp,beta,seed)
       if(mod(it,interval)==0) then
          iit=iit+1
          time=it*dt
          ensxt(iit)=0.0d0
          ensxp(iit)=0.0d0
          do i=1,Nc
             sxt(i)=sin(xt(i))
             cxp(i)=cos(xp(i))
             ensxt(iit)=ensxt(iit)+sxt(i)
             ensxp(iit)=ensxp(iit)+cxp(i)
          enddo
          ensxt(iit)=ensxt(iit)/real(Nc)
          ensxp(iit)=ensxp(iit)/real(Nc)
          if(iprint==1.and.mod(it,interval*10000)==0) write(6,*) sxt(1),cxp(1),time         
       endif
    enddo    

!evaluating amplitudes
    
    amax=0.0d0
    amin=0.0d0
    do j=1,nslide
       do i=1,nsamp
          vsamp(i)=ensxt(i+nwidth*(j-1))
       enddo
       call max0(vsamp,nmax,xmax)
       amax=amax+xmax
    enddo
    amax=amax/real(nslide)
    do j=1,nslide
       imin=imin+1
       do i=1,nsamp
          vsamp(i)=ensxt(i+nwidth*(j-1))
       enddo
       call min0(vsamp,nmin,xmin)
       amin=amin+xmin
    enddo
    amin=amin/real(nslide)
    ampt=amax-amin

    
    amax=0.0d0
    amin=0.0d0
    do j=1,nslide
       do i=1,nsamp
          vsamp(i)=ensxp(i+nwidth*(j-1))
       enddo
       call max0(vsamp,nmax,xmax)
       amax=amax+xmax
    enddo
    amax=amax/real(nslide)   
    do j=1,nslide
       imin=imin+1
       do i=1,nsamp
          vsamp(i)=ensxp(i+nwidth*(j-1))
       enddo
       call min0(vsamp,nmin,xmin)
       amin=amin+xmin
    enddo
    amin=amin/real(nslide)
    ampp=amax-amin

    
    if(iprint==1) write(6,*) 'amplitude',ampt,ampp
    write(64,*) beta/Diffp,bp/Diffp,ampt*0.5d0,ampp*0.5d0
    flush(64)

enddo


      
    stop
  end program main
  
!****************************************************************************
  
  subroutine langevinstep(xt,xp,beta,seed)
  use param

    integer::i,j,k,l,m  
    integer::seed
    real(kind=8)::r8_normal_01
    real(kind=8)::r8_uniform_01 
    real(kind=8),dimension(1:Nc)::xt,xp
    real(kind=8)::xa0,rcp,rsp,at,ap,ft,fp,sxtp,cxtp,beta

    rcp=0.0d0
    rsp=0.0d0
    do i=1,Nc
       rcp=rcp+cos(xp(i))
       rsp=rsp+sin(xp(i))
    enddo
    rcp=rcp/real(Nc)
    rsp=rsp/real(Nc)

    do i=1,Nc
       sxtp=sin(xt(i)-xp(i))
       cxtp=cos(xt(i)-xp(i))
       at=st*r8_normal_01(seed)    
       ap=sp*r8_normal_01(seed)
     
       ft=-bt*sxtp+at
       fp=op +bp*sxtp -beta*(sin(xp(i))*rcp-cos(xp(i))*rsp)*cxtp +ap
       
       xt(i)=xt(i)+dt*ft
       xp(i)=xp(i)+dt*fp
    enddo

   
       
return
end subroutine langevinstep
!************************************************************************

subroutine  max(v,kpeak)
use param

  real(kind=8),dimension(1:nkh2)::v
  real(kind=8)::xmax
  integer::i,kpeak


  xmax=v(1)
  kpeak=0
  do i =1,nkh2
     if(xmax < v(i+1))then
        xmax=v(i+1)
        kpeak=i+1
     else if(xmax > v(i+1))then
        xmax=xmax
     endif
  end do


return
end subroutine max

!************************************************************************

subroutine  min0(v,nmin,xmin)
use param

  real(kind=8),dimension(1:nsamp)::v
  real(kind=8)::xmin
  integer::i,nmin


  xmin=v(1)
  nmin=1
  do i =1,nsamp-1
     if(xmin > v(i+1))then
        xmin=v(i+1)
        nmin=i+1
     else if(xmin < v(i+1))then
        xmin=xmin
     endif
  end do


return
end subroutine min0


!************************************************************************

subroutine  max0(v,nmax,xmax)
use param

  real(kind=8),dimension(1:nsamp)::v
  real(kind=8)::xmax
  integer::i,nmax


  xmax=v(1)
  nmax=1
  do i =1,nsamp-1
     if(xmax < v(i+1))then
        xmax=v(i+1)
        nmax=i+1
     else if(xmax > v(i+1))then
        xmax=xmax
     endif
  end do


return
end subroutine max0



!*****************************************************

subroutine rdft(n,isgn,a,ip,w)
      integer:: n, isgn, nw, nc
      integer,dimension(0 : *)::ip
      real(kind=8),dimension(0 : n - 1)::a
      real(kind=8),dimension(0 : *)::w
      real(kind=8)::xi
      
      nw = ip(0)
      if (n .gt. 4 * nw) then
          nw = n / 4
          call makewt(nw, ip, w)
      end if
      nc = ip(1)
      if (n .gt. 4 * nc) then
          nc = n / 4
          call makect(nc, ip, w(nw))
      end if
      if (isgn .ge. 0) then
          if (n .gt. 4) then
              call bitrv2(n, ip(2), a)
              call cftfsub(n, a, w)
              call rftfsub(n, a, nc, w(nw))
          else if (n .eq. 4) then
              call cftfsub(n, a, w)
          end if
          xi = a(0) - a(1)
          a(0) = a(0) + a(1)
          a(1) = xi
      else
          a(1) = 0.5d0 * (a(0) - a(1))
          a(0) = a(0) - a(1)
          if (n .gt. 4) then
              call rftbsub(n, a, nc, w(nw))
              call bitrv2(n, ip(2), a)
              call cftbsub(n, a, w)
          else if (n .eq. 4) then
              call cftfsub(n, a, w)
          end if
      end if
      end
!
!
! -------- initializing routines --------
!
   subroutine makewt(nw, ip, w)
     integer::nw,j,nwh
     integer,dimension(0 : *)::ip
     real(kind=8),dimension(0 : nw - 1)::w
     real(kind=8)::delta, x, y
     
      ip(0) = nw
      ip(1) = 1
      if (nw .gt. 2) then
          nwh = nw / 2
          delta = atan(1.0d0) / nwh
          w(0) = 1
          w(1) = 0
          w(nwh) = cos(delta * nwh)
          w(nwh + 1) = w(nwh)
          if (nwh .gt. 2) then
              do j = 2, nwh - 2, 2
                  x = cos(delta * j)
                  y = sin(delta * j)
                  w(j) = x
                  w(j + 1) = y
                  w(nw - j) = y
                  w(nw - j + 1) = x
              end do
              call bitrv2(nw, ip(2), w)
          end if
      end if
      end
!
  subroutine makect(nc, ip, c)
      integer::nc,j,nch
      integer,dimension(0 : *)::ip 
      real(kind=8),dimension(0 : nc - 1)::c
      real(kind=8)::delta
      
      ip(1) = nc
      if (nc .gt. 1) then
          nch = nc / 2
          delta = atan(1.0d0) / nch
          c(0) = cos(delta * nch)
          c(nch) = 0.5d0 * c(0)
          do j = 1, nch - 1
              c(j) = 0.5d0 * cos(delta * j)
              c(nc - j) = 0.5d0 * sin(delta * j)
          end do
      end if
      end
!
! -------- child routines --------
!
  subroutine bitrv2(n, ip, a)
      integer::n,j,j1,k,k1,l,m,m2
      integer,dimension(0 : *)::ip 
      real(kind=8),dimension(0 : n - 1)::a
      real(kind=8)::xr,xi,yr,yi

      
      ip(0) = 0
      l = n
      m = 1
      do while (8 * m .lt. l)
          l = l / 2
          do j = 0, m - 1
              ip(m + j) = ip(j) + l
          end do
          m = m * 2
      end do
      m2 = 2 * m
      if (8 * m .eq. l) then
          do k = 0, m - 1
              do j = 0, k - 1
                  j1 = 2 * j + ip(k)
                  k1 = 2 * k + ip(j)
                  xr = a(j1)
                  xi = a(j1 + 1)
                  yr = a(k1)
                  yi = a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + m2
                  k1 = k1 + 2 * m2
                  xr = a(j1)
                  xi = a(j1 + 1)
                  yr = a(k1)
                  yi = a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + m2
                  k1 = k1 - m2
                  xr = a(j1)
                  xi = a(j1 + 1)
                  yr = a(k1)
                  yi = a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + m2
                  k1 = k1 + 2 * m2
                  xr = a(j1)
                  xi = a(j1 + 1)
                  yr = a(k1)
                  yi = a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
              end do
              j1 = 2 * k + m2 + ip(k)
              k1 = j1 + m2
              xr = a(j1)
              xi = a(j1 + 1)
              yr = a(k1)
              yi = a(k1 + 1)
              a(j1) = yr
              a(j1 + 1) = yi
              a(k1) = xr
              a(k1 + 1) = xi
          end do
      else
          do k = 1, m - 1
              do j = 0, k - 1
                  j1 = 2 * j + ip(k)
                  k1 = 2 * k + ip(j)
                  xr = a(j1)
                  xi = a(j1 + 1)
                  yr = a(k1)
                  yi = a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + m2
                  k1 = k1 + m2
                  xr = a(j1)
                  xi = a(j1 + 1)
                  yr = a(k1)
                  yi = a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
              end do
          end do
      end if
      end
!
  subroutine bitrv2conj(n, ip, a)
      integer::n,j,j1,k,k1,l,m,m2
      integer,dimension(0 : *)::ip
      real(kind=8),dimension(0 : n - 1)::a
      real(kind=8)::xr,xi,yr,yi
      
      ip(0) = 0
      l = n
      m = 1
      do while (8 * m .lt. l)
          l = l / 2
          do j = 0, m - 1
              ip(m + j) = ip(j) + l
          end do
          m = m * 2
      end do
      m2 = 2 * m
      if (8 * m .eq. l) then
          do k = 0, m - 1
              do j = 0, k - 1
                  j1 = 2 * j + ip(k)
                  k1 = 2 * k + ip(j)
                  xr = a(j1)
                  xi = -a(j1 + 1)
                  yr = a(k1)
                  yi = -a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + m2
                  k1 = k1 + 2 * m2
                  xr = a(j1)
                  xi = -a(j1 + 1)
                  yr = a(k1)
                  yi = -a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + m2
                  k1 = k1 - m2
                  xr = a(j1)
                  xi = -a(j1 + 1)
                  yr = a(k1)
                  yi = -a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + m2
                  k1 = k1 + 2 * m2
                  xr = a(j1)
                  xi = -a(j1 + 1)
                  yr = a(k1)
                  yi = -a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
              end do
              k1 = 2 * k + ip(k)
              a(k1 + 1) = -a(k1 + 1)
              j1 = k1 + m2
              k1 = j1 + m2
              xr = a(j1)
              xi = -a(j1 + 1)
              yr = a(k1)
              yi = -a(k1 + 1)
              a(j1) = yr
              a(j1 + 1) = yi
              a(k1) = xr
              a(k1 + 1) = xi
              k1 = k1 + m2
              a(k1 + 1) = -a(k1 + 1)
          end do
      else
          a(1) = -a(1)
          a(m2 + 1) = -a(m2 + 1)
          do k = 1, m - 1
              do j = 0, k - 1
                  j1 = 2 * j + ip(k)
                  k1 = 2 * k + ip(j)
                  xr = a(j1)
                  xi = -a(j1 + 1)
                  yr = a(k1)
                  yi = -a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + m2
                  k1 = k1 + m2
                  xr = a(j1)
                  xi = -a(j1 + 1)
                  yr = a(k1)
                  yi = -a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
              end do
              k1 = 2 * k + ip(k)
              a(k1 + 1) = -a(k1 + 1)
              a(k1 + m2 + 1) = -a(k1 + m2 + 1)
          end do
      end if
      end
!
  subroutine cftfsub(n, a, w)
      integer::n,j,j1,j2,j3,l
      real(kind=8),dimension(0 : n - 1)::a
      real(kind=8),dimension(0 : *)::w
      real(kind=8)::x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i

      
      l = 2
      if (n .gt. 8) then
          call cft1st(n, a, w)
          l = 8
          do while (4 * l .lt. n)
              call cftmdl(n, l, a, w)
              l = 4 * l
          end do
      end if
      if (4 * l .eq. n) then
          do j = 0, l - 2, 2
              j1 = j + l
              j2 = j1 + l
              j3 = j2 + l
              x0r = a(j) + a(j1)
              x0i = a(j + 1) + a(j1 + 1)
              x1r = a(j) - a(j1)
              x1i = a(j + 1) - a(j1 + 1)
              x2r = a(j2) + a(j3)
              x2i = a(j2 + 1) + a(j3 + 1)
              x3r = a(j2) - a(j3)
              x3i = a(j2 + 1) - a(j3 + 1)
              a(j) = x0r + x2r
              a(j + 1) = x0i + x2i
              a(j2) = x0r - x2r
              a(j2 + 1) = x0i - x2i
              a(j1) = x1r - x3i
              a(j1 + 1) = x1i + x3r
              a(j3) = x1r + x3i
              a(j3 + 1) = x1i - x3r
          end do
      else
          do j = 0, l - 2, 2
              j1 = j + l
              x0r = a(j) - a(j1)
              x0i = a(j + 1) - a(j1 + 1)
              a(j) = a(j) + a(j1)
              a(j + 1) = a(j + 1) + a(j1 + 1)
              a(j1) = x0r
              a(j1 + 1) = x0i
          end do
      end if
      end
!
  subroutine cftbsub(n, a, w)
      integer::n,j,j1,j2,j3,l
      real(kind=8),dimension(0 : n - 1)::a
      real(kind=8),dimension(0 : *)::w
      real(kind=8)::x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i

      
      l = 2
      if (n .gt. 8) then
          call cft1st(n, a, w)
          l = 8
          do while (4 * l .lt. n)
              call cftmdl(n, l, a, w)
              l = 4 * l
          end do
      end if
      if (4 * l .eq. n) then
          do j = 0, l - 2, 2
              j1 = j + l
              j2 = j1 + l
              j3 = j2 + l
              x0r = a(j) + a(j1)
              x0i = -a(j + 1) - a(j1 + 1)
              x1r = a(j) - a(j1)
              x1i = -a(j + 1) + a(j1 + 1)
              x2r = a(j2) + a(j3)
              x2i = a(j2 + 1) + a(j3 + 1)
              x3r = a(j2) - a(j3)
              x3i = a(j2 + 1) - a(j3 + 1)
              a(j) = x0r + x2r
              a(j + 1) = x0i - x2i
              a(j2) = x0r - x2r
              a(j2 + 1) = x0i + x2i
              a(j1) = x1r - x3i
              a(j1 + 1) = x1i - x3r
              a(j3) = x1r + x3i
              a(j3 + 1) = x1i + x3r
          end do
      else
          do j = 0, l - 2, 2
              j1 = j + l
              x0r = a(j) - a(j1)
              x0i = -a(j + 1) + a(j1 + 1)
              a(j) = a(j) + a(j1)
              a(j + 1) = -a(j + 1) - a(j1 + 1)
              a(j1) = x0r
              a(j1 + 1) = x0i
          end do
      end if
      end
!
  subroutine cft1st(n, a, w)
      integer::n,j,k1,k2
      real(kind=8),dimension(0 : n - 1)::a
      real(kind=8),dimension(0 : *)::w
      real(kind=8)::wk1r,wk1i,wk2r,wk2i,wk3r,wk3i
      real(kind=8)::x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i

      
      x0r = a(0) + a(2)
      x0i = a(1) + a(3)
      x1r = a(0) - a(2)
      x1i = a(1) - a(3)
      x2r = a(4) + a(6)
      x2i = a(5) + a(7)
      x3r = a(4) - a(6)
      x3i = a(5) - a(7)
      a(0) = x0r + x2r
      a(1) = x0i + x2i
      a(4) = x0r - x2r
      a(5) = x0i - x2i
      a(2) = x1r - x3i
      a(3) = x1i + x3r
      a(6) = x1r + x3i
      a(7) = x1i - x3r
      wk1r = w(2)
      x0r = a(8) + a(10)
      x0i = a(9) + a(11)
      x1r = a(8) - a(10)
      x1i = a(9) - a(11)
      x2r = a(12) + a(14)
      x2i = a(13) + a(15)
      x3r = a(12) - a(14)
      x3i = a(13) - a(15)
      a(8) = x0r + x2r
      a(9) = x0i + x2i
      a(12) = x2i - x0i
      a(13) = x0r - x2r
      x0r = x1r - x3i
      x0i = x1i + x3r
      a(10) = wk1r * (x0r - x0i)
      a(11) = wk1r * (x0r + x0i)
      x0r = x3i + x1r
      x0i = x3r - x1i
      a(14) = wk1r * (x0i - x0r)
      a(15) = wk1r * (x0i + x0r)
      k1 = 0
      do j = 16, n - 16, 16
          k1 = k1 + 2
          k2 = 2 * k1
          wk2r = w(k1)
          wk2i = w(k1 + 1)
          wk1r = w(k2)
          wk1i = w(k2 + 1)
          wk3r = wk1r - 2 * wk2i * wk1i
          wk3i = 2 * wk2i * wk1r - wk1i
          x0r = a(j) + a(j + 2)
          x0i = a(j + 1) + a(j + 3)
          x1r = a(j) - a(j + 2)
          x1i = a(j + 1) - a(j + 3)
          x2r = a(j + 4) + a(j + 6)
          x2i = a(j + 5) + a(j + 7)
          x3r = a(j + 4) - a(j + 6)
          x3i = a(j + 5) - a(j + 7)
          a(j) = x0r + x2r
          a(j + 1) = x0i + x2i
          x0r = x0r - x2r
          x0i = x0i - x2i
          a(j + 4) = wk2r * x0r - wk2i * x0i
          a(j + 5) = wk2r * x0i + wk2i * x0r
          x0r = x1r - x3i
          x0i = x1i + x3r
          a(j + 2) = wk1r * x0r - wk1i * x0i
          a(j + 3) = wk1r * x0i + wk1i * x0r
          x0r = x1r + x3i
          x0i = x1i - x3r
          a(j + 6) = wk3r * x0r - wk3i * x0i
          a(j + 7) = wk3r * x0i + wk3i * x0r
          wk1r = w(k2 + 2)
          wk1i = w(k2 + 3)
          wk3r = wk1r - 2 * wk2r * wk1i
          wk3i = 2 * wk2r * wk1r - wk1i
          x0r = a(j + 8) + a(j + 10)
          x0i = a(j + 9) + a(j + 11)
          x1r = a(j + 8) - a(j + 10)
          x1i = a(j + 9) - a(j + 11)
          x2r = a(j + 12) + a(j + 14)
          x2i = a(j + 13) + a(j + 15)
          x3r = a(j + 12) - a(j + 14)
          x3i = a(j + 13) - a(j + 15)
          a(j + 8) = x0r + x2r
          a(j + 9) = x0i + x2i
          x0r = x0r - x2r
          x0i = x0i - x2i
          a(j + 12) = -wk2i * x0r - wk2r * x0i
          a(j + 13) = -wk2i * x0i + wk2r * x0r
          x0r = x1r - x3i
          x0i = x1i + x3r
          a(j + 10) = wk1r * x0r - wk1i * x0i
          a(j + 11) = wk1r * x0i + wk1i * x0r
          x0r = x1r + x3i
          x0i = x1i - x3r
          a(j + 14) = wk3r * x0r - wk3i * x0i
          a(j + 15) = wk3r * x0i + wk3i * x0r
      end do
      end
!
  subroutine cftmdl(n, l, a, w)
      integer::n,l,j,j1,j2,j3,k,k1,k2,m,m2
      real(kind=8),dimension(0 : n - 1)::a
      real(kind=8),dimension(0 : *)::w
      real(kind=8)::wk1r,wk1i,wk2r,wk2i,wk3r,wk3i
      real(kind=8)::x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i

      
      m = 4 * l
      do j = 0, l - 2, 2
          j1 = j + l
          j2 = j1 + l
          j3 = j2 + l
          x0r = a(j) + a(j1)
          x0i = a(j + 1) + a(j1 + 1)
          x1r = a(j) - a(j1)
          x1i = a(j + 1) - a(j1 + 1)
          x2r = a(j2) + a(j3)
          x2i = a(j2 + 1) + a(j3 + 1)
          x3r = a(j2) - a(j3)
          x3i = a(j2 + 1) - a(j3 + 1)
          a(j) = x0r + x2r
          a(j + 1) = x0i + x2i
          a(j2) = x0r - x2r
          a(j2 + 1) = x0i - x2i
          a(j1) = x1r - x3i
          a(j1 + 1) = x1i + x3r
          a(j3) = x1r + x3i
          a(j3 + 1) = x1i - x3r
      end do
      wk1r = w(2)
      do j = m, l + m - 2, 2
          j1 = j + l
          j2 = j1 + l
          j3 = j2 + l
          x0r = a(j) + a(j1)
          x0i = a(j + 1) + a(j1 + 1)
          x1r = a(j) - a(j1)
          x1i = a(j + 1) - a(j1 + 1)
          x2r = a(j2) + a(j3)
          x2i = a(j2 + 1) + a(j3 + 1)
          x3r = a(j2) - a(j3)
          x3i = a(j2 + 1) - a(j3 + 1)
          a(j) = x0r + x2r
          a(j + 1) = x0i + x2i
          a(j2) = x2i - x0i
          a(j2 + 1) = x0r - x2r
          x0r = x1r - x3i
          x0i = x1i + x3r
          a(j1) = wk1r * (x0r - x0i)
          a(j1 + 1) = wk1r * (x0r + x0i)
          x0r = x3i + x1r
          x0i = x3r - x1i
          a(j3) = wk1r * (x0i - x0r)
          a(j3 + 1) = wk1r * (x0i + x0r)
      end do
      k1 = 0
      m2 = 2 * m
      do k = m2, n - m2, m2
          k1 = k1 + 2
          k2 = 2 * k1
          wk2r = w(k1)
          wk2i = w(k1 + 1)
          wk1r = w(k2)
          wk1i = w(k2 + 1)
          wk3r = wk1r - 2 * wk2i * wk1i
          wk3i = 2 * wk2i * wk1r - wk1i
          do j = k, l + k - 2, 2
              j1 = j + l
              j2 = j1 + l
              j3 = j2 + l
              x0r = a(j) + a(j1)
              x0i = a(j + 1) + a(j1 + 1)
              x1r = a(j) - a(j1)
              x1i = a(j + 1) - a(j1 + 1)
              x2r = a(j2) + a(j3)
              x2i = a(j2 + 1) + a(j3 + 1)
              x3r = a(j2) - a(j3)
              x3i = a(j2 + 1) - a(j3 + 1)
              a(j) = x0r + x2r
              a(j + 1) = x0i + x2i
              x0r = x0r - x2r
              x0i = x0i - x2i
              a(j2) = wk2r * x0r - wk2i * x0i
              a(j2 + 1) = wk2r * x0i + wk2i * x0r
              x0r = x1r - x3i
              x0i = x1i + x3r
              a(j1) = wk1r * x0r - wk1i * x0i
              a(j1 + 1) = wk1r * x0i + wk1i * x0r
              x0r = x1r + x3i
              x0i = x1i - x3r
              a(j3) = wk3r * x0r - wk3i * x0i
              a(j3 + 1) = wk3r * x0i + wk3i * x0r
          end do
          wk1r = w(k2 + 2)
          wk1i = w(k2 + 3)
          wk3r = wk1r - 2 * wk2r * wk1i
          wk3i = 2 * wk2r * wk1r - wk1i
          do j = k + m, l + (k + m) - 2, 2
              j1 = j + l
              j2 = j1 + l
              j3 = j2 + l
              x0r = a(j) + a(j1)
              x0i = a(j + 1) + a(j1 + 1)
              x1r = a(j) - a(j1)
              x1i = a(j + 1) - a(j1 + 1)
              x2r = a(j2) + a(j3)
              x2i = a(j2 + 1) + a(j3 + 1)
              x3r = a(j2) - a(j3)
              x3i = a(j2 + 1) - a(j3 + 1)
              a(j) = x0r + x2r
              a(j + 1) = x0i + x2i
              x0r = x0r - x2r
              x0i = x0i - x2i
              a(j2) = -wk2i * x0r - wk2r * x0i
              a(j2 + 1) = -wk2i * x0i + wk2r * x0r
              x0r = x1r - x3i
              x0i = x1i + x3r
              a(j1) = wk1r * x0r - wk1i * x0i
              a(j1 + 1) = wk1r * x0i + wk1i * x0r
              x0r = x1r + x3i
              x0i = x1i - x3r
              a(j3) = wk3r * x0r - wk3i * x0i
              a(j3 + 1) = wk3r * x0i + wk3i * x0r
          end do
      end do
      end
!
      subroutine rftfsub(n, a, nc, c)
      integer n, nc, j, k, kk, ks, m
      real*8 a(0 : n - 1), c(0 : nc - 1), wkr, wki, xr, xi, yr, yi
      m = n / 2
      ks = 2 * nc / m
      kk = 0
      do j = 2, m - 2, 2
          k = n - j
          kk = kk + ks
          wkr = 0.5d0 - c(nc - kk)
          wki = c(kk)
          xr = a(j) - a(k)
          xi = a(j + 1) + a(k + 1)
          yr = wkr * xr - wki * xi
          yi = wkr * xi + wki * xr
          a(j) = a(j) - yr
          a(j + 1) = a(j + 1) - yi
          a(k) = a(k) + yr
          a(k + 1) = a(k + 1) - yi
      end do
      end
!
      subroutine rftbsub(n, a, nc, c)
      integer n, nc, j, k, kk, ks, m
      real*8 a(0 : n - 1), c(0 : nc - 1), wkr, wki, xr, xi, yr, yi
      a(1) = -a(1)
      m = n / 2
      ks = 2 * nc / m
      kk = 0
      do j = 2, m - 2, 2
          k = n - j
          kk = kk + ks
          wkr = 0.5d0 - c(nc - kk)
          wki = c(kk)
          xr = a(j) - a(k)
          xi = a(j + 1) + a(k + 1)
          yr = wkr * xr + wki * xi
          yi = wkr * xi - wki * xr
          a(j) = a(j) - yr
          a(j + 1) = yi - a(j + 1)
          a(k) = a(k) + yr
          a(k + 1) = yi - a(k + 1)
      end do
      a(m + 1) = -a(m + 1)
      end
!
      subroutine dctsub(n, a, nc, c)
      integer n, nc, j, k, kk, ks, m
      real*8 a(0 : n - 1), c(0 : nc - 1), wkr, wki, xr
      m = n / 2
      ks = nc / n
      kk = 0
      do j = 1, m - 1
          k = n - j
          kk = kk + ks
          wkr = c(kk) - c(nc - kk)
          wki = c(kk) + c(nc - kk)
          xr = wki * a(j) - wkr * a(k)
          a(j) = wkr * a(j) + wki * a(k)
          a(k) = xr
      end do
      a(m) = c(0) * a(m)
      end
!
  subroutine dstsub(n, a, nc, c)
      integer::n,nc,j,k,kk,ks,m
      real(kind=8),dimension(0 : n - 1)::a
      real(kind=8),dimension(0 : nc - 1)::c
      real(kind=8)::wkr,wki,xr

      
      m = n / 2
      ks = nc / n
      kk = 0
      do j = 1, m - 1
          k = n - j
          kk = kk + ks
          wkr = c(kk) - c(nc - kk)
          wki = c(kk) + c(nc - kk)
          xr = wki * a(k) - wkr * a(j)
          a(k) = wkr * a(k) + wki * a(j)
          a(j) = xr
      end do
      a(m) = c(0) * a(m)
      end
!




  
!************************************************************************
!
      function r8_normal_01 (seed)
!
!    The normal probability distribution function (PDF) is sampled,
!    with mean A and standard deviation B.

!    Input/output, integer SEED, a seed for the random number generator.
!    Output, double precision R8_NORMAL_AB, a sample of the normal PDF.

      implicit none

      double precision r1
      double precision r2
      double precision r8_normal_01
      double precision r8_pi
      parameter ( r8_pi = 3.141592653589793D+00 )
      double precision r8_uniform_01
      integer seed
      double precision x

      r1 = r8_uniform_01 ( seed )
      r2 = r8_uniform_01 ( seed )
      r8_normal_01 = sqrt ( -2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * r8_pi * r2 )


      return
      end

!*********************************************************************
!
      function r8_uniform_01 ( seed )
!
!  R8_UNIFORM_01 returns a unit pseudorandom R8.

!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2^31 - 1 )
!      r8_uniform_01 = seed / ( 2^31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702

!  Parameters:
!
!    Input/output, integer SEED, the "seed" value, which should NOT be 0.
!    On output, SEED has been updated.
!
!    Output, double precision R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
        

        
      implicit none

      integer i4_huge 
      parameter ( i4_huge = 2147483647 )
      integer k
      double precision r8_uniform_01
      integer seed

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + i4_huge
      end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!

      r8_uniform_01 = dble ( seed ) * 4.656612875D-10


      return
      end
    


