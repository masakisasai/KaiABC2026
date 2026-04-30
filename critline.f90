  module param

  real(kind=8),parameter::pi=4.0d0*ATAN(1.0d0)
  real(kind=8),parameter::tpi=2.0d0*pi
  
  integer,parameter::kfile=50
  integer,parameter::kstart=10
  integer,parameter::nbin=10
  real(kind=8),parameter::width=2.0d0
  
  integer,parameter ::Ndx=500
  real(kind=8),parameter::rdx=pi/real(Ndx)
  integer,parameter ::kb=1000
  integer,parameter ::kbeta=1000
  
  real(kind=8),parameter::op=tpi/12.0d0        
  real(kind=8),parameter::Difft=0.01d0
  real(kind=8),parameter::Diffp=0.01d0
  real(kind=8),parameter::Diffx=Difft+Diffp 
  real(kind=8),parameter::offx=op/Diffx  

  
  end module param

  !#################################################################################

  program main
    use param
    
    integer::i,j,k,l,m,n,ik 
    real(kind=8)::betaD,bD,res,bxD,xip,gxx,beta0
    real(kind=8),dimension(kb,kbeta)::betac
    real(kind=8),dimension(kbeta)::betax    
    real(kind=8),dimension(kb)::cxtp,bdy,zpcx   
    real(kind=8),dimension(-Ndx+1:Ndx-1)::Pcx   
    open(unit=70, file='crit.dat') 
  
    tbin=1.0d0/real(2*nbin+1)
       


    betac=0.0d0
    betax=0.0d0
    cxtp=0.0d0
    bdy=0.0d0
    do i=1,kb
       bxD=30.0d0+real(i)*0.03
       bdy(i)=bxD
       
       res=0.0d0
       call integrate_simpson(hi, 0.0d0, tpi, 10000, res, bxD)
       gxpi=( exp(tpi*offx)-1.0d0 )/res
       
       zpcx(i)=0.0d0
       cxtp(i)=0.0d0
       Pcx=0.0d0                   
       do k=-Ndx+1,Ndx-1
          xip=real(k)*rdx
          if(i<0) xip=xip+tpi
          gxx=0.0d0
          call integrate_simpson(hi, 0.0d0, xip,10000, gxx, bxD)    
          Pcx(k)=exp(-offx*xip+bxD*cos(xip) )*(1.0d0+gxx*gxpi)   
          zpcx(i)=zpcx(i)+Pcx(k)
       enddo
       do k=-Ndx+1,Ndx-1
          Pcx(k)=Pcx(k)/zpcx(i)
          xip=real(k)*rdx
          cxtp(i)=cxtp(i)+Pcx(k)*cos(xip)
       enddo
    enddo

    do i=1,kb
       do j=1,kbeta
           betax(j)=1.7d0+real(j)*0.003    
           beta0=betax(j)*cxtp(i)-2.0d0
           if(abs(beta0)<0.005) then
              write(70,*)  betax(j),bdy(i),betac(i,j)
           endif
       enddo
    enddo



    
    stop

    
contains
  real(kind=8) function hi(x, bxD)
     use param
     real(kind=8),intent(in)::x, bxD
        hi= exp( (offx*x-bxD*cos(x)) )
     return
end function hi

    
  end program main

!************************************************************************

subroutine integrate_simpson(hi, a, b, n, res, bxD)

  
    implicit none
    interface
        real(kind=8) function hi(x, bxD)
        real(kind=8),intent(in) :: x, bxD
        end function hi
    end interface

    real(kind=8),intent(in) :: a, b, bxD
    integer,intent(in) :: n
    real(kind=8),intent(out) :: res

    real(kind=8) :: h, x
    integer :: i

    if (mod(n,2) /= 0) then
        print *, "n must be even for Simpson rule"
        stop
    end if

    
    h = (b - a) / (real(n))
    res = hi(a, bxD) + hi(b, bxD)

    do i = 1, n-1
        x = a + h*real(i)
        if (mod(i,2) == 0) then
            res = res + 2.0d0*hi(x, bxD)
        else
            res = res + 4.0d0*hi(x, bxD)
        end if
    end do

    res = res*h/3.0d0
return
end subroutine integrate_simpson


