module chebys
use arreglos
  implicit none
 
  contains

  Real (kind=8) FUNCTION cheby(m,r)
  implicit none

  integer  j,m
  real(kind=8)  r
  real(kind=8)  cheby0,cheby1,chebyn
  
  cheby0 = 1.0d0
  cheby1 = r 

  if(m.eq.0) then
    cheby = cheby0
  else if(m.eq.1) then 
    cheby = cheby1 
  else if(m.ge.2) then
    do j=2,m
       chebyn = 2.0d0*r*cheby1 - cheby0
       cheby0 = cheby1
       cheby1 = chebyn
    end do
       cheby = chebyn
  end if
 
 end FUNCTION cheby

! --------------------------------------
! defining firts Chebychev derivatives
 
  Real (kind=8) FUNCTION dcheby(m,r)

  implicit none

  integer j,m

  real(kind=8) r
  real(kind=8) dcheby0,dcheby1,dchebyn


  dcheby0 = 0.0d0
  dcheby1 = 1.0d0 

  if(m.eq.0) then
    dcheby = dcheby0
  else if(m.eq.1) then 
    dcheby = dcheby1 
  else if(m.ge.2) then
    do j=2,m
       dchebyn = 2.0*cheby(j-1,r) + 2.0*r*dcheby1 - dcheby0
       dcheby0 = dcheby1
       dcheby1 = dchebyn
    end do
       dcheby = dchebyn
  end if

 end FUNCTION dcheby

! --------------------------------------
! defining second Chebychev derivatives
 
  Real (kind=8) FUNCTION ddcheby(m,r)

  implicit none

  integer j,m

  real(kind=8) r
  real(kind=8) ddcheby0,ddcheby1,ddchebyn 

  ddcheby0 = 0.0d0
  ddcheby1 = 0.0d0 

  if(m.eq.0) then
    ddcheby = ddcheby0
  else if(m.eq.1) then 
    ddcheby = ddcheby1 
  else if(m.ge.2) then
    do j=2,m
       ddchebyn = 4.0*dcheby(j-1,r)+2.0*r*ddcheby1 - ddcheby0;
       ddcheby0 = ddcheby1
       ddcheby1 = ddchebyn
    end do
       ddcheby = ddchebyn
  end if

 end FUNCTION ddcheby

! ======================

  Real(kind=8) FUNCTION exact(p,q)

  implicit none
 
  REAL(kind=8) p,q,a0,b0,pii
 
  pii = -4.0d0*ATAN(1.0d0)
  a0 = 4.0d0
  b0 = 4.0d0

  exact = sin(a0*pii*p)*sin(b0*pii*q)
 
  end FUNCTION exact

! ======================

  Real(kind=8) FUNCTION f(p,q)

  implicit none
 
  ! p->x q->y 

  REAL(kind=8) p,q,pii,a0,b0
  
  pii = -4.0d0*ATAN(1.0d0)
  a0 = 4.0d0
  b0 = 4.0d0

  f = -(a0**2 + b0**2)*pii**2*sin(a0*pii*p)*sin(b0*pii*q)
!   f = exp(-(p**2+q**2)/0.1**2)

  end FUNCTION f

end module chebys
