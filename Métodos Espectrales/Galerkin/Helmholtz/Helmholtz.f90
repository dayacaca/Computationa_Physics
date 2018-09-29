 
!**********************************************************************!         
!   F.D.Lora-Clavijo and A. Cruz-Osorio - IFM - UMSNH                  ! 
!                                                                      !  
!   THIS PROGRAM SOLVE Helmholtz EQUATION  -y''(x)+lamda*y(x)=f(x)     ! 
!   WITH BOUNDARY CONDITIONS y(1)=y(-1)=0 USING THE                    ! 
!   SPECTRAL COLOCATION METHOD with THE CHEBYSHEV POLYNOMIALS          ! 
!                                                           !
!   This equation describes mass transfer processes                    !
!   with volume chemical reactions of the first order, lambda > 0                  !
!                                                                     
!**********************************************************************!
!=======================================================================
! Edited program 												                                    =
! Purpose: Program modudalize for Darly Y. Castro to learn Galerkin Method  =

!=======================================================================


  program Helmholtz

  use lu
  !use elements

  implicit none 

  integer i,l,m,j
  integer N,d,code,Nx

  REAL(kind=8) pii
  REAL(kind=8) suma
  REAL(kind=8) exact
  REAL(kind=8) lambda
  REAL(kind=8) dx
  REAL(kind=8) nm2, nm1
  real(kind=8) phi, dphi                                
  REAL(kind=8) xmax 
  REAL(kind=8) xmin   
  REAL(kind=8) integrate
  REAL(kind=8) intf
 ! REAL(kind=8) Phie
 ! REAL(kind=8) dPhie

  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: A  ! matriz
  !REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: x    ! puntos de colocación
  REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: b    ! vector fuente 
  REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: aa   ! vector solucion, coeficientes de la expansión
  integer, ALLOCATABLE, DIMENSION(:) :: indx
  REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: XX
  REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: error

  pii = 4.0d0*ATAN(1.0d0)

  print *, 'order of Phi'
  read(*,*) N
  print *
  
  print *, 'Number of grid points'
  read(*,*) Nx
  print *

  print *, 'lambda = 4.0 default'
  read(*,*) lambda
  print * 

  ALLOCATE(A(1:N+1,1:N+1))
 ! ALLOCATE(x(0:N))
  ALLOCATE(b(1:N+1))
  ALLOCATE(indx(1:N+1))
  ALLOCATE(aa(0:N))
  ALLOCATE(XX(0:Nx))
  ALLOCATE(error(0:Nx))

  xmax = 1.0d0
  xmin = -1.0d0  

 
! :::::::::::::::::::::::::::::::::::::::::::::::

! Matrix de coeficientes-- 
! i -> corresponde al orden del chebychev
! l -> corresponde al orden del chebychev
 
 
 do i=0,N
     do l=0,N
      A(i,l) = integrate(j,i)
     end do
  end do 

! :::::::::::::::::::::::::::::::::::::::::::::::::::::

! Aqui vamos a definir el vector b -> fuentes de la ecuación

  do i=0,N
   !     intf =(( xmax-  xmin)/6)*(f(i,xmin) +4*f(i,0.5*(xmax+xmin)) +f(i,xmax))
        b(i) = intf(i)
  end do              
 
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::
 
  !print*, '=========================================='

 
  !print*,'bvector',b(:)
  !print*, '=========================================='
 ! print*,'A1l',A(1,:)
 ! print*,'A2l',A(2,:) 
 ! print*,'A3l',A(3,:)
  !print*,'A4l',A(4,:)
 ! print*,'A5l',A(5,:) 
  !print*, '=========================================='
  
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::

! Descomposición LU
  
  call LUDCMP(A,N+1,INDX,D,CODE)
  call LUBKSB(A,N+1,INDX,B)

  !print*, 'b',b(:)
  
! Después de la inversión de la matriz utilizando la descomposición LU
! el vector b que era de fuentes ahora es el vector solución

  do l=0,N
     aa(l)=b(l)
  end do

 ! print*, 'aa', aa(:)   

! Aquí para no confundir asociamos el vector b solución a un vector aa(i) que son 
! los coeficientes de la expasión de la solución y = suma an Tn(x)

  
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 
! Salvando la solución

         dx = (xmax - xmin)/dble(Nx)

         do i=0,Nx
            XX(i) = xmin + dble(i)*dx
         end do

     open(1,file='Solution.dat') 
!      write(1,*) 'x      ','suma      ','exact       '
     do i=0,Nx 
       suma = 0.0d0 
      ! Phie= 0.0d0 
      ! dPhie= 0.0d0 
     do l=0,N
       suma = suma + aa(l)*phi(l,XX(i))
      ! Phie= Phie+phi(l,XX(i))
      ! dPhie= dPhie + dphi(l,XX(i))
     end do
       error(i) = abs(exact(XX(i)) - suma)
       write(1,*) XX(i),suma,exact(XX(i)),error(i)  !,Phie,dPhie
     end do
     close(1)

! :::::::::::::::::::::::::::::::::::::

     nm1 = 0.0d0
     nm2 = 0.0d0
     do i=1,Nx
      nm1 = nm1 + 0.5D0*(abs(error(i-1)) + abs(error(i)))*dx
      nm2 = nm2 + 0.5D0*(error(i-1)**2 + error(i)**2)*dx
    end do 

    nm2 = sqrt(nm2)

    print*,'nm2=',nm2  
    print*,'nm1=',nm1  

! ::::::::::::::::::::::::::::::::::::::

  end program Helmholtz

  
!===================================================================
!Def. de la nueva base de pol. de tal forma que phi(-1)=phi(1)=0  
!====================================================================

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
 real(kind=8) cheby  

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
  real(kind=8) cheby,dcheby  

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

  Real(kind=8) FUNCTION exact(r)
  implicit none
 
  REAL(kind=8) r

  exact = (-3.0d0/(8.0d0*(exp(2.0d0) + exp(-2.0d0))))*(exp(2.0d0*r) &
  + exp(-2.0d0*r)) +(1.0d0/4.0d0)*r*r + (1.0d0/8.0d0)  
  end FUNCTION exact

! ======================

  Real(kind=8) FUNCTION phi(m,r)
  implicit none 
  integer m
  real(kind=8) r
  real(kind=8) cheby
  
  if(mod(m,2).eq.0) then
    phi = cheby(m+2,r)-cheby(0,r)
  else if (mod(m,2).eq.1) then 
    phi = cheby(m+2,r) - cheby(1,r)
  end if
  
 end FUNCTION phi

 Real(kind=8) FUNCTION dphi(m,r)
  implicit none 
    integer m
    real(kind=8) r
    real(kind=8) dcheby  
  
  if(mod(m,2).eq.0) then
    dphi = dcheby(m+2,r) - dcheby(0,r)
  else if (mod(m,2).eq.1) then 
    dphi = dcheby(m+2,r) - dcheby(1,r)
  end if 
  
 end FUNCTION dphi


     
!===================================================================


 Real(kind=8) FUNCTION g(j,i,r)
    implicit none 
     integer i,j
     real(kind=8) r
     real(kind=8) dphi, phi
     real(kind=8) lambda
  
         g = -dphi(j,r)*dphi(i,r) + lambda*phi(j,r)*phi(i,r)
  
 end FUNCTION g 


 Real(kind=8) FUNCTION integrate(j,i)
 implicit none 
     integer i, j
    real(kind=8) g
     Real(kind=8) xmax, xmin
     xmax = 1.0d0
     xmin = -1.0d0 
   
     integrate =(( xmax -  xmin)/6)*(g(j,i,xmin) +4*g(j,i,0.0d0) +g(j,i,xmax)) 

 end FUNCTION integrate



 Real(kind=8) FUNCTION f(r)
    implicit none
    real(kind=8) r
  
     f = r*r
 end FUNCTION f


 Real(kind=8) FUNCTION intf(m)
     implicit none
     integer m
 REAL(kind=8) f,phi
 Real(kind=8) xmax, xmin
 xmax = 1.0d0
 xmin = -1.0d0 

 intf =((xmax-xmin)/6)*( f(xmin)*phi(m,xmin) + +4*(f(0.0d0)*phi(m,0.0d0)) &
  + f(xmax)*phi(m,xmax) )

end FUNCTION intf


