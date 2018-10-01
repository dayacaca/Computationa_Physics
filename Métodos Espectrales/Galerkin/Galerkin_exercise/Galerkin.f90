 
!**********************************************************************!         
!   F.D.Lora-Clavijo and A. Cruz-Osorio - IFM - UMSNH                  ! 
!                                                                      !  
!   THIS PROGRAM SOLVES THE EQUATION  y''(x)-4*y'(x)+4*y(x)=f(x)       ! 
!   WITH BOUNDARY CONDITIONS y(1)=y(-1)=0 USING THE                    ! 
!   SPECTRAL GALERKIN METHOD with THE CHEBYSHEV POLYNOMIALS            ! 
!                                                                      !
!                                                                      !
!**********************************************************************!

  program Colocation

  use LU

  implicit none 

  integer i,l
  integer N,d,code,Nx

  REAL(kind=8) pii
  REAL(kind=8) cheby
  REAL(kind=8) dcheby
  REAL(kind=8) ddcheby
  REAL(kind=8) phi
  REAL(kind=8) dphi
  REAL(kind=8) lhs, rhs
  REAL(kind=8) suma
  REAL(kind=8) exact
  !REAL(kind=8) f
  REAL(kind=8) xmin, xmax, dx
  REAL(kind=8) nm2, nm1
  REAL(kind=8) integral
  double precision s
  double precision h,z
  integer ni,k

 
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: A  ! matriz
  !REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: x    ! puntos de colocación
  REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: b    ! vector fuente 
  REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: aa   ! vector solucion, coeficientes de la expansión
  integer, ALLOCATABLE, DIMENSION(:) :: indx
  REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: XX
  REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: error

  pii = 4.0d0*ATAN(1.0d0)

  print *, 'Order of the polynomial (N = 4 best fit)'
  read(*,*) N
  print *

  print *, 'Number of subintervals (for the integration)'
  read(*,*) ni
  print *

  
  print *, 'Number of grid points'
  read(*,*) Nx
  print *

  !print *, 'lambda = 4.0 default'
  !read(*,*) lambda
  !print * 

  ALLOCATE(A(1:N+1,1:N+1))
  !ALLOCATE(x(0:N))
  ALLOCATE(b(1:N+1))
  ALLOCATE(indx(1:N+1))
  ALLOCATE(aa(0:N))
  ALLOCATE(XX(0:Nx))
  ALLOCATE(error(0:Nx))
  

! :::::::::::::::::::::::::::::::::::::::::::::::
! Puntos de colocación


 xmax = 1.0d0
 xmin = -1.0d0  

  !do i=0,N
  !  x(i) = cos(pii*i/N)
  !end do

! :::::::::::::::::::::::::::::::::::::::::::::::

! Matrix de coeficientes
! i -> corresponde a los puntos de colocación
! l -> corresponde al orden del chebychev

  do i=0,N
     do l=0,N

        !Integration using a composite Simpson's rule
          integral=0
          if((ni/2)*2.ne.ni) ni=ni+1

          ! loop over ni (number of intervals)
          s = 0.0
          h = (xmax-xmin)/dfloat(ni)
          do k=2, ni-2, 2
            z   = xmin+dfloat(k)*h
            s = s + 2.0d0*lhs(i,l,z) + 4.0d0*lhs(i,l,z+h)
          end do
          integral = (s + lhs(i,l,xmin) + lhs(i,l,xmax) + 4.0d0*lhs(i,l,xmin+h))*h/3.0d0

          A(l+1,i+1)=integral
          !A(i+1,l+1)= (xmax-xmin)/(6.0d0)*(lhs(i,l,xmin)+4.0d0*lhs(i,l,(xmin+xmax)/2.0d0) + lhs(i,l,xmax))
     end do
  end do 

! :::::::::::::::::::::::::::::::::::::::::::::::::::::

! Aqui vamos a definir el vector b -> fuentes de la ecuación

  do i=0,N
    ! integrating
    integral=0
    if((ni/2)*2.ne.ni) ni=ni+1

    ! loop over n (number of intervals)
    s = 0.0
    h = (xmax-xmin)/dfloat(ni)
    do k=2, ni-2, 2
      z   = xmin+dfloat(k)*h
      s = s + 2.0d0*rhs(i,z) + 4.0d0*rhs(i,z+h)
    end do
    integral = (s + rhs(i,xmin) + rhs(i,xmax) + 4.0d0*rhs(i,xmin+h))*h/3.0d0
      b(i+1) = integral
  end do
 
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::
 
!  print*, '=========================================='
!  print*,'bvector',b(:)
!  print*, '=========================================='
!  print*,'A1l',A(1,:)
!  print*,'A2l',A(2,:) 
!  print*,'A3l',A(3,:)
!  print*,'A4l',A(4,:)
!  print*,'A5l',A(5,:) 
!  print*, '=========================================='
  
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::

! Descomposición LU
  
  call LUDCMP(A,N+1,INDX,D,CODE)
  call LUBKSB(A,N+1,INDX,B)

!  print*, 'b',b(:)
  
! Después de la inversión de la matriz utilizando la descomposición LU
! el vector b que era de fuentes ahora es el vector solución

  do l=0,N
     aa(l)=b(l+1)
  end do

!  print*, 'aa', aa(:)   

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
     do l=0,N
       !suma = suma + aa(l)*cheby(l,XX(i)) 
      suma = suma + aa(l)*phi(l,XX(i))
     end do
       error(i) = abs(exact(XX(i)) - suma)
       write(1,*) XX(i),suma,exact(XX(i)),error(i)
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

  end program Colocation

! ----------------------------------------
! defining Chebychev Polynomials 

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
! defining phi functions

  Real (kind=8) FUNCTION phi(m,r)
  
  implicit none

  integer m

  real(kind=8) r
  real(kind=8) cheby

  if(mod(m,2).eq.0) then
    phi=cheby(m+2,r)-cheby(0,r)
  else if(mod(m,2).eq.1) then
    phi=cheby(m+2,r)-cheby(1,r)
  end if
  
  end FUNCTION phi 


! ======================
! defining phi functions first derivatives

  Real (kind=8) FUNCTION dphi(m,r)
  
  implicit none

  integer m

  real(kind=8) r
  real(kind=8) dcheby

  if(mod(m,2).eq.0) then
    dphi=dcheby(m+2,r)-dcheby(0,r)
  else if(mod(m,2).eq.1) then
    dphi=dcheby(m+2,r)-dcheby(1,r)
  end if
  
  end FUNCTION dphi


! ======================
! defining phi functions second derivatives

  Real (kind=8) FUNCTION ddphi(m,r)
  
  implicit none

  integer m

  real(kind=8) r
  real(kind=8) ddcheby

  if(mod(m,2).eq.0) then
    ddphi=ddcheby(m+2,r)-ddcheby(0,r)
  else if(mod(m,2).eq.1) then
    ddphi=ddcheby(m+2,r)-ddcheby(1,r)
  end if
  
  end FUNCTION ddphi


 !========================================
  Real(kind=8) FUNCTION exact(r)

  implicit none
 
  REAL(kind=8) r

  exact =  exp(r)-((sinh(1.0d0))/(sinh(2.0d0)))*exp(2.0d0*r)-((4.0d0*exp(1.0d0))/(1+exp(2.0d0)))/4 
 
  end FUNCTION exact

! ======================
  ! Defining the new RHS to be integrated

  Real(kind=8) FUNCTION rhs(m,r)

  implicit none

  integer m

  REAL(kind=8) r
  REAL(kind=8) phi

  rhs = (exp(r)-(4.0d0*exp(1.0d0))/(1+exp(2.0d0)))*phi(m,r)  
 
  end FUNCTION rhs


! ======================
! Defining the function in the LHS to be integrated

  Real(kind=8) FUNCTION lhs(m,n,r)
  
  implicit none

  integer m
  integer n
  
  REAL(kind=8) r
  REAL(kind=8) phi, dphi

  lhs = -dphi(m,r)*dphi(n,r)-4*dphi(m,r)*phi(n,r)+4*phi(m,r)*phi(n,r)

  end FUNCTION lhs
! ===================================================================
! Defining the function to carry on the integration
!
 ! Real(kind=8) FUNCTION integral(g,a,b)

  !implicit none
  
 ! real(kind=8) g
 ! real(kind=8) a
 ! real(kind=8) b
  




!  end FUNCTION integral


