 
!**********************************************************************!         
!   F.D.Lora-Clavijo and A. Cruz-Osorio - IFM - UMSNH                  ! 
!                                                                      !  
!   Poisson 2D cartesian                                               ! 
!   First Version. 16/07/12                                            !
!   Colombo-Mexican Elliptic-Solver using Spectral Methods            
!                                                                      !
!**********************************************************************!

  program Poisson_2D

  use LU

  implicit none 

  integer i,j,k,l
  integer N,d,code

  integer Nx,Ny

  REAL(kind=8) pii
  REAL(kind=8) cheby
  REAL(kind=8) dcheby
  REAL(kind=8) ddcheby
  REAL(kind=8) suma
  REAL(kind=8) exact
  REAL(kind=8) lambda
  REAL(kind=8) f

  REAL(kind=8) dx,dy,nm2
  REAL(kind=8) xmin,xmax,ymin,ymax
 
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: A  ! matriz
  REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: x    ! puntos de colocación x
  REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: y    ! puntos de colocación y
  REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: b    ! vector fuente 
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: aa ! coeficientes de la expansión
  INTEGER, ALLOCATABLE, DIMENSION(:) :: indx
  REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: XX ! Malla 2D donde se va a evaluar la función
  REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: YY
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: error

  pii = 4.0d0*ATAN(1.0d0)
  xmin = -1.0d0
  xmax = 1.0d0
  ymin = -1.0d0
  ymax = 1.0d0 

  print *, 'Number of colocation points'
  read(*,*) N
  print *

  print *, 'Number of points on the x direction'
  read(*,*) Nx
  print *

  print *, 'Number of points on the y direction'
  read(*,*) Ny
  print *

!  print *, 'lambda = 4.0 default'
!  read(*,*) lambda
!  print * 

  ALLOCATE(A(1:(N+1)*(N+1),1:(N+1)*(N+1)))
  ALLOCATE(x(0:N))
  ALLOCATE(y(0:N))
  ALLOCATE(b(1:(N+1)*(N+1)))
  ALLOCATE(indx(1:(N+1)*(N+1)))
  ALLOCATE(aa(0:N,0:N))
  ALLOCATE(XX(0:Nx))
  ALLOCATE(YY(0:Ny))
  ALLOCATE(error(0:Nx,0:Ny))

! :::::::::::::::::::::::::::::::::::::::::::::::
! Puntos de colocación en cada dirección

  do i=0,N
    x(i) = cos(pii*i/N)
    y(i) = cos(pii*i/N)
  end do

! :::::::::::::::::::::::::::::::::::::::::::::::

! Matrix de coeficientes
! i -> corresponde a los puntos de colocación x
! j -> corresponde a los puntos de colocación y
! l,k -> corresponden al orden del chebychev

  do i=0,N
     do j=0,N
        do l=0,N
           do k=0,N 
               if (i.eq.0.or.i.eq.N.or.j.eq.0.or.j.eq.N) then
                  A(i*(N+1)+(j+1),l*(N+1)+(k+1)) = cheby(l,x(i))*cheby(k,y(j)) !Boundary conditions
               else 
                  A(i*(N+1)+(j+1),l*(N+1)+(k+1)) = ddcheby(l,x(i))*cheby(k,y(j)) &
                                                 + cheby(l,x(i))*ddcheby(k,y(j))
               end if 
           end do
        end do 
     end do
  end do 

!  print*,'A4:',A(4,:)

! :::::::::::::::::::::::::::::::::::::::::::::::::::::

! Aqui vamos a definir el vector b -> fuentes de la ecuación

   do i=0,N
     do j=0,N
        if (i.eq.0.or.i.eq.N.or.j.eq.0.or.j.eq.N) then
          b(i*(N+1)+(j+1)) = 0.0d0  !Boundary conditions
        else 
          b(i*(N+1)+(j+1)) = f(x(i),y(j)) 
        end if
     end do
   end do 
 
!  print*,'bvector',b(:)

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
  
  call LUDCMP(A,(N+1)*(N+1),INDX,D,CODE)
  call LUBKSB(A,(N+1)*(N+1),INDX,B)

!  print*, 'b',b(:)
  
! Después de la inversión de la matriz utilizando la descomposición LU
! el vector b que era de fuentes ahora es el vector solución

  do l=0,N
    do k=0,N
      aa(l,k)=b(l*(N+1)+(k+1))
    end do
  end do

!  print*, 'aa', aa(:,:)   

! Aquí para no confundir asociamos el vector b solución a las componentes anm(i,j) que son 
! los coeficientes de la expasión de la solución y = suma anm Tn(x) Tm(y)


! ========================================================  
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! ======================================================== 

! Salvando la solución

! malla donde vamos a salvar datos

         dx = (xmax - xmin)/dble(Nx)
         dy = (ymax - ymin)/dble(Ny)

          do i=0,Nx
             XX(i)  = xmin + dble(i)*dx
          end do

          do j=0,Ny
             YY(j)  = ymin + dble(j)*dy
          end do

     open(1,file='Solution.dat') 
!     write(1,*) 
     do i=0,Nx
     write(1,*) 
     do j=0,Ny 
     suma = 0.0d0 
     do l=0,N
     do k=0,N
         suma = suma + aa(l,k)*cheby(l,XX(i))*cheby(k,YY(j)) 
         error(i,j) = exact(XX(i),YY(j)) - suma 
     end do 
!         write(1,*) XX(i),YY(i),suma
     end do
         write(1,*) XX(i),YY(j),suma,exact(XX(i),YY(j)),error(i,j)
     end do
     end do  
     close(1) 

! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

     nm2 = 0.0d0
     do i=0,Nx-1
     do j=0,Ny-1
      nm2 = nm2 + 0.25D0*( (error(i,j)**2 + error(i+1,j)**2) &
                         + (error(i,j+1)**2 + error(i+1,j+1)**2) )*dx*dy
    end do
    end do

    nm2 = sqrt(nm2)

    print*,N,nm2

! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

! ========================================================  
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! ======================================================== 

  end program Poisson_2D

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


! ===================================================================
