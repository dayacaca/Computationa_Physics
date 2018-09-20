 
!**********************************************************************!         
!   F.D.Lora-Clavijo and A. Cruz-Osorio - IFM - UMSNH                  ! 
!                                                                      !  
!   Poisson 3D cartesian                                               ! 
!   First Version. 17/07/12                                            !
!   Colombo-Mexican Elliptic-Solver using Spectral Methods (CoMeSESo)  !                                            
!                                                                      !
!**********************************************************************!

  program Poisson_3D

  use LU

  implicit none 

  integer i,j,k,l,m,s
  integer N,d,code

  integer Nx,Ny,Nz

  REAL(kind=8) pii
  REAL(kind=8) cheby
  REAL(kind=8) dcheby
  REAL(kind=8) ddcheby
  REAL(kind=8) suma
  REAL(kind=8) exact
  REAL(kind=8) lambda
  REAL(kind=8) f

  REAL(kind=8) dx,dy,dz
  REAL(kind=8) xmin,xmax,ymin,ymax,zmin,zmax
 
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: A  ! matriz
  REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: x    ! puntos de colocación x
  REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: y    ! puntos de colocación y
  REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: z    ! puntos de colocación z 
  REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: b    ! vector fuente 
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:,:) :: aa ! coeficientes de la expansión
  INTEGER, ALLOCATABLE, DIMENSION(:) :: indx
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:,:) :: XX ! Malla 3D donde se va a evaluar la función
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:,:) :: YY
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:,:) :: ZZ

  pii = 4.0d0*ATAN(1.0d0)
  xmin = -1.0d0
  xmax = 1.0d0
  ymin = -1.0d0
  ymax = 1.0d0 
  zmin = -1.0d0
  zmax = 1.0d0 

  print *, 'Number of colocation points'
  read(*,*) N
  print *

  print *, 'Number of points on the x direction'
  read(*,*) Nx
  print *

  print *, 'Number of points on the y direction'
  read(*,*) Ny
  print *

  print *, 'Number of points on the z direction'
  read(*,*) Nz
  print *

!  print *, 'lambda = 4.0 default'
!  read(*,*) lambda
!  print * 

  ALLOCATE(A(1:(N+1)*(N+1)*(N+1),1:(N+1)*(N+1)*(N+1)))
  ALLOCATE(x(0:N))
  ALLOCATE(y(0:N))
  ALLOCATE(z(0:N))
  ALLOCATE(b(1:(N+1)*(N+1)*(N+1)))
  ALLOCATE(indx(1:(N+1)*(N+1)*(N+1)))
  ALLOCATE(aa(0:N,0:N,0:N))
  ALLOCATE(XX(0:Nx,0:Ny,0:Nz))
  ALLOCATE(YY(0:Nx,0:Ny,0:Nz))
  ALLOCATE(ZZ(0:Nx,0:Ny,0:Nz))

! :::::::::::::::::::::::::::::::::::::::::::::::
! Puntos de colocación en cada dirección

  do i=0,N
    x(i) = cos(pii*i/N)
    y(i) = cos(pii*i/N)
    z(i) = cos(pii*i/N)
  end do

! :::::::::::::::::::::::::::::::::::::::::::::::

! Matrix de coeficientes
! i -> corresponde a los puntos de colocación x
! j -> corresponde a los puntos de colocación y
! k -> corresponde a los puntos de colocación z
! l,m,s -> corresponden al orden del chebychev

  do i=0,N
     do j=0,N
        do k=0,N
           do l=0,N
              do m=0,N 
                 do s=0,N
                   if (i.eq.0.or.i.eq.N.or.j.eq.0.or.j.eq.N.or.k.eq.0.or.k.eq.N) then
                      A(i*(N+1)*(N+1)+j*(N+1)+(k+1),l*(N+1)*(N+1)+m*(N+1)+(s+1)) = &
                      + cheby(l,x(i))*cheby(m,y(j))*cheby(s,z(k)) 
                   else 
                      A(i*(N+1)*(N+1)+j*(N+1)+(k+1),l*(N+1)*(N+1)+m*(N+1)+(s+1)) = &
                      + ddcheby(l,x(i))*cheby(m,y(j))*cheby(s,z(k)) &
                      + cheby(l,x(i))*ddcheby(m,y(j))*cheby(s,z(k)) &
                      + cheby(l,x(i))*cheby(m,y(j))*ddcheby(s,z(k))                 
                   end if
                 end do
               end do 
           end do
        end do 
     end do
  end do 

!  print*,'A1:',A(1,:)

! :::::::::::::::::::::::::::::::::::::::::::::::::::::

! Aqui vamos a definir el vector b -> fuentes de la ecuación

   do i=0,N
     do j=0,N
       do k=0,N
          if (i.eq.0.or.i.eq.N.or.j.eq.0.or.j.eq.N.or.k.eq.0.or.k.eq.N) then
            b(i*(N+1)*(N+1)+j*(N+1)+(k+1)) = 0.0d0  !Boundary conditions
          else 
            b(i*(N+1)*(N+1)+j*(N+1)+(k+1)) = f(x(i),y(j),z(k)) 
          end if
       end do 
     end do
   end do 
 
!  print*,'bvector',b(:)

! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::

! Descomposición LU
  
  call LUDCMP(A,(N+1)*(N+1)*(N+1),INDX,D,CODE)
  call LUBKSB(A,(N+1)*(N+1)*(N+1),INDX,B)

!  print*, 'b',b(:)
  
! Después de la inversión de la matriz utilizando la descomposición LU
! el vector b que era de fuentes ahora es el vector solución

  do l=0,N
    do m=0,N
      do s=0,N
        aa(l,m,s)=b(l*(N+1)*(N+1)+m*(N+1)+(s+1))
      end do
    end do
  end do

!  print*, 'aa', aa(:,:,:)   

! Aquí para no confundir asociamos el vector b solución a las componentes anml(i,j,k) que son 
! los coeficientes de la expasión de la solución y = suma anml Tn(x) Tm(y) Tl(z)


! ========================================================  
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! ======================================================== 

! Salvando la solución

! malla donde vamos a salvar datos

         dx = (xmax - xmin)/dble(Nx)
         dy = (ymax - ymin)/dble(Ny)
         dz = (zmax - zmin)/dble(Nz)

         do i=0,Nx
             do j=0,Ny
                do k=0,Nz
                 XX(i,j,k)  = xmin + dble(i) * dx
                 YY(i,j,k)  = ymin + dble(j) * dy
                 ZZ(i,j,k)  = zmin + dble(k) * dz
                end do
             end do
          end do 
          

     open(1,file='Solution_xy.dat') 
     write(1,*) 
     do i=0,Nx
     write(1,*) 
     do j=0,Ny 
     suma = 0.0d0 
     do l=0,N
     do m=0,N
     do s=0,N
         suma = suma + aa(l,m,s)*cheby(l,XX(i,j,Nz/16))*cheby(m,YY(i,j,Nz/16))*cheby(s,ZZ(i,j,Nz/16)) 
     end do 
     end do
     end do
         write(1,*) XX(i,j,Nz/16),YY(i,j,Nz/16),suma,exact(XX(i,j,Nz/16),YY(i,j,Nz/16),ZZ(i,j,Nz/16))
     end do
     end do  
     close(1) 
 
     open(1,file='Solution_xz.dat') 
     write(1,*) 
     do i=0,Nx
     write(1,*) 
     do k=0,Nz 
     suma = 0.0d0 
     do l=0,N
     do m=0,N
     do s=0,N
         suma = suma + aa(l,m,s)*cheby(l,XX(i,Ny/16,k))*cheby(m,YY(i,Ny/16,k))*cheby(s,ZZ(i,Ny/16,k)) 
     end do 
     end do
     end do
         write(1,*) XX(i,Ny/16,k),ZZ(i,Ny/16,k),suma,exact(XX(i,Ny/16,k),YY(i,Ny/16,k),ZZ(i,Ny/16,k))
     end do
     end do  
     close(1) 
  
     open(1,file='Solution_yz.dat') 
     write(1,*) 
     do j=0,Ny
     write(1,*) 
     do k=0,Nz 
     suma = 0.0d0 
     do l=0,N
     do m=0,N
     do s=0,N
         suma = suma + aa(l,m,s)*cheby(l,XX(Nx/16,j,k))*cheby(m,YY(Nx/16,j,k))*cheby(s,ZZ(Nx/16,j,k)) 
     end do 
     end do
     end do
         write(1,*) YY(Nx/16,j,k),ZZ(Nx/16,j,k),suma,exact(XX(Nx/16,j,k),YY(Nx/16,j,k),ZZ(Nx/16,j,k))
     end do
     end do  
     close(1) 
   

! ========================================================  
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! ======================================================== 

  end program Poisson_3D

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

  Real(kind=8) FUNCTION exact(p,q,t)

  implicit none
 
  REAL(kind=8) p,q,t,a0,b0,c0,pii

  a0 = 4.0d0
  b0 = 4.0d0
  c0 = 4.0d0 
  pii = -4.0d0*ATAN(1.0d0)

  exact = sin(a0*pii*p)*sin(b0*pii*q)*sin(c0*pii*t)  
 
  end FUNCTION exact

! ======================

  Real(kind=8) FUNCTION f(p,q,t)

  implicit none
 
  ! p->x q->y t->z

  REAL(kind=8) p,q,pii,t
  REAL(kind=8) a0,b0,c0

  a0 = 4.0d0
  b0 = 4.0d0
  c0 = 4.0d0
  pii = -4.0d0*ATAN(1.0d0)

!  f = -exp(-(p**2+q**2+t**2)/0.1**2)  
 
  f = -(a0**2 + b0**2 + c0**2)*pii**2*sin(a0*pii*p)*sin(b0*pii*q)*sin(c0*pii*t) 

  end FUNCTION f


! ===================================================================
