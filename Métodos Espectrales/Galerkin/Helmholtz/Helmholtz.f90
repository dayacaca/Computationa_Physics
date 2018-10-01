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
! Edited program 	for Darly Y. Castro					                                   
! Purpose: Galerkin Method   

!=======================================================================


  program Helmholtz
 
    use LU
    !use simpson
    !use simpsong
  
    implicit none 
  
    integer i,l,m,k
    integer N,d,code,Nx,ni
  
    REAL(kind=8) pii
    REAL(kind=8) suma
    REAL(kind=8) exact
    REAL(kind=8) lambda
    REAL(kind=8) dx
    REAL(kind=8) nm2, nm1
    real(kind=8) phi, dphi                                
    REAL(kind=8) xmax 
    REAL(kind=8) xmin   
    REAL(kind=8) g  
    REAL(kind=8) f
    REAL(kind=8) s,h,x
   !Phie
   ! REAL(kind=8) dPhie
    REAL(kind=8) integral
    REAL(kind=8) integralg
  
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
  
    print *, 'intervalos de las integrales'
    read(*,*) ni
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
        integralg=0.0d0

        if((ni/2)*2.ne.ni) ni=ni+1

        ! loop over n (number of intervals)
        s = 0.0
        h = (xmax-xmin)/dfloat(ni)
        do k=2, ni-2, 2
           x   = xmin+dfloat(k)*h
           s = s + 2.0d0*g(i,l,x) + 4.0d0*g(i,l,x+h)
        end do

        integralg = (s + g(i,l,xmin) + g(i,l,xmax) + 4.0*g(i,l,xmin+h))*h/3.0d0
        A(l+1,i+1) = integralg
       end do
    end do 
  
  ! :::::::::::::::::::::::::::::::::::::::::::::::::::::
  
  ! Aqui vamos a definir el vector b -> fuentes de la ecuación
    
    do m=0,N
      integral=0.0d0
      if((ni/2)*2.ne.ni) ni=ni+1
      s = 0.0
      h = (xmax-xmin)/dfloat(ni)
        do k=2, ni-2, 2
                x   = xmin+dfloat(k)*h
                s = s + 2.0d0*f(m,x) + 4.0d0*f(m,x+h)
         end do
      integral = (s + f(m,xmin) + f(m,xmax) + 4.0d0*f(m,xmin+h))*h/3.0d0
      b(m+1) = integral
    end do              
    

  ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::
   
    !print*, '=========================================='
  
   
   ! print*,'bvector',b(:)
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
       aa(l)=b(l+1)
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
  
    end program Helmholtz
  


     
  !===================================================================
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
    real(kind=8) dcheby  
  
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
  
  !===================================================================
  !Def. de la nueva base de pol. de tal forma que phi(-1)=phi(1)=0  
  !====================================================================
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


   Real(kind=8) FUNCTION g(i,j,r)
    implicit none 
     integer i,j,lambda
     real(kind=8) r
     real(kind=8) dphi, phi
         g= dphi(i,r)*dphi(j,r) + lambda*phi(i,r)*phi(j,r)
  
  end FUNCTION g 
  
  Real(kind=8) FUNCTION f(m,r)
    implicit none
    integer m
    real(kind=8) r,phi
  
     f = r*r*phi(m,r)
  end FUNCTION f
