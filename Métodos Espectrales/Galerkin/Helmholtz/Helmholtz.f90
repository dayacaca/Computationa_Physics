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
    use simpson
    !use simpsong
  
    implicit none 
  
    integer i,l,m
    integer N,d,code,Nx,ni
  
    REAL(kind=8) pii
    REAL(kind=8) suma
    !REAL(kind=8) exact
    REAL(kind=8) lambda
    REAL(kind=8) dx
    REAL(kind=8) nm2, nm1
    !real(kind=8) phi, dphi                                
    REAL(kind=8) xmax 
    REAL(kind=8) xmin   
    !REAL(kind=8) g  
    !REAL(kind=8) f
  
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
        call simpsong(g,i,l,xmin, xmax,integralg,ni)
        A(l+1,i+1) = integralg
       end do
    end do 
  
  ! :::::::::::::::::::::::::::::::::::::::::::::::::::::
  
  ! Aqui vamos a definir el vector b -> fuentes de la ecuación
    
    do m=0,N
      call simpsonf(f,m,xmin,xmax,integral,ni)
          b(m+1) = integral
    end do              
    

  ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::
   
    !print*, '=========================================='
  
   
   ! print*,'bvector',b(:)
    !print*, '=========================================='
    print*,'A1l',A(1,:)
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
        ! ge=0.0d0
        ! Phie= 0.0d0 
        ! dPhie= 0.0d0 
       do l=0,N
         suma = suma + aa(l)*phi(l,XX(i))
        ! Phie= Phie+phi(l,XX(i))
        ! dPhie= dPhie + dphi(l,XX(i))
       !  ge = ge+g(l,XX(i))
       end do
         error(i) = abs(exact(XX(i)) - suma)
         write(1,*) XX(i),suma,exact(XX(i)),error(i) !,Phie,dPhie
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
  
