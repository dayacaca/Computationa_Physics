 !**********************************************************************!         
!   F.D.Lora-Clavijo and A. Cruz-Osorio - IFM - UMSNH                  ! 
!                                                                      !  
!   THIS PROGRAM SOLVE Helmholtz EQUATION  -y''(x)+lamda*y(x)=f(x)     ! 
!   WITH BOUNDARY CONDITIONS y(1)=y(-1)=0 USING THE                    ! 
!   SPECTRAL COLOCATION METHOD with THE CHEBYSHEV POLYNOMIALS          ! 
!                                                                      !
!   This equation describes mass transfer processes                    !
!   with volume chemical reactions of the first order, lambda > 0                  !
!                                                                      !
!**********************************************************************!
!=======================================================================
! Edited program 														|
! Purpose: Program modudalize for Darly Y. Castro to learn to use		|
! subroutines and modules.												|
! Structure: 															|
!***** Modules:  														|
! 			Chebys: contiene las funciones que defines los polinónmios  |
!					 de Chebychev con sus 1era y 2das derivadas 		|
!					 también define la función exacta y la homogenea	|
!			lu:		descomposición LU  invertir la matriz				|
!			arreglos : Declara todas las varibles y arreglos globales	|
!***** Subroutines:														|
!			  allocate: Da memoria a los arreglos del modulo arreglo	|
!  			  save_solution: Guarda la sln númerica y def error y norma	|
! 																		|
!=======================================================================

 
  program main
    use LU
    use chebys
	use arreglos
  implicit none 

call allocate

! Puntos de colocación

  xmax = 1.0d0
  xmin = -1.0d0  

  do i=0,N
    x(i) = cos(pii*i/N)
  end do


! Matrix de coeficientes
! i -> corresponde a los puntos de colocación
! l -> corresponde al orden del chebychev

  do i=0,N
     do l=0,N
        if (i.eq.0.or.i.eq.N) then
           A(i+1,l+1) = phi(l,x(i))
        else 
           A(i+1,l+1) = ddphi(l,x(i))-4.0d0*dphi(l,x(i)) + 4.0d0*phi(l,x(i))
        end if 
     end do
  end do 

! :::::::::::::::::::::::::::::::::::::::::::::::::::::

! Aqui vamos a definir el vector b -> fuentes de la ecuación

  do i=0,N
     if (i.eq.0.or.i.eq.N) then
       b(i+1) = 0.0d0
     else 
       b(i+1) = f(x(i))
     end if
  end do
 

! Descomposición LU
  
  call LUDCMP(A,N+1,INDX,D,CODE)
  call LUBKSB(A,N+1,INDX,B)
  
! Después de la inversión de la matriz utilizando la descomposición LU
! el vector b que era de fuentes ahora es el vector solución

  do l=0,N
     aa(l)=b(l+1)
  end do
 call save_solution
  
 end program main


!!!!! Guardando la solución 
 subroutine save_solution 
  use arreglos
  use chebys
 implicit none
         dx = (xmax - xmin)/dble(Nx)
         do i=0,Nx
            XX(i) = xmin + dble(i)*dx
         end do
         
     open(1,file='Solution.dat') 
     
     do i=0,Nx 
       suma = 0.0d0 
     do l=0,N
       suma = suma + aa(l)*cheby(l,XX(i)) 
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

end subroutine save_solution 

 
