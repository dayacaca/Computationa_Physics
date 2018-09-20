 
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
  use chebys
  use arreglos


call allocate
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
call save_solution

! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! ======================================================== 

  end program Poisson_2D

! ======================================================== 
subroutine save_solution
use arreglos
use chebys

implicit none
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
end subroutine save_solution
