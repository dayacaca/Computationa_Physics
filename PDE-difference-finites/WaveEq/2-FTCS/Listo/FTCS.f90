! Resuelve la ecuación de onda d_tt phi =v^2 d_xx phi
! la ec. necesita   **IVP**     **BP**
 
program FTCS

  implicit none 

  integer i,j,l 
  integer Nx
  integer Nt

!!$ Grid data

  real(kind=8) xmin 
  real(kind=8) xmax
  real(kind=8) dx
  real(kind=8) dt
  real(kind=8) t

!!$ Wave data

  real(kind=8) amp
  real(kind=8) x0
  real(kind=8) sigma
  real(kind=8) courant
!Declaro las funciones a calcular, en este caso son psi, phi y pi y tengo que usar 2 arreglso para cada función 

  real(kind=8), allocatable, dimension (:) :: x  !mallado 
  real(kind=8), allocatable, dimension (:) :: u_e ! función real
  
  real(kind=8), allocatable, dimension (:) :: pi_n 
  real(kind=8), allocatable, dimension (:) :: pi_n_p
  
  real(kind=8), allocatable, dimension (:) :: psi_n
  real(kind=8), allocatable, dimension (:) :: psi_n_p
  
  real(kind=8), allocatable, dimension (:) :: phi_n
  real(kind=8), allocatable, dimension (:) :: phi_n_p
  

!!$ Reading  parameters of file datos.dat
 open(10,file='datos.dat')
  read(10,*) xmin,xmax,Nx,Nt,amp,sigma,x0,courant
 close(10)
 

!!$ Allocate memory to the arrays

  allocate(x(0:Nx))
  allocate(u_e(0:Nx))
  allocate(pi_n(0:Nx))
  allocate(pi_n_p(0:Nx))  
  allocate(psi_n(0:Nx))  
  allocate(psi_n_p(0:Nx))  
  allocate(phi_n(0:Nx))  
  allocate(phi_n_p(0:Nx))  

!!$ Setting the spatial grid

  dx = ( xmax - xmin ) / dble( Nx ) !!$ Spatial resolution

  dt = courant * dx

  do i=0,Nx

     x(i)  = xmin + dble(i) * dx

  end do

!!$ Initial data

  u_e = amp * exp( - ( x - x0 )**2 / sigma**2 ) + 1e-20
  pi_n = 0d0
  psi_n =  -2*amp*((x-x0)/sigma**2) * exp( - ( x - x0 )**2 / sigma**2 ) 
  phi_n = amp * exp( - ( x - x0 )**2 / sigma**2 ) 
	
!!$ Saving the initial stuff

  t = 0.0d0

  call save1Ddata(Nx,t,x,u_e,'u_e',0)
  call save1Ddata(Nx,t,x,pi_n,'pi_n',0)
  call save1Ddata(Nx,t,x,psi_n,'psi_n',0)
  call save1Ddata(Nx,t,x,phi_n,'phi_n',0)
  
  
!!$ Integration loop

  do l=1,Nt
   t = t + dt
  psi_n_p = psi_n
  pi_n_p  = pi_n
  phi_n_p = phi_n
  

   do i=1,Nx-1
     pi_n(i) = pi_n_p(i)+ 0.5*(dt/dx)*(psi_n_p(i+1)-psi_n_p(i-1)) 
	 psi_n(i) = psi_n_p(i)+ 0.5*(dt/dx)*(pi_n_p(i+1)-pi_n_p(i-1))
   end do 
    !BV periódicas
    !OJO se coloca condiciones a pi y psi, ya que aún no he calculado phi 
  
    pi_n(Nx)=pi_n(1)
    pi_n(0)= pi_n(Nx-1) 	
     	
   psi_n(Nx)=psi_n(1)
   psi_n(0)= psi_n(Nx-1) 	 
   
 
     do j=0,Nx	
     phi_n(j)= phi_n_p(j)+ dt*pi_n(j)
     end do
    u_e(Nx)=u_e(1)
    u_e(0) = u_e(Nx-1) 
    
    
    u_e = 0.5*amp *( exp( - ( x - x0 - t)**2 / sigma**2 ) + exp( - ( x - x0 + t)**2 / sigma**2 ) ) + 1e-20   
    call save1Ddata(Nx,t,x,u_e,'u_e',1)
    call save1Ddata(Nx,t,x,phi_n,'phi_n',1)
	call save1Ddata(Nx,t,x,psi_n,'psi_n',1)
    call save1Ddata(Nx,t,x,pi_n,'pi_n',1)
    
	 
  end do

end program FTCS

!!$ 
!!$

subroutine save1Ddata(Nx_,t_,xval,yval,base_name,first_index)

  implicit none

  character(len=20) filestatus

  character(len=*), intent(IN) :: base_name
  real(kind=8), dimension(0:Nx_), intent(IN) :: xval, yval

  character(len=256) :: filename
  real(kind=8) t_

  integer i, Nx_, first_index

  filename = base_name // '.x'

  if (first_index.eq.0) then
     filestatus = 'replace'
  else
     filestatus = 'old'
  end if

  if (filestatus=='replace') then
     open(1,file=filename,form='formatted',status=filestatus)
  else
     open(1,file=filename,form='formatted',status=filestatus,position='append')
  end if
  write(1,*) ''
  write(1,*) '#"Time = ', t_
  do i=0,Nx_
     write(1,*) xval(i),yval(i)
  end do
  write(1,*)
  close(1)

end subroutine save1Ddata
