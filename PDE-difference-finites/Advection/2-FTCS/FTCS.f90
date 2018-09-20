
program FTCS

  implicit none 

  integer i,l 
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

  real(kind=8), allocatable, dimension (:) :: x
  real(kind=8), allocatable, dimension (:) :: u_e
  real(kind=8), allocatable, dimension (:) :: u_n
  real(kind=8), allocatable, dimension (:) :: u_n_p


!!$ Asking for parameters 

  print *, 'The left end of the spatial domain'
  read(*,*) xmin
  print *

  print *, 'The right end of the spatial domain'
  read(*,*) xmax
  print *

  print *, 'Number of points in the spatial domain'
  read(*,*) Nx
  print *

  print *, 'Number of time iterations'
  read(*,*) Nt
  print *

  print *, 'Aplitude of the gaussian'
  read(*,*) amp
  print *

  print *, 'Width of the initial gaussian'
  read(*,*) sigma
  print *

  print *, 'Center of the gaussian'
  read(*,*) x0
  print *

  print *, 'Courant Factor'
  read(*,*) courant
  print *

!!$ Allocate memory to the arrays

  allocate(x(0:Nx))
  allocate(u_e(0:Nx))
  allocate(u_n(0:Nx))
  allocate(u_n_p(0:Nx))  
   
!!$ Setting the spatial grid

  dx = ( xmax - xmin ) / dble( Nx ) !!$ Spatial resolution

  dt = courant * dx

  do i=0,Nx

     x(i)  = xmin + dble(i) * dx

  end do

!!$ Initial data

  u_e = amp * exp( - ( x - x0 )**2 / sigma**2 ) + 1e-20

  u_n = amp * exp( - ( x - x0 )**2 / sigma**2 ) + 1e-20
  
  
!!$ Saving the initial stuff
  t = 0.0d0
  call save1Ddata(Nx,t,x,u_e,'u_e',0)
  call save1Ddata(Nx,t,x,u_n,'u_n',0)


!!$ Integration loop

  do l=1,Nt
   
    t = t + dt
    u_n_p=u_n
  

  do i=1,Nx-1
    u_n(i)=u_n_p(i)- 0.5*(dt/dx)*(u_n_p(i+1)-u_n_p(i-1))
  end do
  
     u_e = amp * exp( - ( x - x0 - t)**2 / sigma**2 ) + 1e-20

!Coloco condiciones periodicas a las funciones u_e y u_n
   
	u_e(Nx)=u_e(1)
    	u_e(0)=u_e(Nx-1)
	u_n(Nx)=u_n(1)
	u_n(0)=u_n(Nx-1)
	
     call save1Ddata(Nx,t,x,u_e,'u_e',1)
     call save1Ddata(Nx,t,x,u_n,'u_n',1)
	
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
