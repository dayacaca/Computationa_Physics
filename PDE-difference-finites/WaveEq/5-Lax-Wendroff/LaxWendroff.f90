 

program Leadfrog

  implicit none 

  integer i,l
  integer Nx
  integer Nt

!!$ dotos de la malla

  real(kind=8) xmin 
  real(kind=8) xmax
  real(kind=8) dx
  real(kind=8) dt
  real(kind=8) t
  real(kind=8) error

!!$ Datos de la onda Gaussiana

  real(kind=8) amp
  real(kind=8) x0
  real(kind=8) sigma
  real(kind=8) courant

  real(kind=8), allocatable, dimension (:) :: x
  real(kind=8), allocatable, dimension (:) :: u_e  !Función análitica
  real(kind=8), allocatable, dimension (:) :: pi_n
  real(kind=8), allocatable, dimension (:) :: pi_n_p
  real(kind=8), allocatable, dimension (:) :: psi_n
  real(kind=8), allocatable, dimension (:) :: psi_n_p
  real(kind=8), allocatable, dimension (:) :: phi_n
  real(kind=8), allocatable, dimension (:) :: phi_n_p
  real(kind=8), allocatable, dimension (:) :: error
!!$  parametros desde un file


open(10,file='datos.dat')
  read(10,*) xmin,xmax,Nx,Nt,amp,sigma,x0,courant
close(10)
 

!!$ Asignar memoria a los arreglos

  allocate(x(0:Nx))
  allocate(u_e(0:Nx))
  allocate(pi_n(0:Nx))	   !!!!!!!!!!!
  allocate(pi_n_p(0:Nx))    !!!!!!!!!!
  allocate(psi_n(0:Nx))
  allocate(psi_n_p(0:Nx))  
  allocate(phi_n(0:Nx))		!!!!!!!!!!
  allocate(phi_n_p(0:Nx))	!!!!!!!!!!
   allocate(error(0:Nx))
 
!!$ Fijar el espaciado de la malla

  dx = ( xmax - xmin ) / dble( Nx ) !!$ Resolución espacial

  dt = courant * dx

  do i=0,Nx

     x(i)  = xmin + dble(i) * dx

  end do

!!$ Datos iniciales

  u_e = 0.5d0* amp * exp( - ( x - x0 )**2 / sigma**2 ) + 0.5d0*amp * exp( - ( x - x0 )**2 / sigma**2 ) + 1e-20

  pi_n = 0d0

  psi_n=(-2*amp*(x-x0)*exp(-(x-x0)**2/sigma**2))/sigma**2 + 1e-20

  phi_n= amp*exp(-(x-x0)**2/sigma**2) + 1e-20
 error=0d0
!!$ Se guardan los datos iniciales 

  t = 0.0d0

  call save1Ddata(Nx,t,x,u_e,'u_e',0)
  call save1Ddata(Nx,t,x,pi_n,'pi_n',0)
  call save1Ddata(Nx,t,x,psi_n, 'psi_n',0)
  call save1Ddata(Nx,t,x,phi_n, 'phi_n',0)
    call save1Ddata(Nx,t,x,error,'error',0)
!!$  loop de integración 

do l=1,Nt
  t = t + dt

	pi_n_p=pi_n
    psi_n_p=psi_n
    phi_n_p=phi_n   
    
     do i=1,Nx-1
        
        pi_n(i) = pi_n_p(i)-0.5*(dt/dx)*(psi_n_p(i+1)-psi_n_p(i-1))+ 0.5*((dt/dx)**2)*(psi_n_p(i+1)-2*psi_n_p(i) + psi_n_p(i-1))
        pi_n(i) = pi_n_p(i)-0.5*(dt/dx)*(psi_n_p(i+1)-psi_n_p(i-1))+ 0.5*((dt/dx)**2)*(psi_n_p(i+1)-2*psi_n_p(i) + psi_n_p(i-1))
        psi_n(i) = psi_n_p(i)+(pi_n_p(i+1)-pi_n_p(i-1))*(dt/(2*dx))
        
     end do

 do i=1,Nx-1 !Calculo la función con FTCS
  u_n(i)= u_n_p(i)- 0.5*(dt/dx)*(u_n_p(i+1)-u_n_p(i-1))+ 0.5*((dt/dx)**2)*(u_n_p(i+1)-2*u_n_p(i)+u_n_p(i-1))
 end do 
  

 
	u_e(Nx)=u_e(1)
    u_e(0)=u_e(Nx-1)
	u_n(Nx)=u_n(1)
	u_n(0)=u_n(Nx-1)

   u_e = amp * exp( - ( x - x0 - t)**2 / sigma**2 ) + 1e-20
    do i=0,Nx	
       error(i)= u_e(i)-phi_n(i)
    end do
    
   call save1Ddata(Nx,t,x,u_e,'u_e',1)
   call save1Ddata(Nx,t,x,u_n,'u_n',1)
   call save1Ddata(Nx,t,x,error,'error',1)
   
end do

end program Leadfrog

!
!
!
!
!

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
