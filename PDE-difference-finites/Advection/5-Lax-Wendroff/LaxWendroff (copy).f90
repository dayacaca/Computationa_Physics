!Programa que calcula númericamente la ecuación de advección
! por el método de diferencias finitas utilizando el esquema Leapfrog 
!el cual  utiliza dos niveles de tiempo n y n-1 para calcular la 
!función un tiempo posterior n+1, 
!u_n(i)= u_n_old(i) - (dt/dx)*(u_n_p(i-1)+ u_n_p(i+1))
!Como hacer que el programa entienda que ese u_n_old es el de tiempo anteiror
!u_n_old(i) = 0.5*(u_n_p(i+1) + u_n_p(i-1))- 0.5*(dt/dx)*(u_n_p(i+1)-u_n_p(i-1))
! 

program Leadfrog

  implicit none 

  integer i,l, j, m 
  integer Nx
  integer Nt

!!$ dotos de la malla

  real(kind=8) xmin 
  real(kind=8) xmax
  real(kind=8) dx
  real(kind=8) dt
  real(kind=8) t

!!$ Datos de la onda Gaussiana

  real(kind=8) amp
  real(kind=8) x0
  real(kind=8) sigma
  real(kind=8) courant

  real(kind=8), allocatable, dimension (:) :: x
  real(kind=8), allocatable, dimension (:) :: u_e  !Función análitica
  real(kind=8), allocatable, dimension (:) :: u_n  ! Función númerica en n+1 
  real(kind=8), allocatable, dimension (:) :: u_n_o !Fun. numer, en n-1 
  real(kind=8), allocatable, dimension (:) :: u_n_p   
  real(kind=8), allocatable, dimension (:) :: u_n_p_o
  
!!$ Preguntar parametros 


  ! 'The left end of the spatial domain': xmin
  ! 'The right end of the spatial domain' : xmax
  ! 'Number of points in the spatial domain': Nx
 ! 'Number of time iterations' : Nt
 !'Aplitude of the gaussian': amp
 !'Width of the initial gaussian' :sigma
 !'Center of the gaussian' :x0
  !'Courant Factor' : courant


 open(10,file='datos.dat')
  read(10,*) xmin,xmax,Nx,Nt,amp,sigma,x0,courant
 close(10)
 




!!$ Asignar memoria a los arreglos

  allocate(x(0:Nx))
  allocate(u_e(0:Nx))
  allocate(u_n(0:Nx))
  allocate(u_n_p(0:Nx))  
  allocate(u_n_o(0:Nx)) 
 allocate(u_n_p_o(0:Nx))  
 
!!$ Fijar el espaciado de la malla

  dx = ( xmax - xmin ) / dble( Nx ) !!$ Resolución espacial

  dt = courant * dx

  do i=0,Nx

     x(i)  = xmin + dble(i) * dx

  end do

!!$ Datos iniciales

  u_e = amp * exp( - ( x - x0 )**2 / sigma**2 ) + 1e-20

  u_n = amp * exp( - ( x - x0 )**2 / sigma**2 ) + 1e-20
 
  u_n_o = amp * exp( - ( x - x0 )**2 / sigma**2 ) + 1e-20
!!$ Se guardan los datos iniciales 

  t = 0.0d0

  call save1Ddata(Nx,t,x,u_e,'u_e',0)
  call save1Ddata(Nx,t,x,u_n,'u_n',0)
  call save1Ddata(Nx,t,x,u_n,'u_n_o',0)
  
  
!!$  loop de integración 

  do l=1,Nt
  t = t + dt
  u_n_p_o=u_n_o
  u_n_p=u_n 
    
  do i=1,Nx-1 !Calculo la función con FTCS
  u_n_o(i)= 0.5*(u_n_p_o(i+1)+ u_n_p_o(i-1))- 0.5*(dt/dx)*(u_n_p_o(i+1)-u_n_p_o(i-1))
  end do 
  
  
   !2do reciclado de variable

 do j=1,Nx-1
 u_n_p_o(j)=u_n_o(j-1) ! Este paso no estoy segura, ¿La función u_n_o no se, si hace lo que
 ! quiero que haga 
    do i=1,Nx-1
	u_n(i)= u_n_p_o(i)- (dt/dx)*(u_n_p(i+1)-u_n_p(i-1))
	end do 
 end do
 
 
 
   u_e = amp * exp( - ( x - x0 - t)**2 / sigma**2 ) + 1e-20
   call save1Ddata(Nx,t,x,u_e,'u_e',1)
   call save1Ddata(Nx,t,x,u_n,'u_n',1)
   call save1Ddata(Nx,t,x,u_n,'u_n_o',1)
  end do

end program Leadfrog

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
