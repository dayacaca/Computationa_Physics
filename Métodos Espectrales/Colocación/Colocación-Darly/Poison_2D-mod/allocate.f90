subroutine allocate
use arreglos

implicit none 

  
  print *, 'Number of colocation points'
  read(*,*) N
  print *

  print *, 'Number of points on the x direction'
  read(*,*) Nx
  print *

  print *, 'Number of points on the y direction'
  read(*,*) Ny
  print *

  ALLOCATE(A(1:(N+1)*(N+1),1:(N+1)*(N+1)))
  ALLOCATE(x(0:N))
  ALLOCATE(y(0:N))
  ALLOCATE(b(1:(N+1)*(N+1)))
  ALLOCATE(indx(1:(N+1)*(N+1)))
  ALLOCATE(aa(0:N,0:N))
  ALLOCATE(XX(0:Nx))
  ALLOCATE(YY(0:Ny))
  ALLOCATE(error(0:Nx,0:Ny))

end subroutine allocate





