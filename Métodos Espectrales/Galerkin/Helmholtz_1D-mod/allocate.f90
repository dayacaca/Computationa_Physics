subroutine allocate

use arreglos

implicit none

open(10,file='datos.dat')
read(10,*) N,Nx
close(10)

ALLOCATE(A(1:N+1,1:N+1))
ALLOCATE(x(0:N))
ALLOCATE(b(1:N+1))
ALLOCATE(indx(1:N+1))
ALLOCATE(aa(0:N))
ALLOCATE(XX(0:Nx))
ALLOCATE(error(0:Nx))


end subroutine allocate
