
module arreglos

implicit none

  integer i,l,j
  integer N,d,code,Nx
real(kind=8), parameter :: pii = 4.0d0*ATAN(1.0d0)


REAL(kind=8) suma
REAL(kind=8) xmin, xmax, dx
REAL(kind=8) nm2, nm1
REAL(kind=8) lambda

REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: A ! --------------- arrays for matrix ---------------------
REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: x ! ---------arrays for puntos de colocación---------------
REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: b    !--------------- arrays for b vector---------------------
REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: aa  ! ------ arrays for a_n coeficientes de la expansión-------
integer, ALLOCATABLE, DIMENSION(:) :: indx     ! ------ arrays for LU subroutine-------
REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: XX
REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: error


end module arreglos
