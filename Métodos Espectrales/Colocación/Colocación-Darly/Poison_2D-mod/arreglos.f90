module arreglos 


  implicit none 

  integer i,j,k,l
  integer N,d,code

  integer Nx,Ny

  REAL(kind=8) suma
  real(kind=8), parameter :: pii = 4.0d0*ATAN(1.0d0)
  real(kind=8), parameter ::  xmin = -1.0d0
  real(kind=8), parameter ::  xmax = 1.0d0
  real(kind=8), parameter ::  ymin = -1.0d0
  real(kind=8), parameter :: ymax = 1.0d0 

  REAL(kind=8) dx,dy,nm2
 
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: A  ! matriz
  REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: x    ! puntos de colocación x
  REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: y    ! puntos de colocación y
  REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: b    ! vector fuente 
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: aa ! coeficientes de la expansión
  INTEGER, ALLOCATABLE, DIMENSION(:) :: indx
  REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: XX ! Malla 2D donde se va a evaluar la función
  REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: YY
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: error


end module arreglos 




