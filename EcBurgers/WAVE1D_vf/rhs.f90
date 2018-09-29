
subroutine rhs

  use arrays
  use global_numbers

  implicit none

  integer i

  real(kind=8) idr      !paso espacial

  idr  = 1.0D0/dr



f=0.5*u(1,:)**2

  do i=1,Nr-1

     rhs_u(1,i) = 0.5d0*(f(i+1) - f(i-1))*idr
    ! rhs_u(3,i) = 0.5d0*(u(2,i+1) - u(2,i-1))*idr
     
  end do

     !rhs_u(1,:) = u(3,:)

end subroutine rhs


