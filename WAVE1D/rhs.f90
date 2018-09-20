
subroutine rhs

  use arrays
  use global_numbers

  implicit none

  integer i

  real(kind=8) idr

  idr  = 1.0D0/dr

  do i=1,Nr-1

     rhs_u(2,i) = 0.5d0*(u(3,i+1) - u(3,i-1))*idr

     rhs_u(3,i) = 0.5d0*(u(2,i+1) - u(2,i-1))*idr
     
  end do

     rhs_u(1,:) = u(3,:)

end subroutine rhs


