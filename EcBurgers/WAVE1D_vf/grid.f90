
subroutine grid

  use arrays
  use global_numbers

  implicit none

  integer i

  dr  = ( rmax - rmin ) / dble(Nr)

  do i=0,Nr
           r(i)  = rmin + dble(i) * dr
  end do

  dt = courant * dr

end subroutine grid
