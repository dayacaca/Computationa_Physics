
subroutine initial

  use arrays
  use global_numbers

  implicit none

  integer i,j,k

  u(1,:) = amp*exp((-(r-r0)**2)/sigma**2 ) !!$ phi  
  u(2,:) = -2.0d0*(r-r0)*u(1,:)/sigma**2   !!$ psi
  u(3,:) = 0.0d0  !!$ pi

  Left  = 0.50d0*( u(3,:) + u(2,:) ) !!$ Left mode
  Right = 0.50d0*( u(3,:) - u(2,:) ) !!$ Right mode
  
end subroutine initial







