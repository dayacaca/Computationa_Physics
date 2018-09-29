
subroutine BC

  use arrays
  use global_numbers

  implicit none

!!$ Outgoing boundary conditions  

  Left  = 0.50d0*( u(3,:) + u(2,:) ) !!$ Left mode
  Right = 0.50d0*( u(3,:) - u(2,:) ) !!$ Right mode

!!$ Extrapolating de values of the left an right modes
!!$ at the left and right boundaries
  
  Left(0) = 3.0d0*Left(1) - 3.0d0*Left(2) + Left(3)
  Right(0) = 3.0d0*Right(1) - 3.0d0*Right(2) + Right(3)

  Left(Nr) = 3.0d0*Left(Nr - 1) - 3.0d0*Left(Nr - 2) + Left(Nr - 3)
  Right(Nr) = 3.0d0*Right(Nr - 1) - 3.0d0*Right(Nr - 2) + Right(Nr - 3) 

!!$ Defining the boundary conditions
!!$ over the wave variables  
  
  u(3,0) = Left(0)
  u(2,0) = Left(0)

  u(3,Nr) = Right(Nr)
  u(2,Nr) = -Right(Nr)
  
end subroutine BC
