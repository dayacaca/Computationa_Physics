
subroutine exact

  use arrays
  use global_numbers

  implicit none

  phi_exact = 0.50d0*amp*exp( -(r - r0 - t)**2/sigma**2 ) & 
            + 0.50d0*amp*exp( -(r - r0 + t)**2/sigma**2 )

  error_phi = phi_exact - u(1,:)
  
end subroutine exact
