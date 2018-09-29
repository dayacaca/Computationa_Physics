
  subroutine allocate

  use arrays
  use global_numbers 

  implicit none

  allocate( r(0:Nr))

  ! Vector
  allocate(  u(1:nvars,0:Nr))
  allocate(  u_p(1:nvars,0:Nr))
  allocate(f(0:Nr))
  allocate(rhs_u(1:nvars,0:Nr))
  
 ! allocate(Left(0:Nr))
 ! allocate(Right(0:Nr))

  allocate(phi_exact(0:Nr))
  allocate(error_phi(0:Nr))
    
 ! ----- arrays for RK4 ---------------
  
  allocate(K1_u(1:nvars,0:Nr))
  allocate(K2_u(1:nvars,0:Nr))
  allocate(K3_u(1:nvars,0:Nr))
  allocate(K4_u(1:nvars,0:Nr))

! ----- arrays for ICN---------------

  allocate(rhs_u1(1:nvars,0:Nr))


  end subroutine allocate 
