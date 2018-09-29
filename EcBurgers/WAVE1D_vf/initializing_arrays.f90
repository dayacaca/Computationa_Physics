
  subroutine initializing_arrays
 
  use arrays
  use global_numbers

  implicit none 


 ! --------------- arrays for grid ---------------------

   r  = 0.0d0
  
   ! --------------- arrays for wave  -------------------- 
   
   u     = 0.0d0
   u_p   = 0.0d0
   rhs_u = 0.0d0

   !Left  = 0.0d0
   !Right = 0.0d0

   phi_exact = 0.0d0
   error_phi = 0.0d0
   
 ! --------------- arrays for RK4 ----------------------  

   K1_u = 0.0d0
   K2_u = 0.0d0
   K3_u = 0.0d0
   K4_u = 0.0d0

 ! ---------------- arrays for ICN ---------------------
  
   rhs_u1 = 0.0d0  
   
  end subroutine initializing_arrays
