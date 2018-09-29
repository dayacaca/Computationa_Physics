!===================================================
! This is the subroutine of the integration method: 
! Runge-Kutta of 2-order 
!===================================================

subroutine RK2

  use arrays
  use global_numbers


  implicit none

  integer k,n
  real(kind=8) dt_temp

  u_p(:,:) = u(:,:) 

  do k=1,2

     call rhs

     if (k.eq.1) then
        dt_temp = dt

        u(:,:) = u_p(:,:) + dt_temp*rhs_u(:,:)
        
     else if (k.eq.2) then       
        dt_temp = 0.5D0*dt
        
        u(:,:) = 0.5D0*(u_p(:,:) + u(:,:)) + dt_temp*rhs_u(:,:)
        
     end if

     call BC

  end do
  
end subroutine RK2
