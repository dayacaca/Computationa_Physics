!===================================================
! This is the subroutine of the integration method: 
! Runge-Kutta of 4-order 
!===================================================

subroutine RK4

  use arrays
  use global_numbers


  implicit none

  integer k
  real(kind=8) dt_temp

  u_p(:,:) = u(:,:)

  do k=1,4

     call rhs

     if (k.eq.1) then

        K1_u(:,:) = rhs_u(:,:)

        dt_temp = 0.5D0*dt

        u(:,:) = u_p(:,:) + dt_temp*K1_u(:,:)

     else if (k.eq.2) then       

        dt_temp = 0.5D0*dt

        K2_u(:,:) = rhs_u(:,:)

        u(:,:) = u_p(:,:) + dt_temp*K2_u(:,:)

     else if (k.eq.3) then       

        dt_temp = dt

        K3_u(:,:) = rhs_u(:,:)

        u(:,:) = u_p(:,:) + dt_temp*K3_u(:,:)


     else if (k.eq.4) then

        dt_temp = dt

        K4_u(:,:) = rhs_u(:,:)

        u(:,:) = u_p(:,:) + (dt/6.0d0)*(K1_u(:,:) + 2.0d0*K2_u(:,:) &
             + 2.0d0*K3_u(:,:) + K4_u(:,:) )

     end if

     call BC

  end do

end subroutine RK4
