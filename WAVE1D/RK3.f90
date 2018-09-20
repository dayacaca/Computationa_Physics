!===================================================
! This is the subroutine of the integration method: 
! Runge-Kutta of 3-order 
!===================================================

subroutine RK3

  use arrays
  use global_numbers


  implicit none

  integer k
  real(kind=8) dt_temp

  u_p(:,:) = u(:,:)

  do k=1,3

     call rhs

     if (k.eq.1) then
        dt_temp = dt

        u(:,:) = u_p(:,:) + dt_temp*rhs_u(:,:)

     else if (k.eq.2) then       
        dt_temp = 0.25*dt

        u(:,:) = 0.75*u_p(:,:) + 0.25*u(:,:) + dt_temp*rhs_u(:,:)

     else 
        dt_temp = 2.0d0*dt/3.0d0

        u(:,:) = u_p(:,:)/3.0D0 + 2.0D0*u(:,:)/3.0D0 + dt_temp*rhs_u(:,:)

     end if

     call BC

  end do

end subroutine RK3
