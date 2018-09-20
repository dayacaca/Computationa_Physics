
subroutine check_parameters

  use arrays
  use global_numbers

  implicit none


  if ((integrator.ne.'RK2').and. &
       (integrator.ne.'RK3').and. &
       (integrator.ne.'RK4')) then
     print *, '================================================================================'
     print *, 'The integrator you selected --->',integrator,'<--- is not available'
     print *, '================================================================================'
     stop
  end if

end subroutine check_parameters
