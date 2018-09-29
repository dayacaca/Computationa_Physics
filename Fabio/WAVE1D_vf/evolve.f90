
subroutine evolve

  use arrays
  use global_numbers

  implicit none

  integer i,j,k,l
  integer num_steps

  call allocate
  call initializing_arrays
  call grid
  call Info_Screen

  t = 0.0d0

  print *,'----------------------------'
  print *,'|  Time step  |    Time    |'
  print *,'----------------------------'
  write(*,"(A5,I6,A6,ES9.2,A3)") ' |   ',0,'    | ',t,'  |'

  call initial
  call exact

  call save1Ddata_x(u(1,:),'phi',0)
  call save1Ddata_x(u(2,:),'psi',0)
  call save1Ddata_x(u(3,:),'pi',0)

  call save1Ddata_x(Left,'Left',0)
  call save1Ddata_x(Right,'Right',0)

  call save1Ddata_x(phi_exact,'phi_exact',0)
  call save1Ddata_x(error_phi,'error_phi',0)
  
  ! ****************************************************************
  ! ****************************************************************
  ! ******************      Integration Loop        ****************
  ! ****************************************************************
  ! ****************************************************************

  do l=1,Nt
     t = t + dt

     if (mod(l,every_1D).eq.0) then
        call cpu_time(cpu_t)
        write(*,"(A5,I6,A6,ES9.2,A3)") ' |   ',l,'    | ',t,'  |'
        print *,'cpu time = ',cpu_t-cpu_it,'secs'
     end if

     if (integrator.eq.'RK2') then
        call RK2
     else if (integrator.eq.'RK3') then
        call RK3
     else if (integrator.eq.'RK4') then  
        call RK4
     else
        print*, 'The integrator has not been implemented yet'
        stop
     end if

     call exact 

     if (mod(l,every_0D).eq.0) then

     end if

     if (mod(l,every_1D).eq.0) then

        call save1Ddata_x(u(1,:),'phi',1)
        call save1Ddata_x(u(2,:),'psi',1)
        call save1Ddata_x(u(3,:),'pi',1)

        call save1Ddata_x(Left,'Left',1)
        call save1Ddata_x(Right,'Right',1)

        call save1Ddata_x(phi_exact,'phi_exact',1)
        call save1Ddata_x(error_phi,'error_phi',1)

     end if

  end do

end subroutine evolve

