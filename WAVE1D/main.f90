 
program main

  use global_numbers

  implicit none

  integer i,j,k,Nrr,Ntt
  integer every_0Dt, every_1Dt

  Namelist /BHtoy_Input/ rmin, rmax, & 
       res_num, Nrr, courant, Ntt, &
       every_0Dt, every_1Dt, &
       integrator, &
       nvars, r0, Amp, sigma

  open (3, file='input.par', status = 'old' )
  read (3, nml = BHtoy_Input)
  close(3)

  pii   = 4.0d0*atan(1.0d0)

  Nr  = 2**(res_num-1)*Nrr
  Nt  = 2**(res_num-1)*Ntt

  every_0D = 2**(res_num-1)*every_0Dt
  every_1D = 2**(res_num-1)*every_1Dt

  call check_parameters

  call cpu_time(cpu_it)
  call evolve
  call cpu_time(cpu_ft)

  call Info_Screen

  print *, '=================================='
  print *, '======  CAFE_WAVE_1D_code  ======='
  print *, '=================================='

end program main


  
