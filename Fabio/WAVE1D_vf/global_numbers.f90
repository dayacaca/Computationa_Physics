 
module global_numbers

  ! Grid Variables
  real(kind=8) rmin, rmax
  real(kind=8) dr 
  real(kind=8) t, dt, courant
  integer Nr, Nt 
  integer every_0D, every_1D
  integer res_num
  real(kind=8) pii

  ! Hydro Variables
  real(kind=8) amp
  real(kind=8) r0
  real(kind=8) sigma
  integer nvars

  ! Integrators
  character(len=20) :: integrator

  !Parameters to compute the CPU time
  real (kind=8) cpu_it
  real (kind=8) cpu_ft
  real (kind=8) cpu_t

end module global_numbers

