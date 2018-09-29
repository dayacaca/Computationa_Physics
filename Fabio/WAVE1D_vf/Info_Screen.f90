
  subroutine Info_Screen

  use global_numbers

  implicit none

 
  print *, '*****************************************'  
  print *, '***** Some numbers to the Sceen *********'
  print *, '*****************************************'
  print *, 'rmin  =',rmin
  print *, 'rmax  =',rmax
  print *, 'dr  =',dr
  print *, 'dt  =',dt
  print *, 'courant =',courant
  print *, 'Nt  =',Nt
  print *, 'Final time =',dt*dble(Nt)
  print *, 'every_1D =',every_1D
  print *, 'Number of time blocks =',Nt/every_1D
  print *, 'integrator =',integrator
  print *, '*****************************************'
  print *, '*****************************************'

  end subroutine Info_Screen

