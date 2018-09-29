
  subroutine save1Ddata_x(fval,base_name,first_index)

  use arrays
  use global_numbers

  implicit none

  character(len=20) filestatus
  logical firstcall
  data firstcall / .true. /
  save firstcall

  character(len=*), intent(IN) :: base_name
  real(kind=8), dimension(0:Nr), intent(IN) :: fval

  character(len=256) :: filename

  integer i,j,l,first_index, auxcounter,after,before

  if (res_num.eq.1) then
    filename = base_name // '_1.xl'
  else if (res_num.eq.2) then
    filename = base_name // '_2.xl'
  else if (res_num.eq.3) then
    filename = base_name // '_3.xl'
  else if (res_num.eq.4) then
    filename = base_name // '_4.xl'
  else if (res_num.eq.5) then
    filename = base_name // '_5.xl'
  end if

  if (first_index.eq.0) then
     filestatus = 'replace'
  else
     filestatus = 'old'
  end if


  if (filestatus=='replace') then
     open(1,file=filename,status=filestatus)
  else
     open(1,file=filename,status=filestatus,position='append')
  end if
  write(1,*) ''
  write(1,*) '#Time =',t

  do i=0,Nr,2**(res_num-1)
    write(1,*) r(i),fval(i)  
  end do

  write(1,*)
  close(1)

  end subroutine save1Ddata_x

