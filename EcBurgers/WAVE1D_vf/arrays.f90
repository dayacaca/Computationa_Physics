
  module arrays

  implicit none

  ! --------------- arrays for grid ---------------------
  
  real(kind=8), allocatable, dimension (:) :: r

  ! --------------- arrays for hydro --------------------

  ! Vector (phi,psi,pi)
  real(kind=8), allocatable, dimension (:) :: f
  real(kind=8), allocatable, dimension (:,:) :: u
  real(kind=8), allocatable, dimension (:,:) :: u_p
  real(kind=8), allocatable, dimension (:,:) :: rhs_u

  !real(kind=8), allocatable, dimension (:) :: Left
  !real(kind=8), allocatable, dimension (:) :: Right

  real(kind=8), allocatable, dimension (:) :: phi_exact
  real(kind=8), allocatable, dimension (:) :: error_phi
  
  ! --------------- arrays for RK4 ----------------------
  
  real(kind=8), allocatable, dimension (:,:) :: K1_u  
  real(kind=8), allocatable, dimension (:,:) :: K2_u
  real(kind=8), allocatable, dimension (:,:) :: K3_u
  real(kind=8), allocatable, dimension (:,:) :: K4_u

  ! ---------------- arrays for ICN ---------------------
  
  real(kind=8), allocatable, dimension (:,:) :: rhs_u1

  end module arrays
