module mcm_database
!
! Module containing main data base for code
!
! Variables declared in this module should all start mcm_
!  This is to ensure that there are no variable name clashes in the linked version.
!
implicit none
save
!
! define minimum floating point accuracy
integer, parameter :: real_acc = selected_real_kind(P=15,R=50)
integer, parameter :: d=real_acc
!
!---------------------------------------------------------------------------
! Global variables
!
logical :: mcm_ctrlc
!
integer :: mcm_np, mcm_nummat, mcm_ndim, mcm_nstressp, mcm_nvelocp, mcm_max_np
integer :: mcm_max_ngp, mcm_ngp
integer :: mcm_ssp, mcm_esp, mcm_svp, mcm_evp
integer :: mcm_timestep, mcm_status_interval
integer :: mcm_restart_interval, mcm_run_restart
integer :: mcm_axopt, mcm_disctype
integer :: mcm_maxcont, mcm_maxnbr
integer :: mcm_next_restart, mcm_next_run_restart
!
real(kind=real_acc) :: mcm_endtime, mcm_tssf, mcm_itss
real(kind=real_acc) :: mcm_ptime, mcm_dt, mcm_dtold, mcm_critts, mcm_init_ts
real(kind=real_acc) :: mcm_coord_maxmin(2,3)
real(kind=real_acc) :: pi
!
!--------------------------------------------------------------------------
! History variables
!
real(kind=real_acc) :: mcm_thermale, mcm_kinetice, mcm_internale, mcm_totale
!
!---------------------------------------------------------------------------
! Options
!
logical :: mcm_repulsive_force
logical :: mcm_boundary
logical :: mcm_drelax
integer :: mcm_contacttype, mcm_krtype, mcm_init_h_opt, mcm_init_v_opt, mcm_init_rhoe
integer :: mcm_massopt, mcm_tcrit_opt, mcm_veloc_opt, mcm_h_opt
integer :: mcm_nlcur, mcm_nthpx, mcm_nthpy, mcm_nthpz
integer :: mcm_ncontmats, mcm_cont_opt
real(kind=real_acc) :: mcm_drelax_scale
!
!---------------------------------------------------------------------------
! Output
!
integer :: mcm_state_opt
integer :: mcm_thnode, mcm_num_transducer
real(kind=real_acc) :: mcm_stpltime, mcm_thpltime
!
integer :: mcm_istate
real(kind=real_acc) :: mcm_nextsttime, mcm_nextthtime
!
integer, dimension(:), allocatable :: mcm_thnodeid
! transducers
integer, dimension(:), allocatable :: mcm_trmat, mcm_num_tr_nbr
real(kind=real_acc), dimension(:,:), allocatable :: mcm_trans_x
integer, dimension(:,:), allocatable :: mcm_tr_nbr
! MCMGUI variables
integer :: mcm_out_cols(2)
real(kind=real_acc) :: mcm_maxmin(2,30)
character(len=40)   :: mcm_variablename(30)
! Ensight variables
character(len=25)   :: mcm_filename_C
character(len=25)   :: mcm_extension
character(len=25)   :: mcm_filename_PCT
character(len=25)   :: mcm_filename_mass, mcm_filename_Temp, mcm_filename_hoop
character(len=25)   :: mcm_filename_Bound, mcm_filename_neighbours, mcm_filename_eps
character(len=25)   :: mcm_filename_strain, mcm_filename_energy, mcm_filename_smooth
character(len=25)   :: mcm_filename_sig11, mcm_filename_sig22, mcm_filename_sig33,&
                       mcm_filename_sig12, mcm_filename_sig13, mcm_filename_sig23, mcm_filename_dens
character(len=25)   :: mcm_filename_dis, mcm_filename_pres, mcm_filename_vel, mcm_filename_acc, mcm_filename_material, mcm_filename_vtk
character(len=25)   :: mcm_file_time
! Ensight case file variables
character(len=25)   :: mcm_filename_sig, mcm_filename_rho, mcm_filename_snd, mcm_filename_egy
character(len=25)   :: mcm_filename_mas, mcm_filename_bnd, mcm_filename_nbr, mcm_filename_tmp
REAL(kind=real_acc) :: mcm_timenum(2000)
! d3plot values
integer :: mcm_d3buffcnt, mcm_d3count, mcm_d3file
real(4), dimension(512) :: mcm_d3plot_buffer
!
!---------------------------------------------------------------------------
! Neighbours and Neighbour searching
!
real(kind=real_acc) :: mcm_hmax
real(kind=real_acc),dimension(3) :: mcm_gridsize
!
integer, dimension(3) :: mcm_gridlim
integer, dimension(:,:), allocatable :: mcm_nbrlist, mcm_contlist
integer, dimension(:,:,:), allocatable :: mcm_llgrid
!
real(kind=real_acc), dimension(3) :: mcm_gridmin
!
!---------------------------------------------------------------------------
! Symmetry + periodic planes
!
! New master arrays
!
integer, dimension(3) :: mcm_boundary_type
integer, dimension(2,3) :: mcm_boundary_code
real(kind=real_acc), dimension(2,3) :: mcm_boundary_x
!
integer :: mcm_g_maxnbr, mcm_g_maxcont
integer, dimension(:,:), allocatable :: mcm_g_nbrlist, mcm_g_contlist
!
! New multiply and add arrays
!
real(kind=real_acc), dimension(3) :: mcm_xmin_mult, mcm_xmax_mult
real(kind=real_acc), dimension(3) :: mcm_ymin_mult, mcm_ymax_mult
real(kind=real_acc), dimension(3) :: mcm_zmin_mult, mcm_zmax_mult
real(kind=real_acc), dimension(3) :: mcm_xmin_add, mcm_xmax_add
real(kind=real_acc), dimension(3) :: mcm_ymin_add, mcm_ymax_add
real(kind=real_acc), dimension(3) :: mcm_zmin_add, mcm_zmax_add
!
!---------------------------------------------------------------------------
! File handling
!
character(len=12) :: mcm_filein, mcm_fileout
character(len=16) :: mcm_filerestart
character(len=78) :: mcm_title
!
integer :: mcm_filelen(3)
integer :: mcm_openfiles(60)
!
!---------------------------------------------------------------------------
! Load Curves
!
REAL(kind=real_acc), dimension(:), allocatable :: mcm_npc, mcm_plc
!
!---------------------------------------------------------------------------
! Contact Data
!
REAL(kind=real_acc), dimension(:), allocatable :: mcm_k_cont, mcm_n_cont
!
!---------------------------------------------------------------------------
! Temporary variables
!
REAL(kind=real_acc) :: s11, s12, s13, s21, s22, s23, s31, s32, s33 ! rotation matrix - 25-04-02 - Tom - mat52 mat53
REAL(kind=real_acc) :: b11, b12, b13, b21, b22, b23, b31, b32, b33 ! rotation matrix - 25-04-02 - Tom - mat52 mat53
!
!----------------------------------------------------------------------------
! Base accelerations
!
logical :: mcm_baseaccel
real (kind=real_acc), dimension(3) :: mcm_base_a
!
!
!===========================================================================
!
! Define objects as f90 types
!
type particle
  logical :: active                                 ! Particle deletion flag
  integer :: mat                                    ! Material model id
  integer :: dispbc                                 ! Displacement boundary condition
  integer :: llpointer                              ! linked list pointer
  integer :: delpointer                             ! inactive particle linked list pointer
  integer :: nnbr                                   ! number of neighbour particles
  integer :: ncont                                  ! number of neighbour contact particles
  
  integer :: g_nnbr                                 ! number of ghost neighbours
  integer :: g_ncont                                ! number of ghost contact neighbours
  integer :: boundary                               ! flag for boundary particles
  integer :: nsym                                   ! flag for symmetry planes
!  integer :: n_symnbr                               ! number of symmetry neighbours
!  integer :: n_symcont                              ! number of symmetry contact neighbours
  real(kind=real_acc) :: mass                       ! Mass
  real(kind=real_acc) :: h                          ! Smoothing length
  real(kind=real_acc) :: hold                       ! Smoothing length for previous timestep
  real(kind=real_acc) :: h0                         ! Initial smoothing length (fixed during initialisation)
  real(kind=real_acc) :: rho                        ! Density
  real(kind=real_acc) :: rho0                       ! Initial density (fixed during initialisation)
  real(kind=real_acc) :: rhoold
  real(kind=real_acc) :: p                          ! pressure
  real(kind=real_acc) :: e                          ! TOTAL particle internal energy
  real(kind=real_acc) :: etry                       ! trial particle internal energy
  real(kind=real_acc) :: einc                       ! particle internal energy increment
  real(kind=real_acc) :: c                          ! speed of sound
  real(kind=real_acc) :: vabs                       ! Absolute value of velocity
  real(kind=real_acc) :: tracerod                   ! Trace of the rate-of-deformation tensor
  real(kind=real_acc) :: pcut                       ! Cutoff pressure
  real(kind=real_acc) :: p_cut                      ! Cutoff pressure at different time
  real(kind=real_acc) :: efps                       ! Effective plastic strain
  real(kind=real_acc) :: fail                       ! Failure, varies between undamaged (1.0) and failed (0.0)
  real(kind=real_acc) :: temper                     ! Temperature, unused at this time
  real(kind=real_acc) :: mindist                    ! Distance to nearest neighbour
  real(kind=real_acc) :: critts                     ! Particle critical timestep
  real(kind=real_acc) :: capa                       ! Added by M. Topalovic capa for mat model 25
  integer :: ptype                                  ! Added by M. Topalovic ptype for mat model 25 - Tijana added, for control print
  real(kind=real_acc), dimension(6) :: epx          ! material history variables
  real(kind=real_acc) :: alfa(3,3)					! back stress (kinematic hardening)
  ! vectors
  real(kind=real_acc), dimension(3) :: x            ! Current particle coordinates
  real(kind=real_acc), dimension(3) :: xzero        ! Initial particle coordinates
  real(kind=real_acc), dimension(3) :: v            ! Current particle velocity
  real(kind=real_acc), dimension(3) :: smooth_v     ! Interpolated particle velocity, currently used only for XSPH
  real(kind=real_acc), dimension(3) :: a            ! Particle acceleration
  real(kind=real_acc), dimension(3) :: bndnorm      ! Unit normal vector to surface for boundary particles
  real(kind=real_acc), dimension(3) :: repulsion	! contact acceleration
  ! tensors
  real(kind=real_acc), dimension(3,3) :: rod        ! Rate-of-deformation tensor
  real(kind=real_acc), dimension(3,3) :: spin       ! Spin tensor
  real(kind=real_acc), dimension(3,3) :: sigma      ! stress tensor
  real(kind=real_acc), dimension(3,3) :: s	        ! deviatoric stress tensor
  real(kind=real_acc), dimension(3,3) :: q          ! artificial viscosity stress tensor
  real(kind=real_acc), dimension(3,3) :: qold       ! artificial viscosity stress tensor from previous timestep
  real(kind=real_acc), dimension(3,3) :: qq         ! global to current configuration rotation matrix 
  real(kind=real_acc), dimension(3,3) :: qr         ! material axes rotation matrix for orthotropic materials
end type particle
!--------------------------------------------------
type ghost_particle
  integer :: par                                    ! Real particle that ghost is created from
  integer :: mat                                    ! Material model id
  integer :: llpointer                              ! linked list pointer
  real(kind=real_acc) :: mass                       ! Mass
  real(kind=real_acc) :: h                          ! Smoothing length
  real(kind=real_acc) :: hold                       ! Smoothing length for previous timestep
  real(kind=real_acc) :: rho                        ! Density
  ! vectors
  real(kind=real_acc), dimension(3) :: x            ! Current particle coordinates
  real(kind=real_acc), dimension(3) :: xzero        ! 'Initial' particle coordinates
  real(kind=real_acc), dimension(3) :: v            ! Current particle velocity
  ! tensors
  real(kind=real_acc), dimension(3,3) :: sigma      ! stress tensor
  real(kind=real_acc), dimension(3,3) :: q          ! artificial viscosity stress tensor
end type ghost_particle
!--------------------------------------------------
type material
  character(len=6) :: head(2,12)                    ! Material title
  integer :: model                                  ! Material model number
  integer :: eos                                    ! Material equation of state model
  integer :: visc_type                              ! Material artificial viscosity type
  real(kind=real_acc) :: rho                        ! density
  real(kind=real_acc) :: g                          ! shear modulus
  real(kind=real_acc) :: av_l                       ! linear artificial viscosity coefficient
  real(kind=real_acc) :: av_q                       ! quadratic artificial viscosity coefficient
  real(kind=real_acc) :: mass                       ! total material mass, used if mcm_massopt = 0
  real(kind=real_acc) :: h                          ! smoothin length, used if mcm_init_h_opt = 0
  real(kind=real_acc) :: rho_min                    ! minimum density limit factor
  real(kind=real_acc) :: rho_max                    ! maximum density limit factor
  real(kind=real_acc) :: ke                         ! Total kinetic energy of material
  real(kind=real_acc) :: ie                         ! Total internal energy of material
  REAL(kind=real_acc) :: xcf						! Contact force in X-direction
  REAL(kind=real_acc) :: ycf						! Contact force in Y-direction
  REAL(kind=real_acc) :: zcf						! Contact force in Z-direction
  real(kind=real_acc), dimension(3) :: mom          ! Momentun vector
  real(kind=real_acc), dimension(48) :: strinput    ! Constitutive coefficients
  real(kind=real_acc), dimension(48) :: strinput2   ! Additional constitutive coefficients
  real(kind=real_acc), dimension(24) :: eosinput    ! EoS coefficients
end type material

type(particle), dimension(:), allocatable :: par
type(ghost_particle), dimension(:), allocatable :: gpar
type(material), dimension(:), allocatable :: mcm_mat
!
contains
 !
 !=================================================
 !
 subroutine mcm_allocate_memory
 ! allocate memory for particle types
 ! imput number of points only
 implicit none
 integer :: ierr, i
 !
 ! Allocate material memory
 !
 allocate (mcm_mat(mcm_nummat),stat=ierr)
 if(ierr.ne.0) then
  write(*,'("Allocation error for material data")')
  stop
 endif
 !
 ! Allocate particle memory
 !
 allocate (par(mcm_max_np),stat=ierr)
 if(ierr.ne.0) then
  write(*,'("Allocation error for particle data")')
  stop
 endif
 !
 !-------------------------------------------------------
 ! Allocate other variables
 !
 ! Array for list of time history particle ID's
 if(mcm_thnode.gt.0) then
  allocate (mcm_thnodeid(mcm_thnode),stat=ierr)
  if(ierr.ne.0) then
   write(*,'("Allocation error for time history data")')
   stop
  endif
 endif
 ! Transducers
 if(mcm_num_transducer.gt.0) then
  allocate (mcm_trmat(mcm_num_transducer),mcm_trans_x(3,mcm_num_transducer), &
            mcm_num_tr_nbr(mcm_num_transducer),stat=ierr)
  if(ierr.ne.0) then
   write(*,'("Allocation error for pressure transducers")')
   stop
  endif
 endif
 !
 if(mcm_contacttype.gt.0) then
  allocate (mcm_k_cont(mcm_nummat), mcm_n_cont(mcm_nummat), stat=ierr)
  if(ierr.ne.0) then
   write(*,'("Allocation error for contact data")')
   stop
  endif
 endif
 !
 end subroutine mcm_allocate_memory
 !
end module mcm_database