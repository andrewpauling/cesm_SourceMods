 !|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module forcing_hf_interior
! file heavily modified from forcing_s_interior.F90

!BOP
! !MODULE: forcing_hf_interior
!
! !DESCRIPTION:
!  Contains routines and variables used for determining the
!  interior freshwater forcing
!
! !REVISION HISTORY:
!  SVN:$Id: forcing_hf_interior.F90 12674 2008-10-31 22:21:32Z njn01 $
!
! !USES:

   use kinds_mod
   use blocks
   use domain
   use constants
   use io
   use forcing_tools
   use time_management
   use prognostic
   use grid
   use exit_mod

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_hf_interior,      &
              get_hf_interior_data, &
              set_hf_interior

! !PUBLIC DATA MEMBERS:

   real (r8), public ::       &! public for use in restart
      hf_interior_interp_last   ! time last interpolation done

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  internal module variables
!
!-----------------------------------------------------------------------

   real (r8), dimension(:,:,:,:,:), allocatable :: &
      HF_INTERIOR_DATA  ! data to use for interior freshwater flux

   integer (int_kind), dimension(:,:,:), allocatable :: &
      FWF_MAX_LEVEL ! maximum level for applying variable
                    ! interior freshwater flux

   real (r8), dimension(:,:,:), allocatable :: &
      FWF_MAX_LEVEL_REAL ! real number version of FWF_MAX_LEVEL


   real (r8), dimension(12) :: &
      hf_interior_data_time

   real (r8), dimension(20) :: &
      hf_interior_data_renorm   ! factors to convert data to model units

   real (r8) ::               &
      hf_interior_data_inc,    &! time increment between values of forcing data
      hf_interior_data_next,   &! time to be used for the next value of forcing data that is needed
      hf_interior_data_update, &! time when new forcing value to be added to interpolation set
      hf_interior_interp_inc,  &! time increment between interpolation
      hf_interior_interp_next   ! time next interpolation to be done

   integer (int_kind) ::            &
      hf_interior_interp_order,      &! order of temporal interpolation
      hf_interior_data_time_min_loc, &! index of x_data_time with the minimum forcing time
      hf_interior_max_level

   character (char_len) ::    &
      hf_interior_data_type,   &! keyword for period of forcing data
      hf_interior_filename,    &! name of file containing forcing data
      hf_interior_file_fmt,    &! format (bin or nc) of forcing file
      hf_interior_interp_freq, &! keyword for period of temporal interpolation
      hf_interior_interp_type, &!
      hf_interior_data_label,  &!
      hf_interior_formulation   !

   character (char_len), dimension(:), allocatable :: &
      hf_interior_data_names    ! names for required input data fields

   integer (int_kind), dimension(:), allocatable :: &
      hf_interior_bndy_loc,    &! location and field type for ghost
      hf_interior_bndy_type     !    cell updates



!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: init_hf_interior
! !INTERFACE:

 subroutine init_hf_interior

! !DESCRIPTION:
!  Initializes freshwater interior forcing by either calculating or
!  reading in the 3D freshwater.  Also performs initial book-keeping
!  concerning when new data is needed for the temporal interpolation
!  and when the forcing will need to be updated.
!
! !REVISION HISTORY:
!  same as module

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      n,                 &!
      nml_error           ! namelist i/o error flag

   character (char_len) :: &
      forcing_filename,    &! final filename for forcing data
      long_name             ! long name for input data

   type (datafile) :: &
      hf_int_data_file  ! data file descriptor for freshwater interior data

   type (io_field_desc) :: &
      hf_int_data_in    ! io field descriptor for input freshwater interior data

   type (io_dim) :: &
      i_dim, j_dim, &! dimension descriptors for horiz dimensions
      k_dim          ! dimension descriptor  for vertical levels

   namelist /forcing_hf_interior_nml/ hf_interior_formulation,          &
        hf_interior_data_type,        hf_interior_data_renorm,          &
        hf_interior_data_inc,         hf_interior_interp_type,          &
        hf_interior_interp_freq,      hf_interior_interp_inc,           &
        hf_interior_filename,             &
        hf_interior_file_fmt,         hf_interior_max_level

!-----------------------------------------------------------------------
!
!  read hf_interior namelist input after setting default values.
!  do not edit here, instead make changes to namelist in user_nl_pop2
!
!-----------------------------------------------------------------------

   hf_interior_formulation    = 'hosing'
   hf_interior_data_type      = 'none'
   hf_interior_interp_type    = 'linear'
   hf_interior_interp_freq    = 'every-timestep' ! How often to temporally interpolate fresh water flux data to current time.
   hf_interior_filename       = 'unknown-hf_interior'
   hf_interior_file_fmt       = 'nc'
   hf_interior_data_renorm    = c1 !  convert from X to Y

! These next three are not used but keep for now
   hf_interior_data_inc       = 1.e20_r8 ! Increment (hours) between forcing times if hf_interior_data_type='n-hour', which is not allowed. Hence this is not currently avail
   hf_interior_interp_inc     = 1.e20_r8 ! Increment (hours) between interpolation times if hf_interior_interp_freq='n-hour', which is not allowed. Hence this is not currently avail
   hf_interior_max_level      = 0        ! max depth level of input for hf_interior_data_type = 'analytic', which is not allowed


   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old', iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      !*** keep reading until find right namelist
      do while (nml_error > 0)
         read(nml_in, nml=forcing_hf_interior_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   end if

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
     call exit_POP(sigAbort,'ERROR: reading forcing_hf_interior_nml')
   endif

   call broadcast_scalar(hf_interior_formulation,       master_task)
   call broadcast_scalar(hf_interior_data_type,         master_task)
   call broadcast_scalar(hf_interior_interp_type,       master_task)
   call broadcast_scalar(hf_interior_interp_freq,       master_task)
   call broadcast_scalar(hf_interior_filename,          master_task)
   call broadcast_scalar(hf_interior_file_fmt,          master_task)
   call broadcast_array (hf_interior_data_renorm,       master_task)

   call broadcast_scalar(hf_interior_data_inc,          master_task)
   call broadcast_scalar(hf_interior_interp_inc,        master_task)
   call broadcast_scalar(hf_interior_max_level, master_task)

!-----------------------------------------------------------------------
!
!  convert data_type to 'monthly-calendar' if input is 'monthly'
!
!-----------------------------------------------------------------------

   if (hf_interior_data_type == 'monthly') &
       hf_interior_data_type = 'monthly-calendar'

!-----------------------------------------------------------------------
!
!  convert interp_type to corresponding integer value.
!
!-----------------------------------------------------------------------

   select case (hf_interior_interp_type)
   case ('nearest')
      hf_interior_interp_order = 1

   case ('linear')
      hf_interior_interp_order = 2

   case ('4point')
      hf_interior_interp_order = 4

   case default
      call exit_POP(sigAbort, &
         'init_hf_interior: Unknown value for hf_interior_interp_type')

   end select

!-----------------------------------------------------------------------
!
!  set values of the interior hf array (HF_INTERIOR_DATA)
!  depending on the type of hf data.
!
!-----------------------------------------------------------------------

   select case (hf_interior_data_type)

   case ('none')

      !*** no interior forcing, therefore no interpolation in time is
      !*** needed, nor are there any new values to be used.

      hf_interior_data_next = never
      hf_interior_data_update = never
      hf_interior_interp_freq = 'never'

   case ('annual')

      !*** annual mean climatological freshwater (read in from a file)
      !*** constant in time, therefore no new values will be needed

      allocate(HF_INTERIOR_DATA(nx_block,ny_block,km, &
                               max_blocks_clinic,1))

      allocate(hf_interior_data_names(1), &
               hf_interior_bndy_loc  (1), &
               hf_interior_bndy_type (1))

      HF_INTERIOR_DATA = c0
      hf_interior_data_names(1) = 'FRESHWATER' 
      hf_interior_bndy_loc  (1) = field_loc_center
      hf_interior_bndy_type (1) = field_type_scalar

      forcing_filename = hf_interior_filename

      hf_int_data_file = construct_file(hf_interior_file_fmt,            &
                                    full_name=trim(forcing_filename),  &
                                    record_length=rec_type_dbl,        &
                                    recl_words=nx_global*ny_global)

      call data_set(hf_int_data_file, 'open_read')

      i_dim = construct_io_dim('i', nx_global)
      j_dim = construct_io_dim('j', ny_global)
      k_dim = construct_io_dim('k', km)

      hf_int_data_in = construct_io_field(                      &
                         trim(hf_interior_data_names(1)),       &
                         dim1=i_dim, dim2=j_dim, dim3=k_dim,             &
                         field_loc = hf_interior_bndy_loc(1),   &
                         field_type = hf_interior_bndy_type(1), &
                         d3d_array = HF_INTERIOR_DATA(:,:,:,:,1))

      call data_set (hf_int_data_file, 'define', hf_int_data_in)
      call data_set (hf_int_data_file, 'read',   hf_int_data_in)
      call data_set (hf_int_data_file, 'close')
      call destroy_io_field(hf_int_data_in)
      call destroy_file(hf_int_data_file)

      if (hf_interior_data_renorm(1) /= c1) &
         HF_INTERIOR_DATA = HF_INTERIOR_DATA*hf_interior_data_renorm(1)

      hf_interior_data_next = never
      hf_interior_data_update = never
      hf_interior_interp_freq = 'never'

      if (my_task == master_task) then
         write(stdout,blank_fmt)
         write(stdout,'(a30,a)') ' Interior S Annual file read: ', &
                                 trim(forcing_filename)
      endif

   case ('monthly-equal','monthly-calendar')

      !*** monthly mean climatological interior freshwater
      !*** all 12 months are read from a file. interpolation order
      !*** may be specified with namelist input.

      allocate(HF_INTERIOR_DATA(nx_block,ny_block,km,    &
                               max_blocks_clinic,0:12))

      allocate(hf_interior_data_names(12), &
               hf_interior_bndy_loc  (12), &
               hf_interior_bndy_type (12))

      HF_INTERIOR_DATA = c0

      call find_forcing_times(           hf_interior_data_time,         &
                 hf_interior_data_inc,    hf_interior_interp_type,       &
                 hf_interior_data_next,   hf_interior_data_time_min_loc, &
                 hf_interior_data_update, hf_interior_data_type)

      forcing_filename = hf_interior_filename
      hf_int_data_file = construct_file(hf_interior_file_fmt,            &
                                    full_name=trim(forcing_filename),  &
                                    record_length=rec_type_dbl,        &
                                    recl_words=nx_global*ny_global)

      call data_set(hf_int_data_file, 'open_read')

      i_dim = construct_io_dim('i', nx_global)
      j_dim = construct_io_dim('j', ny_global)
      k_dim = construct_io_dim('k', km)

      do n = 1, 12
         if (n<10) then
           write(hf_interior_data_names(n),'(a11,i1)') 'FRESHWATER_',n
         else
           write(hf_interior_data_names(n),'(a11,i2)') 'FRESHWATER_',n
         endif
 
         hf_interior_bndy_loc (n) = field_loc_center
         hf_interior_bndy_type(n) = field_type_scalar

         hf_int_data_in = construct_io_field( &
                         trim(hf_interior_data_names(n)),       &
                         dim1=i_dim, dim2=j_dim, dim3=k_dim,             &
                         field_loc = hf_interior_bndy_loc(n),   &
                         field_type = hf_interior_bndy_type(n), &
                         d3d_array = HF_INTERIOR_DATA(:,:,:,:,n))

         call data_set (hf_int_data_file, 'define', hf_int_data_in)
         call data_set (hf_int_data_file, 'read',   hf_int_data_in)
         call destroy_io_field(hf_int_data_in)
      enddo

      call data_set (hf_int_data_file, 'close')
      call destroy_file(hf_int_data_file)

      if (hf_interior_data_renorm(1) /= c1) &
         HF_INTERIOR_DATA = HF_INTERIOR_DATA*hf_interior_data_renorm(1)

      if (my_task == master_task) then
         write(stdout,blank_fmt)
         write(stdout,'(a31,a)') ' Interior S Monthly file read: ', &
                                 trim(forcing_filename)
      endif

   case default

     call exit_POP(sigAbort, &
              'init_hf_interior: Unknown value for hf_interior_data_type')

   end select

!-----------------------------------------------------------------------
!
!  now check interpolation period (hf_interior_interp_freq) to set the
!    time for the next temporal interpolation (hf_interior_interp_next).
!
!  if no interpolation is to be done, set next interpolation time
!    to a large number so the interior S update test in routine
!    get_forcing_data will always be false.
!
!
!  if interpolation is to be done every timestep, set next interpolation
!    time to a large negative number so the interior S update
!    test in routine set_surface_forcing will always be true.
!
!-----------------------------------------------------------------------

   select case (hf_interior_interp_freq)

   case ('never')

     hf_interior_interp_next = never
     hf_interior_interp_last = never
     hf_interior_interp_inc  = c0

   case ('every-timestep')

     hf_interior_interp_next = always
     hf_interior_interp_inc  = c0

   case default

     call exit_POP(sigAbort, &
        'init_hf_interior: Unknown value for hf_interior_interp_freq')

   end select

   if(nsteps_total == 0) hf_interior_interp_last = thour00

!-----------------------------------------------------------------------
!
!  allocate and read in arrays used for variable interior hf
!    if necessary.
!
!-----------------------------------------------------------------------


      allocate(FWF_MAX_LEVEL(nx_block,ny_block,max_blocks_clinic), &
          FWF_MAX_LEVEL_REAL(nx_block,ny_block,max_blocks_clinic))

      forcing_filename = hf_interior_filename

      hf_int_data_file = construct_file(hf_interior_file_fmt,    &
                                    full_name=trim(forcing_filename),  &
                                    record_length=rec_type_dbl,        &
                                    recl_words=nx_global*ny_global)

      call data_set(hf_int_data_file, 'open_read')

      i_dim = construct_io_dim('i', nx_global)
      j_dim = construct_io_dim('j', ny_global)

      hf_int_data_in = construct_io_field('FWF_MAX_LEVEL',  &
                       dim1=i_dim, dim2=j_dim,                             &
                       field_loc = field_loc_center,             &
                       field_type = field_type_scalar,           &
                       d2d_array = FWF_MAX_LEVEL_REAL)      

      call data_set (hf_int_data_file, 'define', hf_int_data_in)
      call data_set (hf_int_data_file, 'read',   hf_int_data_in)
      FWF_MAX_LEVEL = nint(FWF_MAX_LEVEL_REAL)  
      call destroy_io_field(hf_int_data_in)

!-----------------------------------------------------------------------
!
!  echo forcing options to stdout.
!
!-----------------------------------------------------------------------

   hf_interior_data_label = 'Interior Freshwater Forcing'
   if (my_task == master_task) &
      write(stdout,"('Time-Dependent Interior Freshwater Forcing enabled')")
   call echo_forcing_options(hf_interior_data_type,                     &
                       hf_interior_formulation, hf_interior_data_inc,    &
                       hf_interior_interp_freq, hf_interior_interp_type, &
                       hf_interior_interp_inc, hf_interior_data_label)

!-----------------------------------------------------------------------
!EOC

 end subroutine init_hf_interior

!***********************************************************************
!BOP
! !IROUTINE: get_hf_interior_data
! !INTERFACE:

 subroutine get_hf_interior_data

! !DESCRIPTION:
!  Determines whether new interior forcing data is required and
!  reads the data if necessary.  Also interpolates data to current
!  time if required.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  check if new data is necessary for interpolation.  if yes, then
!    shuffle indices in SSS and hf_interior_data_time arrays. also
!    increment values of hf_interior_data_time_min_loc,
!    hf_interior_data_next and hf_interior_data_update. 
!
!-----------------------------------------------------------------------

   select case(hf_interior_data_type)

   case ('monthly-equal','monthly-calendar')

      hf_interior_data_label = 'HF_INTERIOR Monthly'
      if (thour00 >= hf_interior_data_update) then
         call update_forcing_data(             hf_interior_data_time,   &
                 hf_interior_data_time_min_loc, hf_interior_interp_type, &
                 hf_interior_data_next,         hf_interior_data_update, &
                 hf_interior_data_type,         hf_interior_data_inc,    &
                 HF_INTERIOR_DATA(:,:,:,:,1:12),hf_interior_data_renorm, &
                 hf_interior_data_label,        hf_interior_data_names,  &
                 hf_interior_bndy_loc,          hf_interior_bndy_type,   &
                 hf_interior_filename,          hf_interior_file_fmt)
      endif

      if (thour00 >= hf_interior_interp_next .or. nsteps_run==0) then
         call interpolate_forcing(HF_INTERIOR_DATA(:,:,:,:,0),          &
                                  HF_INTERIOR_DATA(:,:,:,:,1:12),       &
                 hf_interior_data_time,         hf_interior_interp_type, &
                 hf_interior_data_time_min_loc, hf_interior_interp_freq, &
                 hf_interior_interp_inc,        hf_interior_interp_next, &
                 hf_interior_interp_last,       nsteps_run)

         if (nsteps_run /= 0) hf_interior_interp_next = &
                              hf_interior_interp_next + &
                              hf_interior_interp_inc
      endif

   end select

!-----------------------------------------------------------------------
!EOC

 end subroutine get_hf_interior_data

!***********************************************************************
!BOP
! !IROUTINE: set_hf_interior
! !INTERFACE:

 subroutine set_hf_interior(k,this_block,HF_SOURCE)

! !DESCRIPTION:
!  Computes interior freshwater forcing term (DHF\_INTERIOR) using
!  updated and interpolated data.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      k     ! vertical level index

   type (block), intent(in) :: &
      this_block   ! block information for this block

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), intent(inout) :: &
      HF_SOURCE    ! source terms for freshwater forcing 

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::  &
      bid,                &! local block address for this block
      now                  ! index for location of interpolated data

   real (r8), dimension(nx_block,ny_block) :: &
      DHF_INTERIOR    ! interior freshwater restoring term for this block
                     ! and level

!-----------------------------------------------------------------------
!
!  do interior restoring if required (no restoring at surface for any)
!
!-----------------------------------------------------------------------

   if (hf_interior_data_type /= 'none' .and. k > 1) then

!-----------------------------------------------------------------------
!
!     set index for location of interpolated data.
!
!-----------------------------------------------------------------------

      select case(hf_interior_data_type)

      case('annual')
         now = 1

      case ('monthly-equal','monthly-calendar')
         now = 0

      end select

      bid = this_block%local_id

      DHF_INTERIOR = merge(HF_INTERIOR_DATA(:,:,k,bid,now),     &
                            c0, k <= FWF_MAX_LEVEL(:,:,bid))    

      !*** add restoring to any other source terms

      HF_SOURCE = HF_SOURCE - DHF_INTERIOR/salinity_factor*latent_heat_fusion &
			    /cp_sw/dz(k)/10/rho_fw*1/17.0*(tyear+26.0)

   endif ! k=1

!-----------------------------------------------------------------------
!EOC

 end subroutine set_hf_interior

!***********************************************************************

 end module forcing_hf_interior

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
