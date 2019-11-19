!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module preformed_tracers_mod

!BOP
! !MODULE: preformed_tracers_mod
!
! !DESCRIPTION: tracks preformed PO4 and ALK (once tracer arrives at the surface it is considered "preformed"; total tracer - preformed tracer = regenerated tracer)
!
! !REVISION HISTORY:
!  SVN:$Id$

! !USES:

   use POP_KindsMod
   use POP_IOUnitsMod
   use POP_ErrorMod

   use blocks, only: nx_block, ny_block
   use domain_size, only: max_blocks_clinic, km, nt
   use domain, only: nblocks_clinic
   use exit_mod, only: sigAbort, exit_POP
   use communicate, only: my_task, master_task
   use kinds_mod
   use constants, only: c0, c1, char_blank, delim_fmt
   use io, only: data_set
   use io_types, only: stdout, nml_in, nml_filename
   use io_tools, only: document
   use passive_tracer_tools, only: ind_name_pair, rest_read_tracer_block
   use passive_tracer_tools, only: file_read_tracer_block, tracer_read
   implicit none
   private

! !PUBLIC MEMBER FUNCTIONS:

   public :: preformed_tracer_cnt,        &
             preformed_tracer_init,              &
             preformed_tracer_reset

!EOP
!BOC

!-----------------------------------------------------------------------
!  module variables required by passive_tracer
!-----------------------------------------------------------------------

   integer(int_kind), parameter :: &
      preformed_tracer_cnt = 3

!-----------------------------------------------------------------------
!  relative tracer indices
!-----------------------------------------------------------------------

   integer(int_kind), parameter :: &
      preformed_alk_ind = 1     ! preformed alk index

   integer(int_kind), parameter :: &
      preformed_po4_ind = 2     ! preformed po4 index

   integer(int_kind), parameter :: &
      preformed_o2_ind = 3     ! preformed o2 index
!-----------------------------------------------------------------------
!  derived type & parameter for tracer index lookup
!-----------------------------------------------------------------------

   type(ind_name_pair), dimension(preformed_tracer_cnt) :: &
      ind_name_table = (/ ind_name_pair(preformed_alk_ind, 'ALK_preformed'), ind_name_pair(preformed_po4_ind, 'PO4_preformed'), &
      ind_name_pair(preformed_o2_ind, 'O2_preformed') /)

!EOC
!*****************************************************************************

contains

!*****************************************************************************
!BOP
! !IROUTINE: preformed_tracer_init
! !INTERFACE:

 subroutine preformed_tracer_init(preformed_tracers_ind_begin, init_ts_file_fmt, read_restart_filename, &
                                  tracer_d_module, TRACER_MODULE, errorCode)

! !DESCRIPTION:
!  Initialize preformed tracer module. This involves setting metadata, reading
!  the module namelist and setting initial conditions.
!
! !REVISION HISTORY:
!  same as module

! !USES:

   use broadcast, only: broadcast_scalar
   use prognostic, only: curtime, oldtime, tracer_field
   use grid, only: KMT, n_topo_smooth, fill_points
   use passive_tracer_tools, only : file_read_single_tracer

! !INPUT/OUTPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      preformed_tracers_ind_begin         ! starting index of preformed tracers in global tracer array
                                          ! passed through to rest_read_tracer_block

   character (*), intent(in) ::  &
      init_ts_file_fmt,          & ! format (bin or nc) for input file
      read_restart_filename        ! file name for restart file

   type (tracer_field), dimension(preformed_tracer_cnt), intent(inout) :: &
      tracer_d_module   ! descriptors for each tracer

   real(r8), dimension(nx_block,ny_block,km,preformed_tracer_cnt,3,max_blocks_clinic), &
      intent(inout) :: TRACER_MODULE

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character(*), parameter :: subname = 'preformed_tracer_mod:preformed_tracer_init'

   character(char_len) :: &
      init_preformed_tracers_option,          & ! option for initialization of preformed tracers
      init_preformed_tracers_init_file

   logical(log_kind) :: &
      lnml_found             ! Was preformed_tracers_nml found ?

   integer(int_kind) :: &
      n,                   & ! index for looping over tracers
      k,                   & ! index for looping over depth levels
      iblock,              & ! index for looping over blocks
      nml_error              ! namelist i/o error flag

   type(tracer_read), dimension(preformed_tracer_cnt) :: tracer_inputs   ! metadata about file to read

   character (char_len) ::  &
      preformed_tracer_restart_filename  ! modified file name for restart file

   namelist /preformed_tracers_nml/init_preformed_tracers_option, init_preformed_tracers_init_file


!-----------------------------------------------------------------------
!  initialize tracer_d values
!-----------------------------------------------------------------------

   errorCode = POP_Success

   tracer_d_module(preformed_alk_ind)%short_name = 'ALK_preformed'
   tracer_d_module(preformed_alk_ind)%long_name  = 'Preformed Alkalinity'
   tracer_d_module(preformed_alk_ind)%units      = 'mmol/m3'
   tracer_d_module(preformed_alk_ind)%tend_units = 'mmol/m3/s'
   tracer_d_module(preformed_alk_ind)%flux_units = 'cm mmol/m3/s'

   tracer_d_module(preformed_po4_ind)%short_name = 'PO4_preformed'
   tracer_d_module(preformed_po4_ind)%long_name  = 'Preformed PO4'
   tracer_d_module(preformed_po4_ind)%units      = 'mmol/m3'
   tracer_d_module(preformed_po4_ind)%tend_units = 'mmol/m3/s'
   tracer_d_module(preformed_po4_ind)%flux_units = 'cm mmol/m3/s'

   tracer_d_module(preformed_o2_ind)%short_name = 'O2_preformed'
   tracer_d_module(preformed_o2_ind)%long_name  = 'Preformed O2'
   tracer_d_module(preformed_o2_ind)%units      = 'mmol/m3'
   tracer_d_module(preformed_o2_ind)%tend_units = 'mmol/m3/s'
   tracer_d_module(preformed_o2_ind)%flux_units = 'cm mmol/m3/s'
!-----------------------------------------------------------------------
!  default namelist settings
!-----------------------------------------------------------------------

   init_preformed_tracers_option = 'unknown'
   init_preformed_tracers_init_file = 'unknown'

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      do while (nml_error > 0)
         read(nml_in, nml=preformed_tracers_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
      call document(subname, 'init_preformed_tracers_option : ' // init_preformed_tracers_option)
      call document(subname, 'init_preformed_tracers_init_file : ' // init_preformed_tracers_init_file)
      call exit_POP(sigAbort, 'stopping in ' /&
                           &/ subname)
   endif

!-----------------------------------------------------------------------
!  broadcast all namelist variables
!-----------------------------------------------------------------------

   call broadcast_scalar(init_preformed_tracers_option , master_task)
   call broadcast_scalar(init_preformed_tracers_init_file , master_task)

!-----------------------------------------------------------------------
!  initialize tracer read
!-----------------------------------------------------------------------

   do n=1, preformed_tracer_cnt
      tracer_inputs(n)%mod_varname  = tracer_d_module(n)%short_name
      tracer_inputs(n)%file_varname = tracer_d_module(n)%short_name
      tracer_inputs(n)%scale_factor = c1
      tracer_inputs(n)%default_val  = c0
      tracer_inputs(n)%filename = init_preformed_tracers_init_file
      tracer_inputs(n)%file_fmt = 'nc'
   end do

!-----------------------------------------------------------------------
!  initialize tracers
!-----------------------------------------------------------------------

   select case (init_preformed_tracers_option)

   case ('restart', 'ccsm_continue', 'ccsm_branch', 'ccsm_hybrid' )

      if (read_restart_filename == 'undefined') then
         call document(subname, 'no restart file to read preformed tracers from')
         call exit_POP(sigAbort, 'stopping in ' /&
                               &/ subname)
         endif

      call rest_read_tracer_block(preformed_tracers_ind_begin, &
                                  init_ts_file_fmt,            &
                                  read_restart_filename,       &
                                  tracer_d_module,             &
                                  TRACER_MODULE)
   case ('ccsm_startup', 'ccsm_startup_spunup', 'zero')

      TRACER_MODULE = c0

   case ('file')

      do n=1, preformed_tracer_cnt
         if (n /= preformed_o2_ind) then
            call file_read_single_tracer(tracer_inputs, TRACER_MODULE, n)
         else
            TRACER_MODULE(:, :, :, n, :, :) = c0
         end if
      end do

      if (n_topo_smooth > 0) then
         do k=1,km
            call fill_points(k, TRACER_MODULE(:, :, k, n, oldtime, :), errorCode)

            if (errorCode /= POP_Success) then
               call POP_ErrorSet(errorCode, &
                    'preformed_tracer_init: error in fill points for tracers(oldtime)')
               return
            endif

            call fill_points(k, TRACER_MODULE(:, :, k, n, curtime, :), errorCode)

            if (errorCode /= POP_Success) then
               call POP_ErrorSet(errorCode, &
                    'preformed_tracer_init: error in fill points for tracers(newtime)')
               return
            endif
         enddo
      endif

   case DEFAULT
     call document(subname, 'unknown init_preformed_tracers_option = ', init_preformed_tracers_option)
     call exit_POP(sigAbort, 'stopping in ' /&
                          &/ subname)
   end select

!-----------------------------------------------------------------------
!  apply land mask to tracers
!-----------------------------------------------------------------------

   do iblock=1,nblocks_clinic
      do n = 1,preformed_tracer_cnt
         do k = 1,km
            where (k > KMT(:,:,iblock))
               TRACER_MODULE(:,:,k,n,curtime,iblock) = c0
               TRACER_MODULE(:,:,k,n,oldtime,iblock) = c0
            end where
         end do
      end do
   enddo

!-----------------------------------------------------------------------
!EOC

 end subroutine preformed_tracer_init

!***********************************************************************
!***********************************************************************
!BOP
! !IROUTINE: preformed_tracer_reset
! !INTERFACE:

 subroutine preformed_tracer_reset(TRACER_RESET, preformed_reset_to_ind, preformed_tracers_ind_begin, bid)

! !DESCRIPTION:
!  reset surface value for preformed tracers: PO4 and ALK
!
! !REVISION HISTORY:
!  same as module

! !USES:

   use time_management, only: mix_pass, c2dtt
   use prognostic, only: PSURF, newtime
   use constants, only: grav
   use grid, only: dz

! !INPUT PARAMETERS:

   integer(int_kind), dimension(preformed_tracer_cnt), intent(in) :: preformed_reset_to_ind
   integer(int_kind), intent(in) :: preformed_tracers_ind_begin
   integer(int_kind), intent(in) :: bid

! !INPUT/OUTPUT PARAMETERS:

   real(r8), dimension(nx_block,ny_block,km,nt), intent(inout) :: &
      TRACER_RESET      ! preformed tracers: ALK_preformed, PO4_preformed,  O2_preformed

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer :: tracer_ind, preformed_tracer_ind

!-----------------------------------------------------------------------

   do tracer_ind = 1, preformed_tracer_cnt

      preformed_tracer_ind = preformed_tracers_ind_begin + (tracer_ind - 1)
      TRACER_RESET(:,:,1,preformed_tracer_ind) = TRACER_RESET(:,:,1,preformed_reset_to_ind(tracer_ind))

   end do

!-----------------------------------------------------------------------
!EOC

 end subroutine preformed_tracer_reset

!***********************************************************************

end module preformed_tracers_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
