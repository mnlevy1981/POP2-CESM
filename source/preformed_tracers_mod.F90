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
      preformed_tracer_cnt = 2

!-----------------------------------------------------------------------
!  relative tracer indices
!-----------------------------------------------------------------------

   integer(int_kind), parameter :: &
      preformed_alk_ind = 1     ! preformed alk index

   integer(int_kind), parameter :: &
      preformed_po4_ind = 2     ! preformed po4 index

!-----------------------------------------------------------------------
!  derived type & parameter for tracer index lookup
!-----------------------------------------------------------------------

   type(ind_name_pair), dimension(preformed_tracer_cnt) :: &
      ind_name_table = (/ ind_name_pair(preformed_alk_ind, 'ALK_preformed'), ind_name_pair(preformed_po4_ind, 'PO4_preformed')  /)

!EOC
!*****************************************************************************

contains

!*****************************************************************************
!BOP
! !IROUTINE: preformed_tracer_init
! !INTERFACE:

 subroutine preformed_tracer_init(tracer_d_module, TRACER_MODULE, errorCode)

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

! !INPUT/OUTPUT PARAMETERS:

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

   logical(log_kind) :: &
      lnml_found             ! Was preformed_tracers_nml found ?    

   integer(int_kind) :: &
      n,                   & ! index for looping over tracers
      k,                   & ! index for looping over depth levels
      iblock,              & ! index for looping over blocks
      nml_error              ! namelist i/o error flag

!     l,                   & ! index for looping over time levels

   character (char_len) ::  &
      preformed_tracer_restart_filename  ! modified file name for restart file


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


!-----------------------------------------------------------------------
!  initialize tracers
!-----------------------------------------------------------------------

   TRACER_MODULE = c0
   ! TODO: read from a restart if necessary! (put init_ts_option back in)

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

 subroutine preformed_tracer_reset(TRACER_OLD, preformed_reset_to_ind, TRACER_RESET, bid)

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

   integer(int_kind), intent(in) :: bid
   integer(int_kind), dimension(preformed_tracer_cnt), intent(in) :: preformed_reset_to_ind
   real(r8), dimension(nx_block,ny_block,km,nt), intent(in) :: &          !for getting the ALK and PO4 concentrations at the surface
      TRACER_OLD

! !INPUT/OUTPUT PARAMETERS:

   real(r8), dimension(nx_block,ny_block,km,preformed_tracer_cnt), intent(inout) :: &
      TRACER_RESET      ! preformed tracers: ALK_preformed, PO4_preformed

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer :: tracer_ind

!-----------------------------------------------------------------------

   do tracer_ind = 1, preformed_tracer_cnt                                !resetting preformed tracers to total PO4/ALK at the surface

      TRACER_RESET(:,:,1,tracer_ind) = TRACER_OLD(:,:,1,preformed_reset_to_ind(tracer_ind))

   end do

!-----------------------------------------------------------------------
!EOC

 end subroutine preformed_tracer_reset

!***********************************************************************

end module preformed_tracers_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
