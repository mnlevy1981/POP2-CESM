!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module mcog_forcing_fields

!BOP
! !MODULE: mcog_forcing_fields

! !DESCRIPTION:
!  Contains the forcing fields necessary for supporting MCOG

! !REVISION HISTORY:
!  SVN:$Id$
!
! !USES:

   use kinds_mod

   implicit none
   save

!EOP
!BOC
! !PUBLIC DATA MEMBERS:

   real (r8), allocatable, dimension(:,:,:,:), public, target ::  &
      IFRAC_MCOG_ALL,    &! open water fraction of grid cell for ALL categories/columns
      IFRAC_MCOG,        &! open water fraction of grid cell for multiple categories/columns (in bins)
      SWPEN_MCOG,        &! sea ice shortwave penetrating flux for multiple categories/columns (in bins)
      SHF_QSW_MCOG,      &! penetrating shortwave for multiple categories/columns (in bins)
      SHF_QSW_SAVE_MCOG, &! penetrating shortwave for all categories from coupler (in bins)
      SHF_QSW_RAW_MCOG    ! SHF_QSW_RAW for all catetories (in bins)

   real (r8), allocatable, dimension(:,:,:), public, target ::  &
      DIFRAC,         &! difference: IFRAC - sum(IFRAC_MCOG)
      DSWPEN           ! difference: SWPEN - sum(SWPEN_MCOG)

   real (r8), allocatable, dimension(:,:,:,:), public, target ::  &
      PAR_out_MCOG     

   real (r8), allocatable, dimension(:,:,:,:), public, target ::  &
      PAR_in_MCOG,    &
      PAR_avg_MCOG,   &
      NITRIF_MCOG,    &
      light_lim_MCOG   


!***********************************************************************

 end module mcog_forcing_fields

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
