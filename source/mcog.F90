!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module mcog

!BOP
! !MODULE: mcog
!
! !REVISION HISTORY:
! SVN:$Id$

! !USES

   use kinds_mod
   use domain_size
   use domain
   use blocks
   use io
   use io_types
   use constants
   use exit_mod
   use grid
   use communicate
   use global_reductions
   use broadcast
   use tavg
   use time_management
   use mcog_forcing_fields

   implicit none
   private
   save

!-----------------------------------------------------------------------
! Overview comments for the MCOG CPT (Climate Process Team) project:
!
!    Ocean Mixing Processes Associated with High Spatial Heterogeneity
!           in Sea Ice and the Implications for Climate Models
!
! or otherwise known as MCOG (Multi-Column Ocean Grid) from the term
! used in the proposal description. All code changes for this project
! are denoted by "! MCOG +" and "! MCOG -" comments (for starting and
! ending of modified code).
!
! Here we give a commented overview of the changes made to implement
! MCOG for the ocean component. See the sea ice component and coupler
! for descriptions of fields sent from the sea ice component and ice
! fraction weighted averaging over the coupling interval, as well as
! averaging over the open ocean fraction for the appropriate fluxes.
! (Search for "Overview" in the modified code sub-directories to
! find descriptive comments).
!
! Here we start with the category ice fractions and various ice/ocn
! fluxes and stresses received by the ocean component.
!
! The category fields received by the ocean component from the coupler
! in MCOG are:
!
!            a_n              n_th category sea ice concentration
!            FW_n             n_th category fresh water flux
!            FQ_n             n_th category heat flux
!            FSW_n            n_th category shortwave flux
!            taux_n           n_th category ice/ocn zonal stress
!            tauy_n           n_th category ice/ocn meridional stress
!
! for {n=1,Ncat} categories. All fields are averaged over the coupling
! interval weighted by the time varying ice category fractions, and the
! open ocean (category n=0) fields of shortwave fluxes and ice/ocn
! stresses are computed and sent to the ocean.
!
! Please note that in the sea ice and coupler descriptions of the ice
! fraction weighting for fluxes, equations are given showing the full
! normalized averages. But actually, the sea ice component sends category
! fluxes to the coupler which are already multiplied by category ice
! fraction. The normalization is completed in the ocean component where
! the coupling interval averaged category ice fractions are divided into
! the ice fraction weighted fluxes received.
!
! Specifically, here are the averages:
!                         Nstp
!      a_n = [ a_nm ] = Sum ( a_nm ) / Nstp                            (1)
!                         m=1
!
!                         Nstp                         Ncat
!      A   = [ A_m ] = Sum ( A_m ) / Nstp      A_m = Sum a_nm          (1a)
!                         m=1                          n=1
!                                    Nstp
!      F_n = [ a_nm F_nm ] / a_n = Sum ( a_nm F_nm ) / (a_n Nstp)      (2)
!                                    m=1
!                                     Nstp
!  F_0 = [(1-A_m) Fatm_m] / (1-A) = Sum ((1-A_m) Fatm_m) / (1-A)Nstp   (3)
!                                     m=1
!
! As just said, the sea ice component multiplies each category and
! time step flux (F_nm) by the corresponding category and time step
! sea ice concentration (a_nm). The coupler then performs the usual
! averaging over Nstp, and sends the resulting fluxes to the ocean
! component (i.e. makes averages as in Eq. 1). For the open ocean category
! fluxes which are formed in the coupler, the coupler multiplies by the
! open ocean fraction (i.e. term 1-A_m in Eq.3) before averaging. The ocean
! component then completes the average by dividing with a_n from Eq 1 in
! Eq 2, and by (1-A) as in Eq 3 for category 0 fluxes, where A is the
! total ice fraction given by summing over the coupling interval A_m, which
! is the mth time step total ice fraction, as given in Eq 1a above.
!
! In the ocean, the KPP boundary layer parameterization requires as inputs
! surface forcing and profiles of temperature, salinity and tracers in the
! full ocean column. The temperature, salinity and tracer profiles are for
! the full ocean grid box. Specific surface forcings are the surface friction
! velocity (computed from surface stress), solar and non-solar buoyancy fluxes
! (evaluated from the surface shortwave flux and the sensible/latent heat
! fluxes plus the longwave flux, and finally the fresh water flux), and the
! kinematic surface tracer fluxes for both heat and virtual salt flux. The
! output of KPP are the full column diffusivity and viscosity, which are then
! input to the vertical diffusion solver to update the temperature, salinity
! and tracer profiles.
!
! We note some specifics about the category kinematic surface tracer fluxes
! for heat and virtual salt flux. Contributions to the fluxes for heat
! include atmosphere to ocean fluxes (sensible, latent, up and down longwave)
! which are segregated into the open ocean category (0), snow and ice runoff
! which is included in both open ocean and all cateogries (since there is
! no category specific partition that makes physical sense here), and the
! usual sea ice to ocean fluxes placed in the appropriate category. The
! virtual salt flux surface tracer flux segregates the atmosphere precipitation
! and ocean evaporation into the open ocean category forcing, with land
! and ice runoff placed again in every category, with the sea ice to ocean
! category melt fluxes placed in the appropriate category forcing.
!
! Finally, the coupler multiplies the atmosphere to ocean heat and water
! fluxes by the open ocean fraction, which must be divided back out when
! forming category specific (like open ocean) forcing fluxes.
!
! For the open ocean shortwave absorbed, we note that the ocean model
! typically uses a diurnal cycle partition of the received daily shortwave
! in its time integration. We take the open ocean shortwave absorbed from
! the coupler and partition it diurnally into the open ocean shortwave
! category forcing used by the ocean model. The remainder of the sea ice
! to ocean penetrating shortwave fluxes are placed in the appropriate
! category shortwave forcing.
!
! To apply MCOG to KPP, and minimize overhead of always doing all Ncat+1
! columns, we will run KPP over only those ice-covered columns for which
! {a_n > 0}. Surface forcing is defined using the category fluxes as just
! explained. These will include surface friction velocity u*_n, solar buoyancy
! flux BS_n, non-solar buoyancy flux BNS_n, and kinematic surface tracer fluxes.
! Then, KPP will be run up to N+1 times (depending on ice concentration, which
! sets 0 <= N <= Ncat) producing the multi-column diffusivities k_n and
! viscosities mu_n for {n=0,N}.
!
! We proceed by homogenizing k_n and mu_n for the grid-cell as:
!                      N                         N
!               k = Sum k_n a_n          mu = Sum mu_n a_n                  (4)
!                     n=0                       n=0
! and then run the vertical diffusion solver once to produce modified temperature,
! salinity and tracer profiles. Note that other forcing arrays used outside
! of KPP are also evaluated column by column and aggregrated as in Eq.4, in
! particular the boundary layer depth.
!
! In the implementation here, we have placed a category loop in KPP which runs
! from 0 (i.e. open ocean) through all five (currently, Ncat = 5) sea ice
! categories, and then finish with the originally forced single column ocean grid
! (termed SCOG). When MCOG is prognostic, the code uses the MCOG aggregated terms
! as in Eq 1 (and others as noted), while if it is diagnostic the original SCOG
! terms are used. In this way, we can compare diagnostically MCOG and SCOG from
! history field information regardless of whether MCOG is prognostic or not. We
! save both MCOG and SCOG diffusivities and viscosities (Eq. 1) to history file.
! We also save the individual category terms as well.
!
!   Bruce P. Briegleb   September 2011
!
!-----------------------------------------------------------------------

! !PUBLIC MEMBER FUNCTIONS:
 
   public :: init_mcog,        &
             init_mcog_ecosys, &
             tavg_mcog

! !PUBLIC DATA MEMBERS:

   logical (log_kind), public ::  &
      lmcog,                      &! namelist variable; if true, mcog is on  
      lmcog_debug                  ! namelist variable; if true, print mcog debugging stmts

   integer, parameter, public ::  &
      max_bins_MCOG = _NCOL_MCOG   !###### KLUDGE

   integer (int_kind), dimension(0:max_bins_MCOG), public ::  & 
      MCOG_OCN_bins

   integer, public ::  &
      ncols_MCOG,      &
      nbins_MCOG 

   integer (int_kind), dimension(:), allocatable, public ::  &
      tavg_IFRAC_MCOG,    &! tavg ids for sea ice fraction categories (binned)
      tavg_SHF_QSW_MCOG,  &! tavg ids for ocn SHF_QSW   categories 0-ncols_MCOG (binned)
      tavg_SWPEN_MCOG      ! tavg ids for sea ice penetrating shortwave categories 0-ncols_MCOG (binned)

   integer (int_kind), public ::  &
      tavg_DIFRAC,        &! tavg id for difference aggregate sea ice fraction and 1
      tavg_DSWPEN          ! tavg id for difference aggregate category swpen and original

   integer (int_kind), public ::  &
      tavg_NBINS_IFRAC,   &! inter-model consistency check for IFRAC
      tavg_NBINS_SWPEN     ! inter-model consistency check for SWPEN

   real (r8), parameter, public :: &
      minicefrac = 1.0E-13_r8

!EOP
!BOC

!EOC
!***********************************************************************

   contains

!***********************************************************************
!BOP
! !IROUTINE: init_mcog
! !INTERFACE:

 subroutine init_mcog

! !DESCRIPTION:
!  Initializes the logical on/off switch for the multi-column ocean grid
!  (mcog) representation in the vertical ocean mixing.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
!     input namelist variables (for other public namelist variables, see above)
!
!-----------------------------------------------------------------------

   namelist /mcog_nml/ lmcog, lmcog_debug, ncols_MCOG, MCOG_OCN_bins

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

   character (char_len) :: &  
      string,              &! for defining history fields
      id_string             ! id string

   integer (int_kind) ::   &
      nbin,                &! bin index
      ncol,                &! column index
      nml_error,           &! namelist i/o error flag
      nu                    ! i/o unit number

   type (block) ::         &
      this_block            ! block information for current block

!-----------------------------------------------------------------------
!
!     set defaults for niw parameters, then read them from namelist
!
!-----------------------------------------------------------------------

   lmcog       = .false.
   lmcog_debug = .false.
   ncols_MCOG  = 0
   do nbin = 0, max_bins_MCOG
    MCOG_OCN_bins(nbin) = nbin
   enddo
 
!-----------------------------------------------------------------------
!
!  read namelist input and broadcast variables
!
!-----------------------------------------------------------------------

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      do while (nml_error > 0)
         read(nml_in, nml=mcog_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar (nml_error, master_task)
   if (nml_error /= 0) then
     call exit_POP (SigAbort, 'ERROR reading mcog_nml')
   endif

   if (my_task == master_task) then
      write(stdout,blank_fmt)
      write(stdout,ndelim_fmt)
      write(stdout,blank_fmt)
      write(stdout,*) ' MCOG information'
      write(stdout,blank_fmt)
      write(stdout,*) ' mcog_nml namelist settings:'
      write(stdout,blank_fmt)
      write(stdout,mcog_nml)
      write(stdout,blank_fmt)
      write(stdout,*) ' lmcog                   = ',  lmcog
      write(stdout,*) ' ncols_MCOG              = ',  ncols_MCOG
      write(stdout,*) ' MCOG_OCN_bins           = ',  MCOG_OCN_bins
      call POP_IOUnitsFlush(stdout)
   endif

   call broadcast_scalar (lmcog,         master_task)
   call broadcast_scalar (ncols_MCOG,    master_task)
   call broadcast_array  (MCOG_OCN_bins, master_task)

!-----------------------------------------------------------------------
!
!  exit if mcog is not enabled
!
!-----------------------------------------------------------------------

    if (.not. lmcog) return   


!-----------------------------------------------------------------------
!  account for "binning" of ice categories
!-----------------------------------------------------------------------

    !*** determine number of independent bins. Easier for now to 
    !    require increasing bins only, but that is really not 
    !    necessary.  BUT, it makes counting bins easier
    do ncol = 1, ncols_MCOG-1
      if (MCOG_OCN_bins(ncol+1) < MCOG_OCN_bins(ncol)) then
       call document ('init_mcog', 'ncol     ', ncol)
       call document ('init_mcog', 'bin number     ', MCOG_OCN_bins(ncol))
       call document ('init_mcog', 'FATAL ERROR: bin number not increasing')
       call exit_POP (SigAbort, 'FATAL ERROR: bin number not increasing')
      endif
    enddo

    nbins_MCOG = 1  
    do ncol = 0, ncols_MCOG-1
      if (MCOG_OCN_bins(ncol+1) > MCOG_OCN_bins(ncol)) then
         nbins_MCOG = nbins_MCOG + 1
      endif
    enddo

    call document ('init_mcog', 'number of bins', nbins_MCOG)


!-----------------------------------------------------------------------
!
!  allocate and initialize MCOG arrays
!  allocate and define time-averaged MCOG field metadata ("tavg fields")
!
!-----------------------------------------------------------------------

!------------------
!  IFRAC_MCOG
!------------------
   allocate (IFRAC_MCOG_ALL(nx_block,ny_block,0:ncols_MCOG,max_blocks_clinic))
   allocate (IFRAC_MCOG    (nx_block,ny_block,0:nbins_MCOG-1,max_blocks_clinic))
   allocate (DIFRAC        (nx_block,ny_block,max_blocks_clinic))
   IFRAC_MCOG_ALL  = c0
   IFRAC_MCOG      = c0
   DIFRAC          = c0

   allocate (tavg_IFRAC_MCOG(0:nbins_MCOG-1))
   do nbin = 0,nbins_MCOG-1
     if (nbin == 0) then
       string = 'open water fraction'
     else
       write(string,'(a,i2.2)') 'sea ice fraction bin number ',nbin
     endif
     write(id_string,'(a,i2.2)') 'IFRAC_',nbin
  
     call define_tavg_field(tavg_IFRAC_MCOG(nbin),trim(id_string),2, &
                            long_name=trim(string),                  &
                            units='fraction',                        &
                            grid_loc='2110',                         &
                            coordinates  ='TLONG TLAT  time'         )
   enddo ! nbin

   string = 'sea ice fraction difference with 1'
   call define_tavg_field(tavg_DIFRAC,'DIFRAC',2,             &
                          long_name=trim(string),             &
                          units='fraction',                   &
                          grid_loc='2110',                    &
                          coordinates  ='TLONG TLAT  time'    )

!------------------
!  SWPEN_MCOG
!------------------
   allocate (SWPEN_MCOG(nx_block,ny_block,0:nbins_MCOG-1,max_blocks_clinic))
   allocate (DSWPEN    (nx_block,ny_block,max_blocks_clinic))
   SWPEN_MCOG = c0
   DSWPEN     = c0        

   allocate (tavg_SWPEN_MCOG(0:nbins_MCOG-1))
   do nbin = 0,nbins_MCOG-1
     write(string,   '(a,i2.2)') 'sea ice penetrating shortwave  bin/group number ',nbin
     write(id_string,'(a,i2.2)') 'SWPEN_',nbin
  
     call define_tavg_field(tavg_SWPEN_MCOG(nbin),trim(id_string),2, &
                            long_name=trim(string),                  &
                            units='W m-2',                           &
                            grid_loc='2110',                         &
                            coordinates  ='TLONG TLAT  time'         )
   enddo ! nbin

   string = 'sea ice penetrating shortwave aggregate diff with original'
   call define_tavg_field(tavg_DSWPEN,'DSWPEN',2,             &
                          long_name=trim(string),             &
                          units='W m-2',                      &
                          grid_loc='2110',                    &
                          coordinates  ='TLONG TLAT  time'    )


   string = 'SUM_NBINS_IFRAC'
   call define_tavg_field(tavg_NBINS_IFRAC,trim(string) ,2,   &
                          long_name=trim(string),             &
                          units='W m-2  ',                    &
                          grid_loc='2110',                    &
                          coordinates  ='TLONG TLAT  time'    )
   string = 'SUM_NBINS_SWPEN'
   call define_tavg_field(tavg_NBINS_SWPEN,trim(string) ,2,   &
                          long_name=trim(string),             &
                          units='W m-2  ',                    &
                          grid_loc='2110',                    &
                          coordinates  ='TLONG TLAT  time'    )

!------------------
!  SHF_QSW_MCOG
!------------------
   allocate (SHF_QSW_MCOG     (nx_block,ny_block,0:nbins_MCOG-1,max_blocks_clinic))
   allocate (SHF_QSW_SAVE_MCOG(nx_block,ny_block,0:nbins_MCOG-1,max_blocks_clinic))
   allocate (SHF_QSW_RAW_MCOG (nx_block,ny_block,0:nbins_MCOG-1,max_blocks_clinic))

   SHF_QSW_MCOG      = c0
   SHF_QSW_SAVE_MCOG = c0
   SHF_QSW_RAW_MCOG  = c0

   if (lmcog_debug)  call document('pop_set_coupled_forcing','initialize SHF_QSW_MCOG and SHF_QSW_SAVE_MCOG')

   allocate (tavg_SHF_QSW_MCOG(0:nbins_MCOG-1))

   do nbin = 0,nbins_MCOG-1
     write(string,   '(a,i2.2)') 'SHF_QSW  bin/group ',nbin
     write(id_string,'(a,i2.2)') 'SHF_QSW_',nbin
  
     call define_tavg_field(tavg_SHF_QSW_MCOG(nbin),trim(id_string),2, &
                            long_name=trim(string),                    &
                            units='W m-2  ',                           &
                            grid_loc='2110',                           &
                            coordinates  ='TLONG TLAT  time'           )
   enddo ! nbin


!-----------------------------------------------------------------------
!EOC

 end subroutine init_mcog

 subroutine init_mcog_ecosys

! !DESCRIPTION:
!  mcog ecosys initialization
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  exit if mcog is not enabled
!
!-----------------------------------------------------------------------

    if (.not. lmcog) return   

!-----------------------------------------------------------------------
!
!  allocate and initialize MCOG arrays
!  allocate and define time-averaged MCOG field metadata ("tavg fields")
!
!-----------------------------------------------------------------------

!---------------------
!  ecosystem MCOG vars
!---------------------
   allocate (PAR_out_MCOG(nx_block,ny_block,0:nbins_MCOG-1,max_blocks_clinic))
   PAR_out_MCOG = c0

   allocate (PAR_in_MCOG    (nx_block,ny_block,0:nbins_MCOG-1,max_blocks_clinic))
   allocate (PAR_avg_MCOG   (nx_block,ny_block,0:nbins_MCOG-1,max_blocks_clinic))
   allocate (NITRIF_MCOG    (nx_block,ny_block,0:nbins_MCOG-1,max_blocks_clinic))
   allocate (light_lim_MCOG (nx_block,ny_block,0:nbins_MCOG-1,max_blocks_clinic))
   PAR_in_MCOG    = c0
   PAR_avg_MCOG   = c0
   NITRIF_MCOG    = c0
   light_lim_MCOG = c0

!-----------------------------------------------------------------------

 end subroutine init_mcog_ecosys

!***********************************************************************
!BOP
! !IROUTINE: tavg_mcog
! !INTERFACE:

 subroutine tavg_mcog

! !DESCRIPTION:
!   calls tavg accumulation routines for all MCOG variables
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

   character (char_len) :: &  
      string,              &! for defining history fields
      id_string             ! id string

   integer (int_kind) ::   &
      iblock,              &
      nbin,                &! bin index
      ncol                  ! column index

   logical (log_kind), dimension(nx_block,ny_block) :: &
      LTEST_IFRAC

   real (r8), dimension(nx_block,ny_block) :: &
      WORK                  ! temporary storage

   type (block) ::         &
      this_block            ! block information for current block


!-----------------------------------------------------------------------
!  accumulate mcog fields by aggregating over ice fraction
!-----------------------------------------------------------------------

   if (lmcog) then 
      
      do iblock = 1, nblocks_clinic
         do nbin = 0, nbins_MCOG-1
            
            if (nbin == 0) then ! assumes open-ice is bin0 
               LTEST_IFRAC = IFRAC_MCOG(:,:,nbin,iblock) < c1
            else
               LTEST_IFRAC = IFRAC_MCOG(:,:,nbin,iblock) > c0
            endif
            
            if (accumulate_tavg_now(tavg_SHF_QSW_MCOG(nbin))) then
               WORK = c0
               where( LTEST_IFRAC ) WORK(:,:) = SHF_QSW_MCOG(:,:,nbin,iblock) / hflux_factor
               call accumulate_tavg_field(WORK,tavg_SHF_QSW_MCOG(nbin),iblock,1)
            endif

            call accumulate_tavg_field(IFRAC_MCOG(:,:,nbin,iblock),tavg_IFRAC_MCOG(nbin),iblock,1)
            call accumulate_tavg_field(SWPEN_MCOG(:,:,nbin,iblock),tavg_SWPEN_MCOG(nbin),iblock,1)

         enddo ! nbin loop

         call accumulate_tavg_field(DIFRAC(:,:,iblock),tavg_DIFRAC,iblock,1)
         call accumulate_tavg_field(DSWPEN(:,:,iblock),tavg_DSWPEN,iblock,1)

         if (accumulate_tavg_now(tavg_NBINS_IFRAC)) then
            WORK(:,:) = c0
            do nbin = 1,nbins_MCOG-1
               WORK(:,:) = WORK(:,:) + IFRAC_MCOG(:,:,nbin,iblock)
            enddo ! nbin
            call accumulate_tavg_field(WORK(:,:),tavg_NBINS_IFRAC,iblock,1)
         endif

         if (accumulate_tavg_now(tavg_NBINS_SWPEN)) then
            WORK(:,:) = c0
            do nbin = 1,nbins_MCOG-1
               WORK(:,:) = WORK(:,:) + SWPEN_MCOG(:,:,nbin,iblock)*IFRAC_MCOG(:,:,nbin,iblock)
            enddo ! nbin
            call accumulate_tavg_field(WORK(:,:),tavg_NBINS_SWPEN,iblock,1)
         endif

      enddo ! iblock loop
   endif  ! lmcog 

 end subroutine tavg_mcog
end module mcog

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
