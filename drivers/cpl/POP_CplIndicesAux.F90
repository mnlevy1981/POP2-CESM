module POP_CplIndicesAux
  
  use seq_flds_mod
  use mct_mod
  use POP_KindsMod

  implicit none

  SAVE
  public                               ! By default make data private

  type :: indices_x2o                  ! coupler indices
    character (POP_charLength) :: name ! eg, 'PFioi_taux1'
    integer :: index                   ! set by calling mct_avect_indexra
  end type
    

!------------------------------------------------------------------------------
!  drv -> ocn  MCOG
!------------------------------------------------------------------------------


  type (indices_x2o),dimension(:),allocatable ::  &
     indices_x2o_ifrac,  & ! fractional ice,                  categories/columns 0:max
     indices_x2o_swpen     ! sw: net penetrating ice,         categories/columns 1:max

contains

  subroutine POP_CplIndicesSetAux( )

    type(mct_aVect) :: o2x      ! temporary
    type(mct_aVect) :: x2o      ! temporary

!-----------------------------------------------------------------------
!   local variables
!-----------------------------------------------------------------------

    integer ncol
    integer ncols_MCOG
    logical lmcog
    
    ! create temporary attribute vectors
    call mct_aVect_init(x2o, rList=seq_flds_x2o_fields, lsize=1)
    call mct_aVect_init(o2x, rList=seq_flds_o2x_fields, lsize=1)

!-----------------------------------------------------------------------
!   MCOG 
!-----------------------------------------------------------------------
!####################### DEBUG ##########################
!   kludge...
    ncols_MCOG = _NCOL_MCOG 
!########################################################

!-----------------------------------------------------------------------
!   define MCOG attribute vector indices
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!   ifrac indices
!-----------------------------------------------------------------------
    allocate (indices_x2o_ifrac(0:ncols_MCOG))
    do ncol = 0, ncols_MCOG
      write(indices_x2o_ifrac(ncol)%name, '(a,i2.2)') 'Si_ifrac_',ncol
      indices_x2o_ifrac(ncol)%index = mct_avect_indexra(x2o,trim(indices_x2o_ifrac(ncol)%name))
    enddo

!-----------------------------------------------------------------------
!   swpen indices
!-----------------------------------------------------------------------
    allocate (indices_x2o_swpen(0:ncols_MCOG))
    do ncol = 0, ncols_MCOG
      if (ncol == 0) then
        indices_x2o_swpen(ncol)%name = 'Foxx_swpen0'
      else
        write(indices_x2o_swpen(ncol)%name, '(a,i2.2)') 'PFioi_swpen',ncol
      endif
      indices_x2o_swpen(ncol)%index = mct_avect_indexra(x2o,trim(indices_x2o_swpen(ncol)%name))
    enddo

  end subroutine POP_CplIndicesSetAux

end module POP_CplIndicesAux
