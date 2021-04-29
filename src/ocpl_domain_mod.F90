module ocpl_domain_mod

!=========================================================================================
!=========================================================================================

!BOP
! !MODULE: ocpl_domain_mod
! 
! !INTERFACE:

! !DESCRIPTION:
!     Setup the ocpl/ocn domain (with which ocpl exchanges data with the CESM coupler)
!
! !REVISION HISTORY:
!      2020-Mar: first version by Brian Kauffman
!
! !USES:

   use shr_kind_mod, only: IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
   use ocpl_data_mod    ! new ocpl aVect data declarations

!  use esmf_mod
   use seq_infodata_mod
   use shr_sys_mod
   use shr_mct_mod
   use mct_mod

   implicit none
   private              ! all data private by default 
   SAVE                 ! save everything

! !PUBLIC MEMBER FUNCTIONS:

   public :: ocpl_domain_init

! !PUBLIC DATA:

! !PRIVATE MODULE DATA:

   integer(IN),parameter :: debug  = 1    ! debug level

!=========================================================================================
contains
!=========================================================================================
!BOP !====================================================================================
!
! !IROUTINE: ocpl_domain_init
!
! !DESCRIPTION:
!     Setup the ocpl/ocn domain with which ocpl exchanges data with the CESM coupler
!     1) input domain specification data
!     2) create an mct GSmap / decompoition for the domain
!     3) init mct gGrid and fill with domain/coordinate data
!
! !INTERFACE:
!
! !REVISION HISTORY:
!    2020 Mar -- Brian Kauffman, initial version
!
! !INTERFACE: ----------------------------------------------------------------------------

subroutine ocpl_domain_init()

   use shr_dmodel_mod, only: shr_dmodel_readgrid

! !INPUT/OUTPUT PARAMETERS:

   implicit none

!  type(mct_gsMap), intent(in)    :: gsMap_o
!  type(mct_ggrid), intent(inout) :: dom_o

!EOP
!BOC

!-----  local variables --------------------------------------------------------

!  use domain_size, only: km  ! # vertical levels, for 3d coupling

   integer(IN) :: k ! curtain index
   integer(IN) :: lsize, gsize, ni_p, nj_p
   integer(IN) :: rank,ierr
   real   (R8), pointer :: data(:)
   integer(IN), pointer :: gIndex(:)

   integer(IN)             :: nk_o  ! depth, assumed to be 1 
   character(*), parameter :: decomp   = "2d1d"
!  HACK/TO-DO: fileName, at least, should be an input namelist param -- BK
!  character(*), parameter :: fileName = "/glade/p/cesm/cseg/inputdata/share/domains/domain.ocn.gom09_gx1v7.200610.nc"
   character(*), parameter :: fileName = "/glade/p/cesm/cseg/inputdata/share/domains/domain.ocn.gst03_tx0.1v3.210108.nc"
   character(*), parameter ::  latName = "yc"
   character(*), parameter ::  lonName = "xc"
   character(*), parameter :: maskName = "mask"
   character(*), parameter :: areaName = "area"
   character(*), parameter :: fracName = "frac"
   logical     , parameter :: readFrac = .false. ! assume (frac = 1.0) iff (mask != 0)

   character(*), parameter :: subName = "(ocpl_domain_init) "

!-----------------------------------------------------------------------------------------
!  
!-----------------------------------------------------------------------------------------

   write(o_logUnit,'(2a)') subname,"Enter" ;  call shr_sys_flush(o_logUnit)

   !--------------------------------------------------------------------------------------
   write(o_logUnit,*) subname,"call read grid spec file, init gGrid_o & gsMap_o"
   !--------------------------------------------------------------------------------------
   call shr_dmodel_readgrid(gGrid_o,gsMap_o, ni_o,nj_o,nk_o, fileName, OCNID_o, mpicom_o, &
        decomp  = decomp  , readFrac = .false. , &
        lonName = lonName , latName  = latName , maskName = maskName , areaName = areaName )

   call seq_infodata_PutData( infodata_o, ocn_nx = ni_o , ocn_ny = nj_o)
   !--------------------------------------------------------------------------------------
   write(o_logUnit,'(2a)') subname,"Exit" ;  call shr_sys_flush(o_logUnit)
   !--------------------------------------------------------------------------------------

end subroutine ocpl_domain_init


!=========================================================================================
!=========================================================================================
end module ocpl_domain_mod
