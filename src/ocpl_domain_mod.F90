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
!      2021-DEC: cleanup, add namelist, use only clause - Jim Edwards
! !USES:

   use shr_kind_mod, only: IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
   use ocpl_data_mod, only: mpicom_o, o_logunit, ggrid_o, gsmap_o, ni_o, nj_o, ocnid_o, infodata_o    ! new ocpl aVect data declarations
   use seq_infodata_mod, only : seq_infodata_getdata, seq_infodata_putdata

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
   use mpi, only : mpi_comm_rank
! !INPUT/OUTPUT PARAMETERS:

   implicit none

!  type(mct_gsMap), intent(in)    :: gsMap_o
!  type(mct_ggrid), intent(inout) :: dom_o

!EOP
!BOC

!-----  local variables --------------------------------------------------------


   integer(IN) :: k ! curtain index
   integer(IN) :: lsize, gsize, ni_p, nj_p
   integer(IN) :: rank,ierr
   real   (R8), pointer :: data(:)
   integer(IN), pointer :: gIndex(:)

   integer(IN)             :: nk_o  ! depth, assumed to be 1
   character(len=8) :: decomp   = "2d1d"
   character(len=256) :: domainfile = "/glade/p/cesm/cseg/inputdata/share/domains/domain.ocn.gst03_tx0.1v3.210514.nc"
   character(len=8) ::  latName = "yc"
   character(len=8) ::  lonName = "xc"
   character(len=8) :: maskName = "mask"
   character(len=8) :: areaName = "area"
   character(len=8) :: fracName = "frac"
   logical      :: readFrac = .false. ! assume (frac = 1.0) iff (mask != 0)
   integer :: nml_in, nml_error
   integer :: my_task
   integer, parameter :: master_task = 0
   character(*), parameter :: nml_filename = 'ocpl_in'
   character(*), parameter :: subName = "(ocpl_domain_init) "
   namelist /ocpl_nml/ domainfile, latName, lonName, maskName, areaName, fracName, readFrac, decomp
!-----------------------------------------------------------------------------------------
!
!-----------------------------------------------------------------------------------------
   call mpi_comm_rank(mpicom_o, my_task, ierr)
   if(my_task == master_task) then
      write(o_logUnit,'(2a)') subname,"Enter"

      !--------------------------------------------------------------------------------------
      write(o_logUnit,*) subname,"call read grid spec file, init gGrid_o & gsMap_o"
      !--------------------------------------------------------------------------------------

      open (newunit=nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      do while (nml_error > 0)
         read(nml_in, nml=ocpl_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call shr_dmodel_readgrid(gGrid_o,gsMap_o, ni_o,nj_o,nk_o, domainfile, OCNID_o, mpicom_o, &
        decomp  = decomp  , readFrac = .false. , &
        lonName = lonName , latName  = latName , maskName = maskName , areaName = areaName )

   call seq_infodata_PutData( infodata_o, ocn_nx = ni_o , ocn_ny = nj_o)
   if(my_task == master_task) then
      !--------------------------------------------------------------------------------------
      write(o_logUnit,'(2a)') subname,"Exit"
      !--------------------------------------------------------------------------------------
   endif
end subroutine ocpl_domain_init


!=========================================================================================
!=========================================================================================
end module ocpl_domain_mod
