module ocpl_data_mod

!================================================================================
! MODULE: ocpl_data
!
! DESCRIPTION:
!    This module holds data for the ocpl ocean driver 
!
! !REVISION HISTORY:
!    2018-Oct: first version by Brian Kauffman 
!================================================================================

   use mct_mod
   use shr_kind_mod, only: IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
   use shr_mct_mod
   use seq_infodata_mod
   use seq_cdata_mod

   implicit none
   public                 ! default is data public, that's the point of this module

!  PUBLIC FUNCTIONS:

   ! none

!  PUBLIC DATA:

!  integer(IN)     :: logunit = 6 ! use stdout in ocn_comp_mct.F90  ?

   type(mct_aVect) :: x2o_r   ! roms  input: x2o_o mapped to r grid 
   type(mct_aVect) :: r2x_r   ! roms ouptut 
   type(mct_aVect) :: r2x_o   ! roms ouptut: r2x_r mapped to o grid

   character(CL)   :: start_type
   character(CS)   :: case_name

   !--- "o" data (global domain coupler knows about, which is pop domain) ----
!  type(seq_cdata)        , pointer :: cdata_o   ! coupler instantiates this
   type(mct_gGrid)        , pointer :: dom_o 
   type(mct_gsMap)        , pointer :: gsMap_o 
   type(seq_infodata_type), pointer :: infodata_o
   integer(IN)                      :: mpicom_o, OCNID_o
   character(16)                    :: name_o

   !--- "r" data (regional domain coupler doesn't know about, which is roms domain) ----
   type(seq_cdata)                  :: cdata_r   ! regional
   type(mct_gGrid)        , target  :: dom_r    
   type(mct_gsMap)        , target  :: gsMap_r  
!  type(seq_infodata_type), target  :: infodata_r   roms shares infodata_o with pop
   integer(IN)                      :: mpicom_r, OCNID_r
   character(16)                    :: name_r

   type(mct_sMatp)        :: sMatp_r2o
   character(*),parameter :: r2o_mapfile = "/glade/p/cesm/cseg/mapping/makemaps/gom3_to_gx1v7_181014/map_gom3_to_gx1v7_bilin_181014.nc "
   character(*),parameter :: r2o_maptype = "Y"
   type(mct_sMatp)        :: sMatp_o2r
   character(*),parameter :: o2r_mapfile = "/glade/p/cesm/cseg/mapping/makemaps/gom3_to_gx1v7_190809/map_gx1v7_to_gom3_bilinex_190809.nc "
   character(*),parameter :: o2r_maptype = "Y"

#ifdef CPP_VECTOR
   logical :: usevector = .true.
#else
   logical :: usevector = .false.
#endif



   save ! save everything.

!================================================================================
end module ocpl_data_mod
