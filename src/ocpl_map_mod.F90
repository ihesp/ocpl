module ocpl_map_mod

!=========================================================================================
!=========================================================================================

!BOP
! !MODULE: ocpl_map_mod
! 
! !INTERFACE:

! !DESCRIPTION:
!     This is the interface between ocpl and roms wrt exporting/importing 
!     3d data for 3d coupling of an embedded regional model
!
! !REVISION HISTORY:
!      2019-Oct: first version by Brian Kauffman
!
! !USES:

   use shr_kind_mod, only: IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
   use ocpl_fields_mod  ! new ocpl aVect fields
   use ocpl_data_mod    ! new ocpl aVect data declarations


   use mct_mod
!  use esmf_mod
   use seq_flds_mod
   use seq_cdata_mod
   use seq_infodata_mod
   use seq_timemgr_mod
   use shr_file_mod 
   use shr_cal_mod, only : shr_cal_date2ymd
   use shr_sys_mod

   implicit none
!  private              ! all data private by default 
   SAVE                 ! save everything

!
! !PUBLIC MEMBER FUNCTIONS:

!  public :: ocpl_map_init
!  public :: ocpl_map_pop2roms
!  public :: ocpl_map_roms2pop

! ! PUBLIC DATA:

   !--- horizontal maps (2d) ---
   type(mct_sMatp)        :: sMatp_r2o   ! maps roms -> ocn/pop
   character(*),parameter :: r2o_mapfile = "/glade/p/cesm/cseg/mapping/makemaps/gom3_to_gx1v7_181014/map_gom3_to_gx1v7_bilin_181014.nc "
   character(*),parameter :: r2o_maptype = "Y"
   type(mct_sMatp)        :: sMatp_o2r   ! maps ocn/pop -> roms
   character(*),parameter :: o2r_mapfile = "/glade/p/cesm/cseg/mapping/makemaps/gom3_to_gx1v7_190809/map_gx1v7_to_gom3_bilinex_190809.nc "
   character(*),parameter :: o2r_maptype = "X"

! !PRIVATE MODULE VARIABLES

   integer         :: ni_p   ! number of longitudes on pop grid
   integer         :: nj_p   ! number of latitudes  on pop grid
   integer         :: ni_r   ! number of longitudes on roms grid
   integer         :: nj_r   ! number of latitudes  on roms grid

   !--- pop -> roms curtain maps ---

   type(mct_sMatp)        :: sMatp_p2r_rc(4)  ! 4-curtains: index 1,2,3,4 => S,E,N,W

   type(mct_sMatp)        :: sMatp_p2rS  ! maps to roms South curtain
   character(*),parameter :: o2rs_mapfile = "/glade/p/cesm/cseg/mapping/makegrids/gom3_191017/gom3_Scurtain.191017.nc"
   character(*),parameter :: o2rs_maptype = "Y"
   type(mct_sMatp)        :: sMatp_p2rE  ! maps to roms East  curtain
   character(*),parameter :: o2re_mapfile = "/glade/p/cesm/cseg/mapping/makegrids/gom3_191017/gom3_Ecurtain.191017.nc"
   character(*),parameter :: o2re_maptype = "X"
   type(mct_sMatp)        :: sMatp_p2rN  ! maps to roms North curtain
   character(*),parameter :: o2rn_mapfile = "/glade/p/cesm/cseg/mapping/makegrids/gom3_191017/gom3_Ncurtain.191017.nc"
   character(*),parameter :: o2rn_maptype = "Y"
   type(mct_sMatp)        :: sMatp_p2rW  ! maps to roms West  curtain
   character(*),parameter :: o2rw_mapfile = "/glade/p/cesm/cseg/mapping/makegrids/gom3_191017/gom3_Wcurtain.191017.nc"
   character(*),parameter :: o2rw_maptype = "X"



!=========================================================================================

!=========================================================================================
!=========================================================================================
end module ocpl_map_mod
