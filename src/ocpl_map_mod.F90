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


!  use esmf_mod
   use seq_flds_mod
   use seq_cdata_mod
   use seq_infodata_mod
   use seq_timemgr_mod
   use shr_file_mod 
   use shr_cal_mod, only : shr_cal_date2ymd
   use shr_sys_mod
   use shr_mct_mod
   use mct_mod

   implicit none
!  private              ! all data private by default 
   SAVE                 ! save everything

! !PUBLIC MEMBER FUNCTIONS:

   public :: ocpl_map_init
   public :: ocpl_map_pop2roms
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

!  integer         :: ni_p   ! number of longitudes on pop grid
!  integer         :: nj_p   ! number of latitudes  on pop grid
!  integer         :: ni_r   ! number of longitudes on roms grid
!  integer         :: nj_r   ! number of latitudes  on roms grid

   integer,parameter :: debug = 1   ! debug level

   !--- pop -> roms curtain maps ---

   type(mct_sMatp) :: sMatp_p2rc(4)  ! 4-curtains: S,E,N,W

!=========================================================================================
contains
!=========================================================================================
!BOP !====================================================================================
!
! !IROUTINE: ocpl_map_init
!
! !DESCRIPTION:
!    initialize maps between pop & roms
!
! !INTERFACE:
!
! !REVISION HISTORY:
!    2019 Oct -- Brian Kauffman, initial version
!
! !INTERFACE: ----------------------------------------------------------------------------

subroutine ocpl_map_init()

! !INPUT/OUTPUT PARAMETERS:

!EOP
!BOC

!-----  local variables --------------------------------------------------------

!  use domain_size, only: km  ! # vertical levels, for 3d coupling

   integer(IN) :: k ! curtain index
!  integer(IN) :: ni_p, nj_p, ni_rc, nj_rc  ! output of shr_mct_sMatPInitnc
   character(CL) :: mapfile,maptype
   character(*),parameter :: configFile = "ocpl_maps.rc"

   character(*), parameter :: subName = "(ocpl_map_init) "

!-----------------------------------------------------------------------------------------
!  
!-----------------------------------------------------------------------------------------

   write(o_logunit,*) subname,"Enter" ;  call shr_sys_flush(o_logunit)

   !--------------------------------------------------------------------------------------
   write(o_logunit,*) subname,"Init pop->roms curtain maps"
   !--------------------------------------------------------------------------------------

   k = k_Scurtain
   call shr_mct_queryConfigFile(mpicom_o,configFile, &
        "pop2roms_Scurtain_file:",mapfile,"pop2roms_Scurtain_type:",maptype)
   write(o_logunit,'(4a)') subName, "file = ",trim(mapfile)
   call shr_mct_sMatPInitnc(sMatp_p2rc(k),gsMap_o,gsMap_rc(k),trim(mapfile),trim(maptype),mpicom_o)

   k = k_Ecurtain
   call shr_mct_queryConfigFile(mpicom_o,configFile, &
        "pop2roms_Ecurtain_file:",mapfile,"pop2roms_Ecurtain_type:",maptype)
   write(o_logunit,'(4a)') subName, "file = ",trim(mapfile)
   call shr_mct_sMatPInitnc(sMatp_p2rc(k),gsMap_o,gsMap_rc(k),trim(mapfile),trim(maptype),mpicom_o)

   k = k_Ncurtain
   call shr_mct_queryConfigFile(mpicom_o,configFile, &
        "pop2roms_Ncurtain_file:",mapfile,"pop2roms_Ncurtain_type:",maptype)
   write(o_logunit,'(4a)') subName, "file = ",trim(mapfile)
   call shr_mct_sMatPInitnc(sMatp_p2rc(k),gsMap_o,gsMap_rc(k),trim(mapfile),trim(maptype),mpicom_o)

   if (Wcurtain) then
      k = k_Wcurtain
      call shr_mct_queryConfigFile(mpicom_o,configFile, &
           "pop2roms_Wcurtain_file:",mapfile,"pop2roms_Wcurtain_type:",maptype)
      write(o_logunit,'(4a)') subName, "file = ",trim(mapfile)
      call shr_mct_sMatPInitnc(sMatp_p2rc(k),gsMap_o,gsMap_rc(k),trim(mapfile),trim(maptype),mpicom_o)
   else
      write(*,*) subname,"DEBUG Wcurtain = false" ;  call shr_sys_flush(o_logunit)
   end if

write(*,*) subname,"DEBUG exit" ;  call shr_sys_flush(o_logunit)

   write(o_logunit,*) subname,"Exit" ;  call shr_sys_flush(o_logunit)

end subroutine ocpl_map_init

!=========================================================================================
!BOP !====================================================================================
!
! !IROUTINE: ocpl_map_pop2roms
!
! !DESCRIPTION:
!     map data from pop to roms 
!
! !INTERFACE:
!
! !REVISION HISTORY:
!    2019 Oct -- Brian Kauffman, initial version
!
! !INTERFACE: ----------------------------------------------------------------------------

subroutine ocpl_map_pop2roms()

! !INPUT/OUTPUT PARAMETERS:

!EOP
!BOC

!-----  local variables --------------------------------------------------------

!  use domain_size, only: km  ! # vertical levels, for 3d coupling

   integer(IN) :: k ! curtain index
   integer(IN) :: kfld ! field index
   real(R8)    :: tmin,tmax
   logical     :: do_mapping

   character(*), parameter :: subName = "(ocpl_map_pop2roms) "

!-----------------------------------------------------------------------------------------
!  
!-----------------------------------------------------------------------------------------

   write(o_logunit,*) subname,"Enter" ;  call shr_sys_flush(o_logunit)

   !--------------------------------------------------------------------------------------
   ! debug/sanity-check on data
   !--------------------------------------------------------------------------------------
   if (debug > 0) then
      write(*,'(2a,5i6)') subName,'<DEBUG> lsize: pop, S,N,W,E = ', &
         mct_aVect_lsize(p2x_2d_p), &
         mct_aVect_lsize(p2x_2d_rc(k_Scurtain)), mct_aVect_lsize(p2x_2d_rc(k_Ncurtain)), & 
         mct_aVect_lsize(p2x_2d_rc(k_Wcurtain)), mct_aVect_lsize(p2x_2d_rc(k_Ecurtain))
   endif

   !--------------------------------------------------------------------------------------
   ! map 2d curtain fields: ubar, vbar, ssh
   !--------------------------------------------------------------------------------------
   do k=1,4 ! four curtains: N,E,S,W
      call mct_aVect_zero(p2x_2d_rc(k))

      do_mapping = .false.
      if (k==k_Scurtain .and. Scurtain==.true.) do_mapping = .true.
      if (k==k_Ecurtain .and. Ecurtain==.true.) do_mapping = .true.
      if (k==k_Ncurtain .and. Ncurtain==.true.) do_mapping = .true.
      if (k==k_Wcurtain .and. Wcurtain==.true.) do_mapping = .true.

      !f (do_mapping) call mct_sMat_avMult(p2x_2d_p,sMatp_p2rc(k),p2x_2d_rc(k),rList=ocpl_fields_p2x_2d_fields)
      if (do_mapping) call mct_sMat_avMult(p2x_2d_p,sMatp_p2rc(k),p2x_2d_rc(k))

      if (debug>0 ) then
         tmin = -999.
         tmax = -999.
         if (mct_aVect_lsize(p2x_2d_rc(k)) > 0 ) then
            kfld = mct_aVect_indexRA(p2x_2d_rc(k),"So_ssh" )
            tmin = minval(p2x_3d_p(1)%rAttr(kfld,:))
            tmax = maxval(p2x_3d_p(1)%rAttr(kfld,:))
         end if
         write(*,'(2a,i2,2es11.3)') subname,"<DEBUG> k, ssh min,max = ",k,tmin,tmax 
      end if
   end do

   !--------------------------------------------------------------------------------------
   ! map 3d curtain fields: salt,temperature, u, v
   !--------------------------------------------------------------------------------------


   write(o_logunit,*) subname,"Exit" ;  call shr_sys_flush(o_logunit)

end subroutine ocpl_map_pop2roms

!=========================================================================================
!=========================================================================================
end module ocpl_map_mod
