module ocpl_roms_mod

!=========================================================================================
!=========================================================================================

!BOP
! !MODULE: roms_comp_mct
! 
! !INTERFACE:

! !DESCRIPTION:
!     This is the interface between ocpl and roms wrt exporting/importing 
!     3d data for 3d coupling of an embedded regional model
!
! !REVISION HISTORY:
!      2019-Sep: first version by Brian Kauffman
!
! !USES:

   use shr_kind_mod, only: IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
   use ocpl_fields_mod  ! new ocpl aVect fields
   use ocpl_data_mod    ! new ocpl aVect data declarations

! ROMS internal data types
   use mod_grid,     only: GRID
   use mod_param,    only: globalISize0 => Lm,              &
                           globalJSize0 => Mm,              &
                           ROMS_levels  => N,               &
                           Ngrids,                          &
                           NtileI,NtileJ,                   &
                           BOUNDS
   use mod_parallel, only: MyRank

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
   private                              ! By default make data private
   SAVE                                 ! save everything

!
! !PUBLIC MEMBER FUNCTIONS:

   public :: ocpl_roms_init
!  public :: ocpl_roms_export
   public :: ocpl_roms_import


! ! PUBLIC DATA:
!
! !REVISION HISTORY:
!    2019 Sep    - B. Kauffman, first version 
!
!EOP
! !PRIVATE MODULE VARIABLES

   integer(IN),parameter  :: nestID = 1    ! roms nest (grid/domain) #1 
   logical      :: iHaveSouth    ! local roms tile includes south boundary
   logical      :: iHaveEast     ! local roms tile includes east  boundary
   logical      :: iHaveNorth    ! local roms tile includes north boundary
   logical      :: iHaveWest     ! local roms tile includes west  boundary
   integer(IN)  ::  localSize,  localISize,  localJSize   ! local tile size
   integer(IN)  :: globalSize, globalISize, globalJSize   ! global nest size

   integer(IN) :: dbug = 1    ! debug level (higher is more


!=========================================================================================
contains
!=========================================================================================

!=========================================================================================
!BOP !====================================================================================
!
! !IROUTINE: ocpl_roms_init  
!
! !DESCRIPTION:
!    init data sent to roms
!
! !INTERFACE:
!
! !REVISION HISTORY:
!    2019 Sep -- Brian Kauffman, initial version
!
! !INTERFACE: ----------------------------------------------------------------------------

subroutine ocpl_roms_init()

! !INPUT/OUTPUT PARAMETERS:

      USE mod_boundary ! internal roms data structures
!EOP
!BOC

!-----  local variables --------------------------------------------------------

#ifdef _OPENMP
   integer, external :: omp_get_max_threads  ! max number of threads that can execute
#endif

   integer(IN)  :: lSize,m,k

   character(*),parameter :: subName =   "(ocpl_roms_init) "
   character(*),parameter :: F03     = '("(ocpl_roms_init) ",a,2es11.3)'
   character(*),parameter :: F09     = '("(ocpl_roms_init) ",a," ",60("="))'

!-------------------------------------------------------------------------------
!  initialize roms curtain BC data types 
!-------------------------------------------------------------------------------

   write(o_logunit,F09) "Enter" ;  call shr_sys_flush(o_logunit)

   ! get local sizes from ROMS parameters
   localISize = BOUNDS(nestID)%IendR(MyRank) - BOUNDS(nestID)%IstrR(MyRank) + 1
   localJSize = BOUNDS(nestID)%JendR(MyRank) - BOUNDS(nestID)%JstrR(MyRank) + 1
   localSize  = localISize * localJSize

   ! get global sizes from ROMS parameters
   globalISize = globalISize0(nestID) + 2
   globalJSize = globalJSize0(nestID) + 2
   globalSize  = globalISize * globalJSize

   nlev_r = ROMS_levels(nestID)

   mpicom_r = mpicom_o ! roms has same mpi comm as ocn comp
   OCNID_r  = OCNID_o  ! roms has same comp ID  as ocn comp

   write(o_logunit,*) subName,"myRank   = ",MyRank
   write(o_logunit,*) subName,"nlev_r   = ",nlev_r
   write(o_logunit,*) subName,"mpicom_r = ",mpicom_r
   write(o_logunit,*) subName,"OCNID_r  = ",OCNID_r 
   write(o_logunit,*) subName,"global ni,nj,n = ",globalISize,globalJSize,globalSize 
   write(o_logunit,*) subName,"local  ni,nj,n = ",localISize , localJSize, localSize
   write(o_logunit,*) subName,"nlev_r         = ",nlev_r 

   !--------------------------------------------------------------------------------------
   write(o_logunit,*) subName,"init curtain gsMaps" ; call shr_sys_flush(o_logunit)
   !--------------------------------------------------------------------------------------
   call ocpl_roms_gsMapInit(mpicom_r, OCNID_r)

   !--------------------------------------------------------------------------------------
   write(o_logunit,*) subName,"init curtain aVects" ; call shr_sys_flush(o_logunit)
   !--------------------------------------------------------------------------------------

   ! allocate pop->roms curtain attribute vectors arrays
   allocate( p2x_2d_rc(4       ) )
   allocate( p2x_3d_rc(4,nlev_r) )

   do k = 1, 4   ! for each of four curtains N,E,S,W
      lsize = 0
      lSize = mct_gsMap_lsize(gsMap_rc(k), mpicom_r) ! local size wrt to curtain gsMap

      !----- init 2d fields specifically for pop/roms coupling -----
      call mct_aVect_init(p2x_2d_rc(k), rList=ocpl_fields_p2x_2d_fields,lsize=lsize)
      call mct_aVect_zero(p2x_2d_rc(k))

      !----- init 3d fields specifically for pop/roms coupling -----
      do m = 1, nlev_r
         call mct_aVect_init(p2x_3d_rc(k,m), rList=ocpl_fields_p2x_3d_fields,lsize=lsize)
         call mct_aVect_zero(p2x_3d_rc(k,m))
      enddo
   enddo

   !--- set aVect field indicies -- this is already done in ocpl_pop_mod.F90
!  k_p2x_2d_So_ssh  = mct_aVect_indexRA(p2x_2d_r   ,"So_ssh" )
!  k_p2x_2d_So_ubar = mct_aVect_indexRA(p2x_2d_r   ,"So_ubar")
!  k_p2x_2d_So_vbar = mct_aVect_indexRA(p2x_2d_r   ,"So_vbar")
!  k_p2x_3d_So_temp = mct_aVect_indexRA(p2x_3d_r(1),"So_temp")
!  k_p2x_3d_So_salt = mct_aVect_indexRA(p2x_3d_r(1),"So_salt")
!  k_p2x_3d_So_uvel = mct_aVect_indexRA(p2x_3d_r(1),"So_uvel")
!  k_p2x_3d_So_vvel = mct_aVect_indexRA(p2x_3d_r(1),"So_vvel")

   !--- DEBUG ---
   if (dbug > 0) then
      m = k_Scurtain 
      k = k_p2x_2d_So_ssh
      write(o_logunit,F03) "<DEBUG> check south curtain..."
      write(o_logunit,F03) "p2x_2d_rc    ssh  min,max: ",minval(p2x_2d_rc(m)  %rAttr(k,:)),maxval(p2x_2d_rc(m)  %rAttr(k,:))
      k = k_p2x_3d_So_temp
      write(o_logunit,F03) "p2x_3d_rc(1) temp min,max= ",minval(p2x_3d_rc(m,1)%rAttr(k,:)),maxval(p2x_3d_rc(m,1)%rAttr(k,:))
      write(o_logunit,F03) "p2x_3d_rc(2) temp min,max= ",minval(p2x_3d_rc(m,2)%rAttr(k,:)),maxval(p2x_3d_rc(m,2)%rAttr(k,:))
      call shr_sys_flush(o_logunit)
   end if

   !--------------------------------------------------------------------------------------
   write(o_logunit,*) subName,"init curtain data arrays" ; call shr_sys_flush(o_logunit)
   !--------------------------------------------------------------------------------------

   BOUNDARY_OCPL(nestID) % zeta_west  = 1.0e30
   BOUNDARY_OCPL(nestID) % zeta_east  = 1.0e30
   BOUNDARY_OCPL(nestID) % zeta_south = 1.0e30
   BOUNDARY_OCPL(nestID) % zeta_north = 1.0e30
   BOUNDARY_OCPL(nestID) % ubar_west  = 1.0e30
   BOUNDARY_OCPL(nestID) % ubar_east  = 1.0e30
   BOUNDARY_OCPL(nestID) % ubar_south = 1.0e30
   BOUNDARY_OCPL(nestID) % ubar_north = 1.0e30
   BOUNDARY_OCPL(nestID) % vbar_west  = 1.0e30
   BOUNDARY_OCPL(nestID) % vbar_east  = 1.0e30
   BOUNDARY_OCPL(nestID) % vbar_south = 1.0e30
   BOUNDARY_OCPL(nestID) % vbar_north = 1.0e30
   BOUNDARY_OCPL(nestID) %    u_west  = 1.0e30
   BOUNDARY_OCPL(nestID) %    u_east  = 1.0e30
   BOUNDARY_OCPL(nestID) %    u_south = 1.0e30
   BOUNDARY_OCPL(nestID) %    u_north = 1.0e30
   BOUNDARY_OCPL(nestID) %    v_west  = 1.0e30
   BOUNDARY_OCPL(nestID) %    v_east  = 1.0e30
   BOUNDARY_OCPL(nestID) %    v_south = 1.0e30
   BOUNDARY_OCPL(nestID) %    v_north = 1.0e30
   BOUNDARY_OCPL(nestID) % temp_west  = 1.0e30
   BOUNDARY_OCPL(nestID) % temp_east  = 1.0e30
   BOUNDARY_OCPL(nestID) % temp_south = 1.0e30 
   BOUNDARY_OCPL(nestID) % temp_north = 1.0e30
   BOUNDARY_OCPL(nestID) % salt_west  = 1.0e30
   BOUNDARY_OCPL(nestID) % salt_east  = 1.0e30
   BOUNDARY_OCPL(nestID) % salt_south = 1.0e30
   BOUNDARY_OCPL(nestID) % salt_north = 1.0e30
   BOUNDARY_OCPL(nestID) % newdata    = .false.

   write(o_logunit,*) subname,"Exit" ;  call shr_sys_flush(o_logunit)

end subroutine ocpl_roms_init

!=========================================================================================
!BOP !====================================================================================
!
! !IROUTINE: ocpl_roms_import
!
! !DESCRIPTION:
!    put lateral BCs into roms
!
! !INTERFACE:
!
! !REVISION HISTORY:
!    2019 Sep -- Brian Kauffman, initial version
!
! !INTERFACE: ----------------------------------------------------------------------------

subroutine ocpl_roms_import()

! !INPUT/OUTPUT PARAMETERS:

      USE mod_boundary ! internal roms data structures
!EOP
!BOC

!-----  local variables --------------------------------------------------------

!  use domain_size, only: km  ! # vertical levels, for 3d coupling

#ifdef _OPENMP
   integer, external :: omp_get_max_threads  ! max number of threads that can execute
#endif
   integer           :: ng = 1               ! roms grid/domain #1 

   character(*), parameter :: subName = "(ocpl_roms_import) "

!-------------------------------------------------------------------------------
!  load roms BC data into roms data types 
!  note: roms uses BC temperature units of Celcius (not Kelvin)
!-------------------------------------------------------------------------------

   write(o_logunit,*) subname,"Enter" ;  call shr_sys_flush(o_logunit)

   BOUNDARY_OCPL(nestID) % newdata    = .false.
   if ( Scurtain) then
      BOUNDARY_OCPL(nestID) % zeta_south = 1.0e30
      BOUNDARY_OCPL(nestID) % ubar_south = 1.0e30
      BOUNDARY_OCPL(nestID) % vbar_south = 1.0e30
      BOUNDARY_OCPL(nestID) %    u_south = 1.0e30
      BOUNDARY_OCPL(nestID) %    v_south = 1.0e30
      BOUNDARY_OCPL(nestID) % temp_south = 5.0
      BOUNDARY_OCPL(nestID) % salt_south = 1.0e30
      BOUNDARY_OCPL(nestID) % newdata    = .true.
   end if
   if ( Ecurtain) then
      BOUNDARY_OCPL(nestID) % zeta_east  = 1.0e30
      BOUNDARY_OCPL(nestID) % ubar_east  = 1.0e30
      BOUNDARY_OCPL(nestID) % vbar_east  = 1.0e30
      BOUNDARY_OCPL(nestID) %    u_east  = 1.0e30
      BOUNDARY_OCPL(nestID) %    v_east  = 1.0e30
      BOUNDARY_OCPL(nestID) % temp_east  = 1.0e30
      BOUNDARY_OCPL(nestID) % salt_east  = 1.0e30
      BOUNDARY_OCPL(nestID) % newdata    = .true.
   end if
   if ( Ncurtain) then
      BOUNDARY_OCPL(nestID) % zeta_north = 1.0e30
      BOUNDARY_OCPL(nestID) % ubar_north = 1.0e30
      BOUNDARY_OCPL(nestID) % vbar_north = 1.0e30
      BOUNDARY_OCPL(nestID) %    u_north = 1.0e30
      BOUNDARY_OCPL(nestID) %    v_north = 1.0e30
      BOUNDARY_OCPL(nestID) % temp_north = 1.0e30
      BOUNDARY_OCPL(nestID) % salt_north = 1.0e30
      BOUNDARY_OCPL(nestID) % newdata    = .true.
   end if
   if ( Wcurtain) then
      BOUNDARY_OCPL(nestID) % zeta_west  = 1.0e30
      BOUNDARY_OCPL(nestID) % ubar_west  = 1.0e30
      BOUNDARY_OCPL(nestID) % vbar_west  = 1.0e30
      BOUNDARY_OCPL(nestID) %    u_west  = 1.0e30
      BOUNDARY_OCPL(nestID) %    v_west  = 1.0e30
      BOUNDARY_OCPL(nestID) % temp_west  = 1.0e30
      BOUNDARY_OCPL(nestID) % salt_west  = 1.0e30
      BOUNDARY_OCPL(nestID) % newdata    = .true.
   end if

   write(o_logunit,*) subname,"Exit" ;  call shr_sys_flush(o_logunit)

end subroutine ocpl_roms_import

!=========================================================================================
!=========================================================================================
!BOP
!
! !IROUTINE: ocpl_roms_gsMapInit
!
! !INTERFACE:
   subroutine ocpl_roms_gsMapInit(comm, compid)
!
! !DESCRIPTION:
!
! !USES:
! use cpl_comm_mod, only: OCN_ID => cpl_comm_mph_cid_ocn
! use OCN_fields_mod
! use NRCM_data_mod, only: compid => compid_p 
! use NRCM_fields_mod
  use shr_sys_mod

! !ARGUMENTS:
   implicit none 
   integer(IN),        intent(in) :: comm
   integer(IN),        intent(in) :: compid
 
! !REVISION HISTORY:
 
!EOP
 
! !LOCAL VARIABLES:
   integer(IN)              :: i, j, ij, k
   integer(IN), allocatable :: indx(:)

   character(*),parameter :: subName = "(ocpl_roms_gsMapInit) "

!-------------------------------------------------------------------------------
! create gsMaps for the four roms BC curtains
! surface forcing aVect gsMaps are created by roms component
!-------------------------------------------------------------------------------

   write(o_logunit,*) subname,"Enter" ;  call shr_sys_flush(o_logunit)

   allocate(gsMap_rc(4)) ! allocate for four BC curtains

   !----------------------------------------------------------------------------
   ! initialize boundary gsMaps for ROMS curtain forcing
   !----------------------------------------------------------------------------
   iHaveWest = .false.
   k = k_Wcurtain
   allocate(indx(localJSize))
   if (BOUNDS(nestID)%IstrR(MyRank).le.1) then
      iHaveWest = .true.
      ij    = 0
      do j  = BOUNDS(nestID)%JstrR(MyRank)+1, BOUNDS(nestID)%JendR(MyRank)+1
         ij = ij + 1
         indx(ij) = j
      enddo
      call mct_gsMap_init(gsMap_rc(k), indx, comm, compid, localJSize, globalJSize)
   else
      call mct_gsMap_init(gsMap_rc(k), indx, comm, compid,        0  , globalJSize)
   endif
   write(o_logunit,*) subName,"west : size global, local = ", &
            mct_gsMap_gsize(gsMap_rc(k)),mct_gsMap_lsize(gsMap_rc(k),comm)
   deallocate(indx)
   write(*,*) subName,"DEBUG west  : lsize = ",mct_gsMap_lsize(gsMap_rc(k),comm)

   !----------------------------------------------------------------------------
   ! NORTH next - create index of local cells
   !----------------------------------------------------------------------------
   iHaveNorth = .false.
   k = k_Ncurtain
   allocate(indx(localISize))
   if (BOUNDS(nestID)%JendR(MyRank)+1.eq.globalJSize) then
      iHaveNorth = .true.
      ij    = 0
      do i  = BOUNDS(nestID)%IstrR(MyRank)+1, BOUNDS(nestID)%IendR(MyRank)+1
         ij = ij + 1
         indx(ij) = i
      enddo
      call mct_gsMap_init(gsMap_rc(k), indx, comm, compid, localISize, globalISize)
   else
      call mct_gsMap_init(gsMap_rc(k), indx, comm, compid,        0  , globalISize)
   endif
   write(o_logunit,*) subName,"north: size global, local = ", &
            mct_gsMap_gsize(gsMap_rc(k)),mct_gsMap_lsize(gsMap_rc(k),comm)
   deallocate(indx)
   write(*,*) subName,"DEBUG north : lsize = ",mct_gsMap_lsize(gsMap_rc(k),comm)

   !----------------------------------------------------------------------------
   ! EAST next - create index of local cells
   !----------------------------------------------------------------------------
   iHaveEast = .false.
   k = k_Ecurtain
   allocate(indx(localJSize))
   if (BOUNDS(nestID)%IendR(MyRank)+1.eq.globalISize) then
      iHaveEast = .true.
      ij    = 0
      do j  = BOUNDS(nestID)%JstrR(MyRank)+1, BOUNDS(nestID)%JendR(MyRank)+1
         ij = ij + 1
         indx(ij) = j
      enddo
      call mct_gsMap_init(gsMap_rc(k), indx, comm, compid, localJSize, globalJSize)
   else
      call mct_gsMap_init(gsMap_rc(k), indx, comm, compid,        0  , globalJSize)
   endif
   write(o_logunit,*) subName,"east : size global, local = ", &
            mct_gsMap_gsize(gsMap_rc(k)),mct_gsMap_lsize(gsMap_rc(k),comm)
   deallocate(indx)
   write(*,*) subName,"DEBUG east  : lsize = ",mct_gsMap_lsize(gsMap_rc(k),comm)

   !----------------------------------------------------------------------------
   ! SOUTH last - create index of local cells
   !----------------------------------------------------------------------------
   iHaveSouth = .false.
   k = k_Scurtain
   allocate(indx(localISize))
   if (BOUNDS(nestID)%JstrR(MyRank).le.1) then
      iHaveSouth = .true.
      ij    = 0
      do i  = BOUNDS(nestID)%IstrR(MyRank)+1, BOUNDS(nestID)%IendR(MyRank)+1
         ij       = ij + 1
         indx(ij) = i
      enddo
      call mct_gsMap_init(gsMap_rc(k), indx, comm, compid, localISize, globalISize)
   else
      call mct_gsMap_init(gsMap_rc(k), indx, comm, compid,        0  , globalISize)
   endif
   write(o_logunit,*) subName,"south: size global, local = ", &
            mct_gsMap_gsize(gsMap_rc(k)),mct_gsMap_lsize(gsMap_rc(k),comm)
   deallocate(indx)
   write(*,*) subName,"DEBUG south : lsize = ",mct_gsMap_lsize(gsMap_rc(k),comm)

   write(o_logunit,*) subname,"Exit" ;  call shr_sys_flush(o_logunit)

end subroutine ocpl_roms_gsMapInit

!=========================================================================================
!=========================================================================================
end module ocpl_roms_mod
