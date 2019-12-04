module ocpl_roms_mod

!=========================================================================================
!=========================================================================================

!BOP
! !MODULE: ocpl_roms_mod
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
   use mod_param,    only: globalISize0 => Lm,  &
                           globalJSize0 => Mm,  &
                           ROMS_levels  => N,   &
                           BOUNDS
   use mod_parallel, only: MyRank
   use communicate,  only: master_task,my_task

   use mct_mod
!  use esmf_mod
   use seq_flds_mod
   use seq_cdata_mod
   use seq_infodata_mod
   use seq_timemgr_mod
   use shr_file_mod 
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

   integer(IN),parameter :: nestID = 1 ! roms nest (grid/domain) #1 
   logical     :: iHave_Scurtain       ! local roms tile includes south boundary
   logical     :: iHave_Ecurtain       ! local roms tile includes east  boundary
   logical     :: iHave_Ncurtain       ! local roms tile includes north boundary
   logical     :: iHave_Wcurtain       ! local roms tile includes west  boundary
   integer(IN) ::  localSize,  localISize,  localJSize   ! local tile size
   integer(IN) :: globalSize, globalISize, globalJSize   ! global nest size

   integer(IN) :: debug = 0    ! debug level (higher is more)

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
   if (debug>0) write(o_logunit,'(2a,i3)') subName,"debug level = ",debug

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
   if (debug>0) then
      m = k_Scurtain
      lSize = mct_gsMap_lsize(gsMap_rc(m), mpicom_r) ! local size wrt to Scurtain
      if (lsize > 0) then
         write(o_logUnit,F03) "check south curtain..."
         k = k_p2x_2d_So_ssh
         write(o_logUnit,F03) "   p2x_2d_rc    ssh  min,max: ",minval(p2x_2d_rc(m)  %rAttr(k,:)),maxval(p2x_2d_rc(m)  %rAttr(k,:))
         k = k_p2x_3d_So_temp
         write(o_logUnit,F03) "   p2x_3d_rc(1) temp min,max= ",minval(p2x_3d_rc(m,1)%rAttr(k,:)),maxval(p2x_3d_rc(m,1)%rAttr(k,:))
         write(o_logUnit,F03) "   p2x_3d_rc(9) temp min,max= ",minval(p2x_3d_rc(m,9)%rAttr(k,:)),maxval(p2x_3d_rc(m,9)%rAttr(k,:))
      else
         write(o_logUnit,F03) "check south curtain... lSize=0 => interior tile"
      end if
      call shr_sys_flush(6)
   end if

   !--------------------------------------------------------------------------------------
   write(o_logunit,*) subName,"init curtain data arrays" ; call shr_sys_flush(o_logunit)
   !--------------------------------------------------------------------------------------

   BOUNDARY_OCPL(nestID) % bypass     = .true.  ! tell roms NOT to use this data
   BOUNDARY_OCPL(nestID) % newdata    = .false. ! tell roms this data has NOT been updated
   BOUNDARY_OCPL(nestID) % debug      = debug   ! debug level for write statements
   BOUNDARY_OCPL(nestID) % debug      = 3       ! debug level for write statements

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

   write(o_logunit,'(2a)') subname,"Exit" ;  call shr_sys_flush(o_logunit)

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
   type(mct_aVect),allocatable,save :: roms2D_BC(:)   ! global gather of 2d curtain data
   type(mct_aVect),allocatable,save :: roms3D_BC(:,:) ! global gather of 3d curtain data
   integer(IN) :: i,j,ij                   ! indicies for grid cell
   integer(IN) :: k,n                      ! indicies for curtain and level
   integer(IN) :: stat                     ! return code
   integer(IN) :: kfld                     ! index of some aVect field
   real(R8)    :: tmin,tmax                ! min/max value of field
   logical     :: first_call = .true.      ! flags one-time initializations

   character(*), parameter :: subName = "(ocpl_roms_import) "

!-------------------------------------------------------------------------------
!  load roms BC data into roms data types 
!  note: roms uses BC temperature units of Celcius (not Kelvin)
!-------------------------------------------------------------------------------

   write(o_logunit,'(2a)') subname,"Enter" ;  call shr_sys_flush(o_logunit)
   if (debug>0) write(o_logunit,'(2a,i3)') subName,"debug level = ",debug

  !-----------------------------------------------------------------------------
  ! create global (not decomposed/distributed) curtain aVects
  !-----------------------------------------------------------------------------
  if (first_call) then
     allocate (roms2D_BC(4))
     allocate (roms3D_BC(4,nlev_r))
  end if
  first_call = .false.

  !-----------------------------------------------------------------------------
  ! global gather/broadcast: non-decomposed roms curtain data
  !-----------------------------------------------------------------------------
   do k = 1,4  ! for each curtain: N,E,S,W
      call    mct_aVect_gather(p2x_2d_rc(k), roms2D_BC(k), &
                                gsMap_rc(k), master_task, mpicom_r, stat)
      call    mct_aVect_bcast (roms2D_BC(k), master_task, mpicom_r, stat)

      do n = 1,nlev_r
         call mct_aVect_gather(p2x_3d_rc(k,n), roms3D_BC(k,n), &
                                gsMap_rc(k)  , master_task, mpicom_r, stat)
         call mct_aVect_bcast (roms3D_BC(k,n), master_task, mpicom_r, stat)
      enddo
   enddo

  !--- check values of global gather/broadcast ---
   if (debug>0 ) then
      do k = 1,4  ! for each curtain: N,E,S,W
         kfld = mct_aVect_indexRA(roms2D_BC(k),"So_ssh" )
         tmin = minval( roms2D_BC(k)%rAttr(kfld,:) )
         tmax = maxval( roms2D_BC(k)%rAttr(kfld,:) )
         write(o_logUnit,'(2a,i2,2es11.3)') subname,"global k, ssh   min,max = ",k,tmin,tmax 

         kfld = mct_aVect_indexRA(roms3D_BC(k,1),"So_temp" )
         tmin = minval( roms3D_BC(k,1)%rAttr(kfld,:) )
         tmax = maxval( roms3D_BC(k,1)%rAttr(kfld,:) )
         write(o_logUnit,'(2a,i2,2es11.3)') subname,"global k, T bot min,max = ",k,tmin,tmax 
         tmin = minval( roms3D_BC(k,nlev_r)%rAttr(kfld,:) )
         tmax = maxval( roms3D_BC(k,nlev_r)%rAttr(kfld,:) )
         write(o_logUnit,'(2a,i2,2es11.3)') subname,"global k, T top min,max = ",k,tmin,tmax 
      end do
   end if

  !-----------------------------------------------------------------------------
  ! unit conversion and vector rotation as per roms internal conventions
  !-----------------------------------------------------------------------------
  ! TO DO

  !-----------------------------------------------------------------------------
  ! put gathered data into roms internal data types 
  !-----------------------------------------------------------------------------
   BOUNDARY_OCPL(nestID) % bypass  = .false. ! tell roms     to use this data
   BOUNDARY_OCPL(nestID) % newdata = .false.
   if ( do_Scurtain) then
      k = k_Scurtain
      do i  = 0, globalISize0(nestID) + 1   ! 0-based roms internal data
         ij = i + 1                         ! 1-based aVect 
         BOUNDARY_OCPL(nestID) % zeta_south(i) = roms2D_BC(k)%rAttr(k_p2x_2d_So_ssh ,ij) *0.01_r8
         BOUNDARY_OCPL(nestID) % ubar_south(i) = roms2D_BC(k)%rAttr(k_p2x_2d_So_ubar,ij) *0.01_r8
         BOUNDARY_OCPL(nestID) % vbar_south(i) = roms2D_BC(k)%rAttr(k_p2x_2d_So_vbar,ij) *0.01_r8
         do n = 1,nlev_r
            BOUNDARY_OCPL(nestID) %    u_south(i,n) = roms3D_BC(k,n)%rAttr(k_p2x_3d_So_uvel,ij) *0.01_r8
            BOUNDARY_OCPL(nestID) %    v_south(i,n) = roms3D_BC(k,n)%rAttr(k_p2x_3d_So_vvel,ij) *0.01_r8
            BOUNDARY_OCPL(nestID) % temp_south(i,n) = roms3D_BC(k,n)%rAttr(k_p2x_3d_So_temp,ij)
            BOUNDARY_OCPL(nestID) % salt_south(i,n) = roms3D_BC(k,n)%rAttr(k_p2x_3d_So_salt,ij) *1000.0_r8
         end do
      end do
      BOUNDARY_OCPL(nestID) % newdata    = .true.
   end if
   if ( do_Ecurtain) then
      k = k_Ecurtain
      do j  = 0, globalJSize0(nestID) + 1   ! 0-based roms internal data
         ij = j + 1                         ! 1-based aVect 
         BOUNDARY_OCPL(nestID) % zeta_east (j) = roms2D_BC(k)%rAttr(k_p2x_2d_So_ssh ,ij) *0.01_r8
         BOUNDARY_OCPL(nestID) % ubar_east (j) = roms2D_BC(k)%rAttr(k_p2x_2d_So_ubar,ij) *0.01_r8
         BOUNDARY_OCPL(nestID) % vbar_east (j) = roms2D_BC(k)%rAttr(k_p2x_2d_So_vbar,ij) *0.01_r8
         do n = 1,nlev_r
            BOUNDARY_OCPL(nestID) %    u_east (j,n) = roms3D_BC(k,n)%rAttr(k_p2x_3d_So_uvel,ij) *0.01_r8
            BOUNDARY_OCPL(nestID) %    v_east (j,n) = roms3D_BC(k,n)%rAttr(k_p2x_3d_So_vvel,ij) *0.01_r8
            BOUNDARY_OCPL(nestID) % temp_east (j,n) = roms3D_BC(k,n)%rAttr(k_p2x_3d_So_temp,ij)
            BOUNDARY_OCPL(nestID) % salt_east (j,n) = roms3D_BC(k,n)%rAttr(k_p2x_3d_So_salt,ij) *1000.0_r8
         end do
      end do
      BOUNDARY_OCPL(nestID) % newdata    = .true.
   end if
   if ( do_Ncurtain) then
      k = k_Ncurtain
      do i  = 0, globalISize0(nestID) + 1   ! 0-based roms internal data
         ij = i + 1                         ! 1-based aVect 
         BOUNDARY_OCPL(nestID) % zeta_north(i) = roms2D_BC(k)%rAttr(k_p2x_2d_So_ssh ,ij) *0.01_r8
         BOUNDARY_OCPL(nestID) % ubar_north(i) = roms2D_BC(k)%rAttr(k_p2x_2d_So_ubar,ij) *0.01_r8
         BOUNDARY_OCPL(nestID) % vbar_north(i) = roms2D_BC(k)%rAttr(k_p2x_2d_So_vbar,ij) *0.01_r8
         do n = 1,nlev_r
            BOUNDARY_OCPL(nestID) %    u_north(i,n) = roms3D_BC(k,n)%rAttr(k_p2x_3d_So_uvel,ij) *0.01_r8
            BOUNDARY_OCPL(nestID) %    v_north(i,n) = roms3D_BC(k,n)%rAttr(k_p2x_3d_So_vvel,ij) *0.01_r8
            BOUNDARY_OCPL(nestID) % temp_north(i,n) = roms3D_BC(k,n)%rAttr(k_p2x_3d_So_temp,ij)
            BOUNDARY_OCPL(nestID) % salt_north(i,n) = roms3D_BC(k,n)%rAttr(k_p2x_3d_So_salt,ij) *1000.0_r8
         end do
      end do
      BOUNDARY_OCPL(nestID) % newdata    = .true.
   end if
   if ( do_Wcurtain) then
      do j  = 0, globalJSize0(nestID) + 1   ! 0-based roms internal data
         ij = j + 1                         ! 1-based aVect 
         BOUNDARY_OCPL(nestID) % zeta_west (j) = roms2D_BC(k)%rAttr(k_p2x_2d_So_ssh ,ij) *0.01_r8
         BOUNDARY_OCPL(nestID) % ubar_west (j) = roms2D_BC(k)%rAttr(k_p2x_2d_So_ubar,ij) *0.01_r8
         BOUNDARY_OCPL(nestID) % vbar_west (j) = roms2D_BC(k)%rAttr(k_p2x_2d_So_vbar,ij) *0.01_r8
         do n = 1,nlev_r
            BOUNDARY_OCPL(nestID) %    u_west (j,n) = roms3D_BC(k,n)%rAttr(k_p2x_3d_So_uvel,ij) *0.01_r8
            BOUNDARY_OCPL(nestID) %    v_west (j,n) = roms3D_BC(k,n)%rAttr(k_p2x_3d_So_vvel,ij) *0.01_r8
            BOUNDARY_OCPL(nestID) % temp_west (j,n) = roms3D_BC(k,n)%rAttr(k_p2x_3d_So_temp,ij)
            BOUNDARY_OCPL(nestID) % salt_west (j,n) = roms3D_BC(k,n)%rAttr(k_p2x_3d_So_salt,ij) *1000.0_r8
         end do
      end do
      BOUNDARY_OCPL(nestID) % newdata    = .true.
   end if

   write(o_logunit,'(2a)') subname,"Exit" ;  call shr_sys_flush(o_logunit)

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

   write(o_logunit,'(2a)') subname,"Enter" ;  call shr_sys_flush(o_logunit)
   if (debug>0) write(o_logunit,'(2a,i3)') subName,"debug level = ",debug

   allocate(gsMap_rc(4)) ! allocate for four BC curtains

   !----------------------------------------------------------------------------
   ! initialize boundary gsMaps for ROMS curtain forcing
   !----------------------------------------------------------------------------
   iHave_Wcurtain = .false.
   k = k_Wcurtain
   allocate(indx(localJSize))
   if (BOUNDS(nestID)%IstrR(MyRank).le.1) then
      iHave_Wcurtain = .true.
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

   !----------------------------------------------------------------------------
   ! NORTH next - create index of local cells
   !----------------------------------------------------------------------------
   iHave_Ncurtain = .false.
   k = k_Ncurtain
   allocate(indx(localISize))
   if (BOUNDS(nestID)%JendR(MyRank)+1.eq.globalJSize) then
      iHave_Ncurtain = .true.
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

   !----------------------------------------------------------------------------
   ! EAST next - create index of local cells
   !----------------------------------------------------------------------------
   iHave_Ecurtain = .false.
   k = k_Ecurtain
   allocate(indx(localJSize))
   if (BOUNDS(nestID)%IendR(MyRank)+1.eq.globalISize) then
      iHave_Ecurtain = .true.
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

   !----------------------------------------------------------------------------
   ! SOUTH last - create index of local cells
   !----------------------------------------------------------------------------
   iHave_Scurtain = .false.
   k = k_Scurtain
   allocate(indx(localISize))
   if (BOUNDS(nestID)%JstrR(MyRank).le.1) then
      iHave_Scurtain = .true.
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

   write(o_logunit,'(2a)') subname,"Exit" ;  call shr_sys_flush(o_logunit)

end subroutine ocpl_roms_gsMapInit

!=========================================================================================
!=========================================================================================
end module ocpl_roms_mod
