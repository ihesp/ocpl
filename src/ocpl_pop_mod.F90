module ocpl_pop_mod

!=========================================================================================
!=========================================================================================

!BOP
! !MODULE: pop_comp_mct
! 
! !INTERFACE:

! !DESCRIPTION:
!     This is the interface between ocpl and pop wrt exporting/importing 
!     3d data for 3d coupling of an embedded regional model
!
! !REVISION HISTORY:
!      2019-Sep: first version by Brian Kauffman
!
! !USES:

!  use shr_kind_mod, only: IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
   use ocpl_fields_mod  ! new ocpl aVect fields
   use ocpl_data_mod    ! new ocpl aVect data declarations

!  use POP_KindsMod
   use POP_ErrorMod
   use POP_CommMod
   use POP_FieldMod
   use POP_GridHorzMod
   use POP_HaloMod
   use POP_IOUnitsMod
   use POP_MCT_vars_mod

   use ocpl_fields_mod

#ifdef NRCMRESTORING
   use forcing_pt_interior, only: PT_INTERIOR_DATA     ! alter interal pop data for 2-way coupling
   use forcing_pt_interior, only: PT_RESTORE_RTAU      ! alter interal pop data for 2-way coupling
   use forcing_pt_interior, only: PT_RESTORE_MAX_LEVEL ! alter interal pop data for 2-way coupling
   use forcing_s_interior , only:  S_INTERIOR_DATA     ! alter interal pop data for 2-way coupling
   use forcing_s_interior , only:  S_RESTORE_RTAU      ! alter interal pop data for 2-way coupling
   use forcing_s_interior , only:  S_RESTORE_MAX_LEVEL ! alter interal pop data for 2-way coupling
#endif

   use mct_mod
!  use esmf_mod
   use seq_flds_mod
   use seq_cdata_mod
   use seq_infodata_mod
   use seq_timemgr_mod
   use shr_file_mod 
   use shr_cal_mod, only : shr_cal_date2ymd
   use shr_sys_mod
   use perf_mod
   use ocn_communicator,  only: mpi_communicator_ocn

   use kinds_mod,         only: int_kind, r8
   use domain_size,       only: km  ! # vertical levels, for 3d coupling
   use POP_CplIndices
   use POP_KindsMod
   use POP_ErrorMod
   use POP_InitMod,       only: POP_Initialize1, POP_Initialize2, &
                                timer_total, cpl_ts 
   use communicate,       only: my_task, master_task
   use constants
   use blocks
   use domain,            only: distrb_clinic, POP_haloClinic
   use exit_mod
   use forcing_shf,       only: SHF_QSW
   use forcing_sfwf,      only: lsend_precip_fact, precip_fact
   use forcing_fields
   use forcing_coupled,   only: ncouple_per_day,  &
                                update_ghost_cells_coupler_fluxes, &
                                rotate_wind_stress, pop_set_coupled_forcing, &
                                pop_init_coupled,  &
                                orb_eccen, orb_obliqr, orb_lambm0, orb_mvelpp
   use ice,               only: tfreez, tmelt, liceform,QFLUX, QICE, AQICE, &
                                tlast_ice
   use grid,              only: TLAT, TLON, KMT
   use global_reductions, only: global_sum_prod
   use io_tools,          only: document
   use named_field_mod,   only: named_field_register, named_field_get_index, &
                                named_field_set, named_field_get
   use prognostic
   use timers,            only: get_timer, timer_start, timer_stop
   use diagnostics,       only: check_KE
   use output,            only: output_driver
   use step_mod,          only: step
   use time_management
   use registry

   implicit none
   private                              ! By default make data private
   SAVE                                 ! save everything

!
! !PUBLIC MEMBER FUNCTIONS:

   public :: ocpl_pop_init
   public :: ocpl_pop_export
!  public :: ocpl_pop_import


! ! PUBLIC DATA:
!
! !REVISION HISTORY:
!    2019 Sep    - B. Kauffman, first version 
!
!EOP
! !PRIVATE MODULE VARIABLES


!  real (r8),allocatable :: SBUFF_SUM(:,:,:,:) ! accum/sum/tavg export quantities
   real (r8) ::  tlast_coupled

   integer (int_kind)  ::   nsend, nrecv

   type(seq_infodata_type), pointer :: infodata   

   integer(IN) :: dbug = 0    ! debug level (higher is more)

!=========================================================================================
contains
!=========================================================================================

!BOP !====================================================================================
!
! !IROUTINE: ocpl_pop_init
!
! !DESCRIPTION:
!    Initialize ocpl/pop 3d interface data
!
! !INTERFACE:
!
! !REVISION HISTORY:
!    2019 Sep -- Brian Kauffman, initial version
!
! !INTERFACE: ----------------------------------------------------------------------------

subroutine ocpl_pop_init( p2x_p, p2x_2d_p, p2x_3d_p)

! !INPUT/OUTPUT PARAMETERS:

   type(mct_aVect)             , intent(in   ) :: p2x_p
   type(mct_aVect)             , intent(inout) :: p2x_2d_p
   type(mct_aVect),pointer     , intent(inout) :: p2x_3d_p(:)

!  use domain_size, only: km  ! # vertical levels, for 3d coupling

   integer(int_kind) :: lsize              ! size of local aVect
   integer(int_kind) :: k                  ! array index
#ifdef _OPENMP
   integer, external :: omp_get_max_threads  ! max number of threads that can execute
#endif

   character(*), parameter :: subName = "(ocpl_pop_init) "

!------------------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------------------

   write(stdout,*) subname,"Enter"

   nlev_p = km
   lsize = mct_aVect_lsize( p2x_p )

   !----- init 2d fields specifically for pop/roms coupling -----
   call mct_aVect_init(p2x_2d_p, rList=ocpl_fields_p2x_2d_fields,lsize=lsize)
   call mct_aVect_zero(p2x_2d_p)

   !----- init 3d fields specifically for pop/roms coupling -----
   write(stdout,*) subName,"     nlev_p    = ",nlev_p
   allocate(p2x_3d_p(nlev_p))

   do k = 1, nlev_p
      call mct_aVect_init(p2x_3d_p(k), rList=ocpl_fields_p2x_3d_fields,lsize=lsize)
      call mct_aVect_zero(p2x_3d_p(k))
   enddo

   !--- set aVect field indicies ---
   p2x_2d_So_ssh  = mct_aVect_indexRA(p2x_2d_p   ,"So_ssh" )
   p2x_2d_So_ubar = mct_aVect_indexRA(p2x_2d_p   ,"So_ubar")
   p2x_2d_So_vbar = mct_aVect_indexRA(p2x_2d_p   ,"So_vbar")
   p2x_3d_So_temp = mct_aVect_indexRA(p2x_3d_p(1),"So_temp")
   p2x_3d_So_salt = mct_aVect_indexRA(p2x_3d_p(1),"So_salt")
   p2x_3d_So_uvel = mct_aVect_indexRA(p2x_3d_p(1),"So_uvel")
   p2x_3d_So_vvel = mct_aVect_indexRA(p2x_3d_p(1),"So_vvel")

   !--- DEBUG ---
   if (dbug > 0) then
      k = p2x_2d_So_ssh
      write(stdout,*) subname,"p2x_2d_p    ssh  min,max: ",minval(p2x_2d_p   %rAttr(k,:)),maxval(p2x_2d_p   %rAttr(k,:))
      k = p2x_3d_So_temp
      write(stdout,*) subname,"p2x_3d_p(1) temp min,max= ",minval(p2x_3d_p(1)%rAttr(k,:)),maxval(p2x_3d_p(1)%rAttr(k,:))
      write(stdout,*) subname,"p2x_3d_p(2) temp min,max= ",minval(p2x_3d_p(2)%rAttr(k,:)),maxval(p2x_3d_p(1)%rAttr(k,:))
      call flushm (stdout)
   end if

   write(stdout,*) subname,"Exit"

end subroutine ocpl_pop_init

!=========================================================================================
!BOP !====================================================================================
!
! !IROUTINE: ocpl_pop_export
!
! !DESCRIPTION:
!    get output data from pop, ocpl/pop 3d interface  
!
! !INTERFACE:
!
! !REVISION HISTORY:
!    2019 Sep -- Brian Kauffman, initial version
!
! !INTERFACE: ----------------------------------------------------------------------------

subroutine ocpl_pop_export( p2x_2d_p, p2x_3d_p)

! !INPUT/OUTPUT PARAMETERS:

   type(mct_aVect)             , intent(inout) :: p2x_2d_p
   type(mct_aVect),pointer     , intent(inout) :: p2x_3d_p(:)

!EOP
!BOC

!-----------------------------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------------------------

!  use domain_size, only: km  ! # vertical levels, for 3d coupling

#ifdef _OPENMP
   integer, external :: omp_get_max_threads  ! max number of threads that can execute
#endif

   integer (int_kind) ::  &
      k,i,j,n,kuse

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) ::  &
      WORK                ! local work arrays

   integer (int_kind) :: &
      iblock              ! block index

   type (block) ::       &
      this_block          ! block information for current block

   real (r8), dimension(nx_block,ny_block) ::   &
      WORK1_nrcm,       & ! local 2d work space for NRCM
      WORK2_nrcm          ! local 2d work space for NRCM

   real (r8), dimension(nx_block,ny_block,km) ::   &
      WORK3_nrcm,       & ! local 3d work space for NRCM
      WORK4_nrcm          ! local 3d work space for NRCM

   character(*), parameter :: subName = "(ocpl_pop_export) "

!------------------------------------------------------------------------------
!  extract data from POP and put into corresponding attribute vectors
!------------------------------------------------------------------------------

   write(stdout,*) subname,"Enter" ; call flushm (stdout)

      ! accumulate variables
      n = 0
      do iblock = 1,nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)

         ! move UPTROP and VBTROP to tgrid
         call ugrid_to_tgrid(WORK1_nrcm(:,:),UBTROP(:,:,curtime,iblock),iblock)
         call ugrid_to_tgrid(WORK2_nrcm(:,:),VBTROP(:,:,curtime,iblock),iblock)

         ! move UVEL and VVEL to tgrid
         do k = 1,km
            call ugrid_to_tgrid(WORK3_nrcm(:,:,k),UVEL(:,:,k,curtime,iblock),iblock)
            call ugrid_to_tgrid(WORK4_nrcm(:,:,k),VVEL(:,:,k,curtime,iblock),iblock)
         enddo

         do j = this_block%jb,this_block%je
         do i = this_block%ib,this_block%ie
            n = n + 1

            p2x_2d_p%rAttr(p2x_2d_So_ssh ,n) = p2x_2d_p%rAttr(p2x_2d_So_ssh,n) &
                                             +      PSURF(i,j,curtime,iblock)/grav
            p2x_2d_p%rAttr(p2x_2d_So_ubar,n) = p2x_2d_p%rAttr(p2x_2d_So_ubar,n) &
                                             +      WORK1_nrcm(i,j)
            p2x_2d_p%rAttr(p2x_2d_So_vbar,n) = p2x_2d_p%rAttr(p2x_2d_So_vbar,n) &
                                             +      WORK2_nrcm(i,j)
            do k = 1,km
               if (KMT(i,j,iblock).gt.0) then
                  kuse = min(k, KMT(i,j,iblock))
                  p2x_3d_p(k)%rAttr(p2x_3d_So_temp,n) = p2x_3d_p(k)%rAttr(p2x_3d_So_temp,n) &
                                                      +      TRACER(i,j,kuse,1,curtime,iblock)
                  p2x_3d_p(k)%rAttr(p2x_3d_So_salt,n) = p2x_3d_p(k)%rAttr(p2x_3d_So_salt,n) &
                                                      +      TRACER(i,j,kuse,2,curtime,iblock)
                  p2x_3d_p(k)%rAttr(p2x_3d_So_uvel,n) = p2x_3d_p(k)%rAttr(p2x_3d_So_uvel,n) &
                                                      +      WORK3_nrcm(i,j,kuse)
                  p2x_3d_p(k)%rAttr(p2x_3d_So_vvel,n) = p2x_3d_p(k)%rAttr(p2x_3d_So_vvel,n) &
                                                      +      WORK4_nrcm(i,j,kuse)
               endif
            enddo
         enddo
         enddo
      enddo

   !--- DEBUG ---
   if (dbug > 0) then
      k = p2x_2d_So_ssh
      write(stdout,*) subname,"p2x_2d_p    ssh  min,max: ",minval(p2x_2d_p   %rAttr(k,:)),maxval(p2x_2d_p   %rAttr(k,:))
      k = p2x_3d_So_temp
      write(stdout,*) subname,"p2x_3d_p(1) temp min,max= ",minval(p2x_3d_p(1)%rAttr(k,:)),maxval(p2x_3d_p(1)%rAttr(k,:))
      write(stdout,*) subname,"p2x_3d_p(2) temp min,max= ",minval(p2x_3d_p(2)%rAttr(k,:)),maxval(p2x_3d_p(1)%rAttr(k,:))
      call flushm (stdout)
   end if

   write(stdout,*) subname,"Exit" ; call flushm (stdout)

end subroutine ocpl_pop_export

!=========================================================================================
!=========================================================================================
end module ocpl_pop_mod
