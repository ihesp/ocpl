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
!  use POP_IOUnitsMod
   use POP_MCT_vars_mod

   use ocpl_fields_mod  
   use ocpl_data_mod

#ifdef OCPLRESTORING
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
   integer(int_kind) :: k,n                ! array index
#ifdef _OPENMP
   integer, external :: omp_get_max_threads  ! max number of threads that can execute
#endif

   character(*), parameter :: subName = "(ocpl_pop_init) "

!------------------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------------------

   write(o_logunit,'(2a)') subname,"Enter"

   nlev_p = km
   lsize = mct_aVect_lsize( p2x_p )

   !----- init 2d fields specifically for pop/roms coupling -----
   call mct_aVect_init(p2x_2d_p, rList=ocpl_fields_p2x_2d_fields,lsize=lsize)
   call mct_aVect_zero(p2x_2d_p)
   p2x_2d_p%rAttr(:,:) = 1.0e30

   !----- init 3d fields specifically for pop->roms coupling -----
   write(o_logunit,'(2a,i4)') subName,"     nlev_p    = ",nlev_p
   allocate(p2x_3d_p(nlev_p))

   do k = 1, nlev_p
      call mct_aVect_init(p2x_3d_p(k), rList=ocpl_fields_p2x_3d_fields,lsize=lsize)
      call mct_aVect_zero(p2x_3d_p(k))
      p2x_3d_p(k)%rAttr(:,:) = 1.0e30
   enddo

   !----- init 3d,2d fields specifically for roms->pop coupling -----
   write(o_logunit,'(2a,i4)') subName,"     nlev_rp    = ",nlev_rp
   call mct_aVect_init(r2x_2d_p, rList=ocpl_fields_r2x_2d_fields, lsize=lsize)
   call mct_aVect_zero(r2x_2d_p)
   allocate(r2x_3d_p(nlev_rp))
   do k = 1, nlev_rp
      call mct_aVect_init(r2x_3d_p(k), rList=ocpl_fields_r2x_3d_fields,lsize=lsize)
      call mct_aVect_zero(r2x_3d_p(k))
      p2x_3d_p(k)%rAttr(:,:) = 1.0e30
   enddo

   !--- set aVect field indicies ---
   k_p2x_2d_So_ssh  = mct_aVect_indexRA(p2x_2d_p   ,"So_ssh" )
   k_p2x_2d_So_ubar = mct_aVect_indexRA(p2x_2d_p   ,"So_ubar")
   k_p2x_2d_So_vbar = mct_aVect_indexRA(p2x_2d_p   ,"So_vbar")
   k_p2x_3d_So_temp = mct_aVect_indexRA(p2x_3d_p(1),"So_temp")
   k_p2x_3d_So_salt = mct_aVect_indexRA(p2x_3d_p(1),"So_salt")
   k_p2x_3d_So_uvel = mct_aVect_indexRA(p2x_3d_p(1),"So_uvel")
   k_p2x_3d_So_vvel = mct_aVect_indexRA(p2x_3d_p(1),"So_vvel")

   !--- DEBUG ---
   if (dbug > 0) then
      k = k_p2x_2d_So_ssh
      write(o_logunit,'(2a,2e11.3)') subname,"p2x_2d_p    ssh  min,max: ",minval(p2x_2d_p   %rAttr(k,:)),maxval(p2x_2d_p   %rAttr(k,:))
      k = k_p2x_3d_So_temp
      write(o_logunit,'(2a,2e11.3)') subname,"p2x_3d_p(1) temp min,max= ",minval(p2x_3d_p(1)%rAttr(k,:)),maxval(p2x_3d_p(1)%rAttr(k,:))
      write(o_logunit,'(2a,2e11.3)') subname,"p2x_3d_p(9) temp min,max= ",minval(p2x_3d_p(9)%rAttr(k,:)),maxval(p2x_3d_p(1)%rAttr(k,:))
      call flushm (o_logunit)
   end if

   write(o_logunit,'(2a)') subname,"Exit"

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

   type(mct_aVect) :: global_p  !  for debug global gather test
!EOP
!BOC

!-----------------------------------------------------------------------------------------
!  local variables

#ifdef _OPENMP
   integer, external :: omp_get_max_threads  ! max number of threads that can execute
#endif

   integer(IN) :: k,i,j,n,kuse
   integer(IN) :: lsize,ier
   integer(IN) :: master_task = 0     ! PID 0 is the master task
   real(R8)    :: x_min,x_max         ! debug test values

   integer(IN) :: iblock              ! block index
   type(block) :: this_block          ! block information for current block

   real(r8) :: WORK      (nx_block,ny_block,max_blocks_clinic) ! local work arrays
   real(R8) :: WORK1_ocpl(nx_block,ny_block)                   ! local 2d work space
   real(R8) :: WORK2_ocpl(nx_block,ny_block)                   ! local 2d work space
   real(R8) :: WORK3_ocpl(nx_block,ny_block,km)                ! local 3d work space
   real(R8) :: WORK4_ocpl(nx_block,ny_block,km)                ! local 3d work space

   character(*), parameter :: subName = "(ocpl_pop_export) "

!------------------------------------------------------------------------------
!  extract data from POP and put into corresponding attribute vectors
!  Q: do pop velocities need to be rotated?
!------------------------------------------------------------------------------

   write(o_logunit,'(2a)') subname,"Enter" ; call flushm (o_logunit)

   n = 0
   do iblock = 1,nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)

      ! move UPTROP and VBTROP to tgrid
      call ugrid_to_tgrid(WORK1_ocpl(:,:),UBTROP(:,:,curtime,iblock),iblock)
      call ugrid_to_tgrid(WORK2_ocpl(:,:),VBTROP(:,:,curtime,iblock),iblock)

      ! move UVEL and VVEL to tgrid
      do k = 1,km
         call ugrid_to_tgrid(WORK3_ocpl(:,:,k),UVEL(:,:,k,curtime,iblock),iblock)
         call ugrid_to_tgrid(WORK4_ocpl(:,:,k),VVEL(:,:,k,curtime,iblock),iblock)
      enddo

      do j = this_block%jb,this_block%je
      do i = this_block%ib,this_block%ie
         n = n + 1 ! local aVect index
         p2x_2d_p%rAttr(k_p2x_2d_So_ssh ,n) =      PSURF(i,j,curtime,iblock)/grav
         p2x_2d_p%rAttr(k_p2x_2d_So_ubar,n) = WORK1_ocpl(i,j)
         p2x_2d_p%rAttr(k_p2x_2d_So_vbar,n) = WORK2_ocpl(i,j)
         do k = 1,km
            if (KMT(i,j,iblock).gt.0) then
               kuse = min(k, KMT(i,j,iblock)) ! fill column from above
               p2x_3d_p(k)%rAttr(k_p2x_3d_So_temp,n) = TRACER(i,j,kuse,1,curtime,iblock)
               p2x_3d_p(k)%rAttr(k_p2x_3d_So_salt,n) = TRACER(i,j,kuse,2,curtime,iblock)
               p2x_3d_p(k)%rAttr(k_p2x_3d_So_uvel,n) = WORK3_ocpl(i,j,kuse)
               p2x_3d_p(k)%rAttr(k_p2x_3d_So_vvel,n) = WORK4_ocpl(i,j,kuse)
            endif
         enddo
      enddo
      enddo
   enddo

   !--- DEBUG ---
   if (dbug > 0) then
      !--- Note: tentitively gsMap_o is the pop gsMap and gsMap_p doesn't exist ---
      call mct_aVect_gather(p2x_2d_p, global_p, gsMap_o, master_task, mpicom_o, ier)

      if (seq_comm_iamroot(OCNID_o)) then
        lsize = mct_aVect_lsize(global_p)
        write(o_logunit,'(2a,i7)') subname,"global lsize = ",lsize
        write(o_logunit,'(2a,i7)') subname,"Note: vertical columns filled from above"
      end if

      if (seq_comm_iamroot(OCNID_o)) then
         k = k_p2x_2d_So_ssh
         x_min =  1.0e30
         x_max = -1.0e30
         do n=1,lsize
            if (global_p%rAttr(k,n) < 1.0e10) then
               x_min = min(x_min,global_p%rAttr(k,n))
               x_max = max(x_max,global_p%rAttr(k,n))
            end if
         end do
         write(o_logunit,'(2a,2es11.3)') subname,"global  ssh  min,max: ",x_min,x_max
      end if

      call mct_aVect_gather(p2x_3d_p(1), global_p, gsMap_o, master_task, mpicom_o, ier)
      if (seq_comm_iamroot(OCNID_o)) then
         k = k_p2x_3d_So_temp
         !---------------------
         x_min =  1.0e30
         x_max = -1.0e30
         do n=1,lsize
            if (global_p%rAttr(k,n) < 1.0e10) then
               x_min = min(x_min,global_p%rAttr(k,n))
               x_max = max(x_max,global_p%rAttr(k,n))
            end if
         end do
         write(o_logunit,'(2a,2es11.3)') subname,"global  sst  min,max: ",x_min,x_max
      end if
      call mct_aVect_gather(p2x_3d_p(9), global_p, gsMap_o, master_task, mpicom_o, ier)
      if (seq_comm_iamroot(OCNID_o)) then
         k = k_p2x_3d_So_temp
         !---------------------
         x_min =  1.0e30
         x_max = -1.0e30
         do n=1,lsize
            if (global_p%rAttr(k,n) < 1.0e10) then
               x_min = min(x_min,global_p%rAttr(k,n))
               x_max = max(x_max,global_p%rAttr(k,n))
            end if
         end do
         write(o_logunit,'(2a,2es11.3)') subname,"global T(9)  min,max: ",x_min,x_max
      end if
         
      call flushm (o_logunit)
   end if

   write(o_logunit,'(2a)') subname,"Exit" ; call flushm (o_logunit)

end subroutine ocpl_pop_export

!=========================================================================================
!=========================================================================================
end module ocpl_pop_mod
