module ocpl_pop_mod

!=========================================================================================
!=========================================================================================

!BOP
! !MODULE: ocpl_pop_mod
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
   public :: ocpl_pop_import


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
   logical :: master
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
   master = (my_task == master_task)
   if(master) write(o_logunit,'(2a)') subname,"Enter"

   nlev_p = km
   lsize = mct_gsMap_lsize(gsMap_p, mpicom_p)
   if(master) write(o_logunit,'(2a,i4)') subName," gsMap_p lsize = ",lsize 
   lsize = mct_aVect_lsize( p2x_p )
   if(master) write(o_logunit,'(2a,i4)') subName,"     lsize     = ",lsize 

   !----- init 2d fields specifically for pop/roms coupling -----
   call mct_aVect_init(p2x_2d_p, rList=ocpl_fields_p2x_2d_fields,lsize=lsize)
   call mct_aVect_zero(p2x_2d_p)
   p2x_2d_p%rAttr(:,:) = 1.0e30

   !----- init 3d fields specifically for pop->roms coupling -----
   if(master) write(o_logunit,'(2a,i4)') subName,"     nlev_p    = ",nlev_p
   allocate(p2x_3d_p(nlev_p))

   do k = 1, nlev_p
      call mct_aVect_init(p2x_3d_p(k), rList=ocpl_fields_p2x_3d_fields,lsize=lsize)
      call mct_aVect_zero(p2x_3d_p(k))
      p2x_3d_p(k)%rAttr(:,:) = 1.0e30
   enddo

   !----- init 3d,2d fields specifically for roms->pop coupling -----
   call mct_aVect_init(r2x_2d_p, rList=ocpl_fields_r2x_2d_fields, lsize=lsize)
   call mct_aVect_zero(r2x_2d_p)

   if(master) write(o_logunit,'(2a,i4)') subName,"     nlev_rp   = ",nlev_rp
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

   if(master) write(o_logunit,'(2a)') subname,"Exit"

end subroutine ocpl_pop_init

!=========================================================================================
!BOP !====================================================================================
!
! !IROUTINE: ocpl_pop_import
!
! !DESCRIPTION:
!    put input data into pop, ocpl/pop 3d interface  
!
! !INTERFACE:
!
! !REVISION HISTORY:
!    2020 Jan -- Brian Kauffman (NCAR), initial version
!
! !INTERFACE: ----------------------------------------------------------------------------

subroutine ocpl_pop_import(p2x_p)   ! p2x_p used to put temporary debug fields onto cpl history file
!ubroutine ocpl_pop_import()

   use blocks
   use forcing_pt_interior, only: PT_INTERIOR_DATA     ! alter interal pop data for 2-way coupling
   use forcing_pt_interior, only: PT_RESTORE_RTAU      ! alter interal pop data for 2-way coupling
   use forcing_pt_interior, only: PT_RESTORE_MAX_LEVEL ! alter interal pop data for 2-way coupling
   use forcing_s_interior , only:  S_INTERIOR_DATA     ! alter interal pop data for 2-way coupling
   use forcing_s_interior , only:  S_RESTORE_RTAU      ! alter interal pop data for 2-way coupling
   use forcing_s_interior , only:  S_RESTORE_MAX_LEVEL ! alter interal pop data for 2-way coupling

! !INPUT/OUTPUT PARAMETERS:

   type(mct_aVect)             , intent(inout) :: p2x_p   

!EOP
!BOC

!  local variables
   type (block)            :: this_block     ! pop block information for current bloc
   integer(IN)             :: i,j,k,n,iblock,nCellsFromCoast
   real(r8)                :: rday ! restoring time-scale
   real(r8)                :: frac ! roms cell fraction mapped to pop grid
   logical                 :: first_call = .true.
   logical                 :: pop_restoring  = .true. ! normally true, false only for DEBUG
   integer(IN)             :: dbug_save
   character(*), parameter :: subName = "(ocpl_pop_import) "
   integer                 :: n_save   ! Jaison, Sep/20/2021, with suggestions from Brian & Abishek.
 
!-----------------------------------------------------------------------------------------
!
!-----------------------------------------------------------------------------------------

   dbug_save = dbug ! debug this routine (only?)
   if (first_call) dbug = 2

   if (master .or. dbug > 0) write(o_logunit,'(2a)') subName,"Enter" 

   if (master .or. dbug > 1) then !----- optional debug info -----
      write(o_logunit,'(2a,2e12.3)') subName,"min/max r2x_2d_p%rAttr(k_r2x_2d_wgts,      :)) = " &
                                    ,minval(r2x_2d_p   %rAttr(k_r2x_2d_wgts   ,:)) &
                                    ,maxval(r2x_2d_p   %rAttr(k_r2x_2d_wgts   ,:))
      write(o_logunit,'(2a,2e12.3)') subName,"min/max r2x_2d_p%rAttr(k_r2x_2d_frac,      :)) = " &
                                    ,minval(r2x_2d_p   %rAttr(k_r2x_2d_frac   ,:)) &
                                    ,maxval(r2x_2d_p   %rAttr(k_r2x_2d_frac   ,:))
      write(o_logunit,'(2a,2e12.3)') subName,"min/max r2x_3d_p(1)%rAttr(k_r2x_3d_So_temp,:)) = " &
                                    ,minval(r2x_3d_p(1)%rAttr(k_r2x_3d_So_temp,:)) &
                                    ,maxval(r2x_3d_p(1)%rAttr(k_r2x_3d_So_temp,:))
      write(o_logunit,'(2a,2e12.3)') subName,"min/max r2x_3d_p(1)%rAttr(k_r2x_3d_So_salt,:)) = " &
                                    ,minval(r2x_3d_p(1)%rAttr(k_r2x_3d_So_salt,:)) &
                                    ,maxval(r2x_3d_p(1)%rAttr(k_r2x_3d_So_salt,:))
      write(o_logunit,'(2a,i3)'    ) subName,"nblocks_clinic = ",nblocks_clinic
   end if

   if (master .and. first_call) write(o_logunit,'(2a,L2)') subName,"pop 3d restoring = ",pop_restoring

   if (pop_restoring) then 
   !! Sep/17/2021: In the previous version of this program, n=0 step was done within the do loop on "iblock".
   !!                Which was okay for a case with nblocks_clinic=1 (total blocks on a processor) but not for 
   !!                cases with nblocks_clinic > 1. The new variable n_save fixes this issue by using correct
   !!                initial value for "n" for each block within the "iblock" do loop.
   !!                Jaison, Abishek and Brian, Sep/20/2021.
   n_save = 0  
   do iblock = 1, nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)

      if (dbug > 1) then !----- optional debug info -----
         write(o_logunit,'(2a,i3)'    ) subName,"iblock = ",iblock
         write(o_logunit,'(2a,2e12.3)') subName,"orig  min/max PT_RESTORE_RTAU = ",minval(PT_RESTORE_RTAU     (:,:  ,iblock)  ) &
                                                                                  ,maxval(PT_RESTORE_RTAU     (:,:  ,iblock)  )
         write(o_logunit,'(2a,2e12.3)') subName,"orig  min/max  S_RESTORE_RTAU = ",minval( S_RESTORE_RTAU     (:,:  ,iblock)  ) &
                                                                                  ,maxval( S_RESTORE_RTAU     (:,:  ,iblock)  )
         write(o_logunit,'(2a,2i12  )') subName,"orig  min/max PT_RESTORE_MAX_L= ",minval(PT_RESTORE_MAX_LEVEL(:,:  ,iblock  )) &
                                                                                  ,maxval(PT_RESTORE_MAX_LEVEL(:,:  ,iblock  ))
         write(o_logunit,'(2a,2i12  )') subName,"orig  min/max  S_RESTORE_MAX_L= ",minval( S_RESTORE_MAX_LEVEL(:,:  ,iblock  )) &
                                                                                  ,maxval( S_RESTORE_MAX_LEVEL(:,:  ,iblock  )) 
         write(o_logunit,'(2a,2e12.3)') subName,"orig  min/max  S_INTERIOR_DATA= ",minval( S_INTERIOR_DATA    (:,:,1,iblock,1)) &
                                                                                  ,maxval( S_INTERIOR_DATA    (:,:,1,iblock,1))
         write(o_logunit,'(2a,2e12.3)') subName,"orig  min/max PT_INTERIOR_DATA= ",minval(PT_INTERIOR_DATA    (:,:,1,iblock,1)) &
                                                                                  ,maxval(PT_INTERIOR_DATA    (:,:,1,iblock,1))
         call shr_sys_flush(o_logunit)
      end if

      !-----------------------------------------------------------------------------------
      ! one-time setup of restoring depth and timescale
      !-----------------------------------------------------------------------------------
      if (first_call) then 
         if(master) write(o_logunit,'(2a)') subName,"first call: set restoring max levels and rtau"

         n = n_save ! Sep/20/2021: Use the maximum value from last iteration
         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1

            if (r2x_2d_p%rAttr(k_r2x_2d_frac,n) > 0.01) then ! roms data exists for this pop cell

               !----- taper-off restoring near edge of restoring region -----
               !----- there's probably a better way of doing this -----
               rday =   1000       ! extremely weak
               if (r2x_2d_p%rAttr(k_r2x_2d_wgts,n) > 0.0 ) rday = 200  ! very weak
               if (r2x_2d_p%rAttr(k_r2x_2d_wgts,n) > 0.3 ) rday = 100
               if (r2x_2d_p%rAttr(k_r2x_2d_wgts,n) > 0.6 ) rday =  50
               if (r2x_2d_p%rAttr(k_r2x_2d_wgts,n) > 0.9 ) rday =  10  ! very strong
               PT_RESTORE_RTAU(i,j,iblock) = 1.0_r8 / (  rday * 86400.0_r8)

               !----- restoring depth limited to surface layer near coast -----
               nCellsFromCoast = PT_RESTORE_MAX_LEVEL(i,j,iblock) ! on input max level identifies coastal cells
               if (nCellsFromCoast < 1) then
                  if(master) write(o_logunit,'(2a,2i6,a   )') subName,"i,j =",i,j,", land cell"
                  PT_RESTORE_MAX_LEVEL(i,j,iblock) = 0
               else if (nCellsFromCoast < 4) then
                  if(master) write(o_logunit,'(2a,2i6,a,i6)') subName,"i,j =",i,j,", coastal cell, distance from coast =",nCellsFromCoast
                  PT_RESTORE_MAX_LEVEL(i,j,iblock) = min(nCellsFromCoast,int(r2x_2d_p%rAttr(k_r2x_2d_reslev,n) ))
                ! PT_RESTORE_MAX_LEVEL(i,j,iblock) = 1 ! limit PT restoring to surface layer only (?)
               else
                  if(master) write(o_logunit,'(2a,2i6,a   )') subName,"i,j =",i,j,", cells from coast > 3"
                  PT_RESTORE_MAX_LEVEL(i,j,iblock) = max(1,int(r2x_2d_p%rAttr(k_r2x_2d_reslev,n) ))
               end if

            else ! no roms data available to restore to

               rday = -1
               PT_RESTORE_MAX_LEVEL(i,j,iblock) = 0
               PT_RESTORE_RTAU     (i,j,iblock) = 0.0_r8  ! indicates no restoring

            end if ! roms data exist for the cell

            !----- temporary debug/valication fields (don't commit, requires non-standard coupler) -----
            p2x_p%rAttr(mct_aVect_indexRA(p2x_p,"So_rday"  ),n) = rday
            p2x_p%rAttr(mct_aVect_indexRA(p2x_p,"So_rtau"  ),n) = PT_RESTORE_RTAU     (i,j,iblock)
            p2x_p%rAttr(mct_aVect_indexRA(p2x_p,"So_maxlev"),n) = PT_RESTORE_MAX_LEVEL(i,j,iblock)
            p2x_p%rAttr(mct_aVect_indexRA(p2x_p,"So_frac"  ),n) = r2x_2d_p%rAttr(k_r2x_2d_frac,n)
            p2x_p%rAttr(mct_aVect_indexRA(p2x_p,"So_wgts"  ),n) = r2x_2d_p%rAttr(k_r2x_2d_wgts,n)

            !----- salinity gets same treatment as temperature ----- 
            S_RESTORE_RTAU     (i,j,iblock) = PT_RESTORE_RTAU     (i,j,iblock) ! use same restoring as PT
            S_RESTORE_MAX_LEVEL(i,j,iblock) = PT_RESTORE_MAX_LEVEL(i,j,iblock)

         enddo ! i
         enddo ! j

         !----- for testing: turns off pop restoring -----
         if (.false.) then 
            write(o_logunit,'(2a)') subName,"DEBUG: pop restoring DEACTIVATED"
            PT_RESTORE_MAX_LEVEL(:,:,iblock) = 0
            PT_RESTORE_RTAU     (:,:,iblock) = 0.0_r8
         end if 

      end if   ! first_call

      !-----------------------------------------------------------------------------------
      ! update PT & S values that pop will restore to (RTAU and MAX_LEVEL are constant)
      !-----------------------------------------------------------------------------------
      n = n_save       ! Sep/20/2021: Use the maximum value from last iteration (do not consider
                       !              the iterations done under "first_call" abvoe, that is just to 
                       !              initialize few arrays).
      do j=this_block%jb,this_block%je
      do i=this_block%ib,this_block%ie
         n = n + 1

         if (r2x_2d_p%rAttr(k_r2x_2d_frac,n) > 0.01) then ! roms data exists for this pop cell
            do k=1,PT_RESTORE_MAX_LEVEL(i,j,iblock)
         !  do k=1,nLev_rp
                PT_INTERIOR_DATA(i,j,k,iblock,1) = r2x_3d_p(k)%rAttr(k_r2x_3d_So_temp,n) - T0_kelvin
                 S_INTERIOR_DATA(i,j,k,iblock,1) = r2x_3d_p(k)%rAttr(k_r2x_3d_So_salt,n) ! * 0.001_r8
            end do
         end if

      enddo ! i
      enddo ! j

      !-----------------------------------------------------------------------------------
      ! optional debug info
      !-----------------------------------------------------------------------------------
      if (master .or. dbug > 1) then
         write(o_logunit,'(2a,2e12.3)') subName,"after min/max PT_RESTORE_RTAU = ",minval(PT_RESTORE_RTAU     (:,:  ,iblock)  ) &
                                                                                  ,maxval(PT_RESTORE_RTAU     (:,:  ,iblock)  )
         write(o_logunit,'(2a,2e12.3)') subName,"after min/max  S_RESTORE_RTAU = ",minval( S_RESTORE_RTAU     (:,:  ,iblock)  ) &
                                                                                  ,maxval( S_RESTORE_RTAU     (:,:  ,iblock)  )
         write(o_logunit,'(2a,2i12  )') subName,"after min/max PT_RESTORE_MAX_L= ",minval(PT_RESTORE_MAX_LEVEL(:,:  ,iblock  )) &
                                                                                  ,maxval(PT_RESTORE_MAX_LEVEL(:,:  ,iblock  ))
         write(o_logunit,'(2a,2i12  )') subName,"after min/max  S_RESTORE_MAX_L= ",minval( S_RESTORE_MAX_LEVEL(:,:  ,iblock  )) &
                                                                                  ,maxval( S_RESTORE_MAX_LEVEL(:,:  ,iblock  ))
         write(o_logunit,'(2a,2e12.3)') subName,"after min/max PT_INTERIOR_DATA= ",minval(PT_INTERIOR_DATA    (:,:,1,iblock,1)) &
                                                                                  ,maxval(PT_INTERIOR_DATA    (:,:,1,iblock,1))
         write(o_logunit,'(2a,2e12.3)') subName,"after min/max  S_INTERIOR_DATA= ",minval( S_INTERIOR_DATA    (:,:,1,iblock,1)) &
                                                                                  ,maxval( S_INTERIOR_DATA    (:,:,1,iblock,1))
      end if

      n_save = n       ! Sep/20/2021: Update n_save for next iteration
   enddo ! iblock
   end if ! pop_restoring

   first_call = .false.

   dbug = dbug_save

end subroutine ocpl_pop_import

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

   !----- local variables -----

   type(mct_aVect) :: global_p  !  for debug global gather test

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

   if(master) write(o_logunit,'(2a)') subname,"Enter" ; call flushm (o_logunit)

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
      !--- Note: tentitively gsMap_p is the pop gsMap and gsMap_p doesn't exist ---
      call mct_aVect_gather(p2x_2d_p, global_p, gsMap_p, master_task, mpicom_p, ier)

      if (seq_comm_iamroot(OCNID_p)) then
        lsize = mct_aVect_lsize(global_p)
        write(o_logunit,'(2a,i7)') subname,"global lsize = ",lsize
        write(o_logunit,'(2a,i7)') subname,"Note: vertical columns filled from above"
      end if

      if (seq_comm_iamroot(OCNID_p)) then
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

      call mct_aVect_gather(p2x_3d_p(1), global_p, gsMap_p, master_task, mpicom_p, ier)
      if (seq_comm_iamroot(OCNID_p)) then
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
      call mct_aVect_gather(p2x_3d_p(9), global_p, gsMap_p, master_task, mpicom_p, ier)
      if (seq_comm_iamroot(OCNID_p)) then
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

   if(master) write(o_logunit,'(2a)') subname,"Exit" 

end subroutine ocpl_pop_export

!=========================================================================================
!=========================================================================================
end module ocpl_pop_mod
