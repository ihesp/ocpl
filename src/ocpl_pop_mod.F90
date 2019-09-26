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
!  public :: ocpl_pop_export
!  public :: ocpl_pop_import


! ! PUBLIC DATA:
!
! !REVISION HISTORY:
!    2019 Sep    - B. Kauffman, first version 
!
!EOP
! !PRIVATE MODULE VARIABLES


   real (r8),allocatable :: SBUFF_SUM(:,:,:,:) ! accum/sum/tavg export quantities
   real (r8) ::  tlast_coupled

   integer (int_kind)  ::   nsend, nrecv

   type(seq_infodata_type), pointer :: infodata   

   integer(IN) :: dbug = 0    ! debug level (higher is more)

!=========================================================================================
contains
!=========================================================================================

!BOP !====================================================================================

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

   nlev_p = km
   lsize = mct_aVect_lsize( p2x_p )

   !----- init 2d fields specifically for pop/roms coupling -----
   call mct_aVect_init(p2x_2d_p, rList=ocpl_fields_p2x_3d_fields,lsize=lsize)
   call mct_aVect_zero(p2x_2d_p)

   !----- init 3d fields specifically for pop/roms coupling -----
   write(stdout,*) subName,"     nlev_p    = ",nlev_p
   allocate(p2x_3d_p(nlev_p))

   do k = 1, nlev_p
      call mct_aVect_init(p2x_3d_p(k), rList=ocpl_fields_p2x_3d_fields,lsize=lsize)
      call mct_aVect_zero(p2x_3d_p(k))
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

end subroutine ocpl_pop_init
!=========================================================================================

 
!=========================================================================================
end module ocpl_pop_mod
