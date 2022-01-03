module ocn_comp_mct

!================================================================================
!================================================================================
!BOP
! !MODULE: ocn_comp_mct
! !INTERFACE:

! !DESCRIPTION:
!  This is the ocean driver that interfaces with the CESM Flux Coupler
!
!  This ocean driver also coordinates the interaction of the underlying
!  prognotic ocean models pop and roms, with pop running globally and
!  roms running regionaly, roms embedded into pop, and the two-way
!  coupling between roms and pop.  This ocean driver interacts with the 
!  Flux Coupler in the same way as any other plug-and-play CESM ocean --
!  from the Flux Coupler's point of view, it cannot detect this ocean
!  component is in fact two ocean models, one embedded into the other.
!
! !REVISION HISTORY:
!    2018-Sep: first version by Brian Kauffman 
!
! !USES:
   use pop_comp_mct , only :  pop_init_mct, pop_run_mct, pop_final_mct
   use roms_comp_mct, only : roms_init_mct,roms_run_mct,roms_final_mct
   use io_types     , only :  pop_logUnit => stdout  ! ocpl redirects logUnit for pop
   use mod_iounits  , only : roms_logUnit => stdout  ! ocpl redirects logUnit for roms

   use mct_mod
   use shr_mct_mod
   use esmf

   use seq_flds_mod
   use seq_cdata_mod
   use seq_infodata_mod
   use seq_timemgr_mod
!  use seq_comm_mct,      only : seq_comm_suffix, seq_comm_inst, seq_comm_name
   use seq_comm_mct

   use shr_file_mod 
   use shr_cal_mod,       only : shr_cal_date2ymd
   use shr_mpi_mod,       only : shr_mpi_bcast
   use shr_kind_mod, only: IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL

   use ocpl_data_mod
   use ocpl_domain_mod
   use ocpl_map_mod
   use ocpl_pop_mod
   use ocpl_roms_mod
   use ocpl_fields_mod ! for temporary debug/validation info

! !PUBLIC MEMBER FUNCTIONS:

  implicit none
  private              ! default: all data is private
  SAVE                 ! save everything

  public :: ocn_init_mct
  public :: ocn_run_mct
  public :: ocn_final_mct

! ! PUBLIC DATA:

! ! PRIVATE DATA:

   integer(IN) :: shrlogunit,shrloglev  ! to save & restore shared log units & levels
   integer(IN) :: debug = 0            ! debug level

!=========================================================================================
contains
!=========================================================================================
!BOP =====================================================================================
!
! !IROUTINE: ocn_init_mct
!
! !DESCRIPTION:
!     initialize ocean component: composite of pop & ROMS 
!     In this version, the ocean grid the coupler knows about is the pop grid
!
! !REVISION HISTORY:
!    2018 Oct -- Brian Kauffman, original code
!
! !INTERFACE:
!-----------------------------------------------------------------------------------------

  subroutine ocn_init_mct( EClock, cdata_o, x2o_o, o2x_o, NLFilename )


   implicit none

! INPUT/OUTPUT PARAMETERS:

   type(ESMF_Clock)            , intent(inout) :: EClock
   type(seq_cdata)             , intent(inout) :: cdata_o
   type(mct_aVect)             , intent(inout) :: x2o_o, o2x_o
   character(len=*), optional  , intent(in)    :: NLFilename ! Namelist filename

 
   !--- local variables ---
   integer(IN) :: ioffset, joffset
   integer(IN) :: ni_o, nj_o
   integer(IN) :: ni_p, nj_p
   integer(IN) :: ni_r, nj_r
   integer(IN) :: rc
   integer(IN) :: k       ! aVect field index 
   integer(IN) :: n       ! aVect cell index 
   integer(IN) :: lSize_o ! aVect local size of *_o grid
   integer(IN) :: lSize_p ! aVect local size of *_p grid
   integer(IN) :: lSize_r ! aVect local size of *_r grid

   logical :: restart

   character(CL) :: start_type
   character(CS) :: case_name

   integer(IN) :: ymd,tod
   integer(IN) :: ncomp

   character(*), parameter :: subName = "(ocn_init_mct) "
   character(*), parameter :: F00 =  "( '(ocn_init_mct) ===== ',a,' ',70('=') )"
   character(*), parameter :: F01 =  "( '(ocn_init_mct) ----- ',a,' ',70('-') )"

!-----------------------------------------------------------------------------------------
! Note: the seq_cdata data type contains 
!    type seq_cdata
!       character(len=16)                :: name               ! user defined name
!       integer                          :: ID                 ! component id
!       integer                          :: mpicom             ! mpi communicator
!       type(mct_gGrid)         ,pointer :: dom => null()      ! domain info
!       type(mct_gsMap)         ,pointer :: gsMap => null()    ! decomp info
!       type(seq_infodata_type) ,pointer :: infodata => null() ! control flags to/from cpl
!    end type seq_cdata
!
! The ocpl, pop, & roms components have own cdata_o, cdata_p, & cdata_r
! which contain their own separate/different name, gGrid, and gsMap (comps alloc/init these)
! but their ID & mpicom must be the same (the are all part of the same ocean component)
! and their infodata is not only the same, but literally pointing to the same infodata data 
! in memory (cpl allocs/inits this data type)
!-----------------------------------------------------------------------------------------

   !--- set stdout/log units ---
   call shr_file_getLogUnit (shrlogunit) ! save log unit
   call shr_file_getLogLevel(shrloglev)  ! save log level
   call seq_cdata_setptrs(cdata_o, ID=OCNID_o)
   if (seq_comm_iamroot(OCNID_o)) then       ! set logUnit to ocn.log.*
      o_logUnit = shr_file_getUnit()        ! get/reserve unused unit number for ocpl's logUnit
      call shr_file_setIO("ocn_modelio.nml",o_logUnit)  ! assign ocpl logUnit to a named file
      call shr_file_setLogUnit (o_logunit)  ! set shared logUnit = ocpl's logUnit
      call shr_file_setLogLevel(1)
      d_logunit = o_logUnit                 ! where debug info goes
   else
      o_logUnit = shrlogunit ! wherever cpl was using for non-root logunit, presumably cesm.log.*
      d_logunit = o_logUnit  ! where debug info goes
      if (debug==0) then     ! send o_logunit to /dev/null
         d_logunit = shr_file_getUnit() ! get/reserve unused unit number
         open(d_logunit,file="/dev/null",status='replace',action='write')
      endif
   endif
   pop_logUnit  = o_logUnit              ! set pop    logUnit = ocpl's logUnit
   roms_logUnit = o_logUnit              ! set roms   logUnit = ocpl's logUnit

   ncomp = OCNID(1)
   if (seq_comm_iamroot(OCNID_o)) then       ! set logUnit to ocn.log.*
      write(o_logunit,F00) " ENTER"
      write(o_logunit,*) subName,"comp  ID = ",                ncomp
      write(o_logunit,*) subName,"comp PID = ",seq_comm_iam   (ncomp)
      write(o_logunit,*) subName,"glob PID = ",seq_comm_gloiam(ncomp)
      write(o_logunit,*) subName,"logUnit  = ",o_logUnit
   endif
   !--- Get data/pointers out of cdata_o ---
   call seq_cdata_setptrs(cdata_o, ID=OCNID_o, dom=gGrid_o, gsMap=gsMap_o, infodata=infodata_o, mpicom=mpicom_o,name=name_o)

   if (seq_comm_iamroot(OCNID_o)) then       ! set logUnit to ocn.log.*
      write(o_logunit,*) subName, 'cdata_o ID     : ' ,  OCNID_o
      write(o_logunit,*) subName, 'cdata_o mpicom : ' ,  mpicom_o
      write(o_logunit,*) subName, 'cdata_o name   : ' // trim(name_o )
   endif
   !--- Get data out of infodata ---
   call seq_infodata_GetData(infodata_o, case_name=case_name, start_type=start_type)
   if (seq_comm_iamroot(OCNID_o)) then       ! set logUnit to ocn.log.*
      write(o_logunit,*) subName, 'case_name  : ' // trim(case_name )
      write(o_logunit,*) subName, 'start_type : ' // trim(start_type)

   !--------------------------------------------------------------------------------------
   ! init ocpl/ocn gsMap and domain
   !--------------------------------------------------------------------------------------
      write(o_logunit,F01) "ocpl_domain_init call"
   endif
   call ocpl_domain_init() 
   if (seq_comm_iamroot(OCNID_o)) then       ! set logUnit to ocn.log.*
      write(o_logunit,F01) "ocpl_domain_init return"
   endif
   call seq_infodata_GetData(infodata_o, ocn_nx = ni_o, ocn_ny = nj_o)

   !--- init x2o_o & o2x_o (now that domain is initialized so we know lsize) ---
   lsize_o = mct_gsMap_lsize(gsMap_o, mpicom_o)
   call mct_aVect_init(x2o_o, rList=seq_flds_x2o_fields, lsize=lsize_o)
   call mct_aVect_zero(x2o_o)
   call mct_aVect_init(o2x_o, rList=seq_flds_o2x_fields, lsize=lsize_o)
   call mct_aVect_zero(o2x_o)

   !--- sanity check on some data ---
   call seq_cdata_setptrs(cdata_o, ID=OCNID_o, mpicom=mpicom_o,name=name_o)
   if (seq_comm_iamroot(OCNID_o)) then       ! set logUnit to ocn.log.*
      write(o_logunit,'(3a    )') subName,'cdata_o name       = ' ,trim(name_o )
      write(o_logunit,'(2a,2i5)') subName,'cdata_o ID         = ' ,OCNID_o
      write(o_logunit,'(2a,2i5)') subName,'cdata_o mpicom     = ' ,mpicom_o
      write(o_logunit,'(2a,2i7)') subName,'ni_o,nj_o          = ' ,ni_o,nj_o
      write(o_logunit,'(2a,2i7)') subname,'x2o_o lsize, nflds = ' ,mct_aVect_lsize(x2o_o),mct_aVect_nRAttr(x2o_o)
      write(o_logunit,'(2a,2i7)') subname,'o2x_o lsize, nflds = ' ,mct_aVect_lsize(o2x_o),mct_aVect_nRAttr(o2x_o)
      write(o_logunit,'(2a,2i7)') subname,'gsMap_o lsize      = ' , mct_gsMap_lsize(gsMap_o, mpicom_o)
   endif
   !--------------------------------------------------------------------------------------
   ! call pop  initialize phase -- note: ocpl, roms & pop share ID, mpicom & infodata
   !--------------------------------------------------------------------------------------

   !--- init cdata_p --- contains ID,mpicom,infdodata,gGrid,gsMap 
   OCNID_p = OCNID_o  
   call seq_cdata_init(cdata_p, ID=OCNID_p, dom=gGrid_p, gsMap=gsMap_p, infodata=infodata_o, name='POP')
   call seq_cdata_setptrs(cdata_p, mpicom=mpicom_p) ! seq_cdata_init sets mpicom based on ID
   !  --- alternative to calling seq_cdata_setptrs ---
   !  mpicom_p = mpicom_o
   !  cdata_p%name      =  "POP"
   !  cdata_p%ID        =  OCNID_p
   !  cdata_p%mpicom    =  mpicom_p
   !  cdata_p%dom       => gGrid_p
   !  cdata_p%gsmap     => gsMap_p
   !  cdata_p%ginfodata => infodata_o

   !--- initialize pop --
   if (seq_comm_iamroot(OCNID_o)) then       ! set logUnit to ocn.log.*
      write(o_logunit,F01) "pop_init_mct call" 
   endif
   call pop_init_mct( EClock, cdata_p, x2o_p, p2x_p, NLFilename ) 
   if (seq_comm_iamroot(OCNID_o)) write(o_logunit,F01) "pop_init_mct return" 
   call seq_timemgr_EClockGetData( EClock, curr_ymd=ymd, curr_tod=tod )
   if (seq_comm_iamroot(OCNID_o)) write(o_logunit,*) subname,'pop  ymd, tod =',ymd,tod

   !--- restore ocpl/ocn dims in infodata, were changed to pop dims by pop_init_mct() ---
   call seq_infodata_GetData(infodata_o, ocn_nx = ni_p, ocn_ny = nj_p)
   call seq_infodata_PutData(infodata_o, ocn_nx = ni_o, ocn_ny = nj_o) ! restore ocpl dims

   !--- init aVect for pop output on ocpl/ocn grid ----
   lsize_o = mct_aVect_lsize(o2x_o) 
   call mct_aVect_init(p2x_o, p2x_p, lsize=lsize_o)
   call mct_aVect_zero(p2x_o)

   !--- sanity check on some data ---
   call seq_cdata_setptrs(cdata_p, ID=OCNID_p, mpicom=mpicom_p,name=name_p)
   if (seq_comm_iamroot(OCNID_o)) then       ! set logUnit to ocn.log.*
      write(o_logunit,'(3a    )') subName,'cdata_p name       = ' ,trim(name_p )
      write(o_logunit,'(2a,2i5)') subName,'cdata_p ID         = ' ,OCNID_p
      write(o_logunit,'(2a,2i5)') subName,'cdata_p mpicom     = ' ,mpicom_p
      write(o_logunit,'(2a,2i5)') subName,'ni_p,nj_p          = ' ,ni_p,nj_p
      write(o_logunit,'(2a,2i7)') subname,'x2o_p lsize, nflds = ' ,mct_aVect_lsize(x2o_p),mct_aVect_nRAttr(x2o_p)
      write(o_logunit,'(2a,2i7)') subname,'p2x_p lsize, nflds = ' ,mct_aVect_lsize(p2x_p),mct_aVect_nRAttr(p2x_p)
      write(o_logunit,'(2a,2i7)') subname,'p2x_o lsize, nflds = ' ,mct_aVect_lsize(p2x_o),mct_aVect_nRAttr(p2x_o)
      write(o_logunit,'(2a,2i7)') subname,'gsMap_p lsize      = ' , mct_gsMap_lsize(gsMap_p, mpicom_p)
   endif
   !--------------------------------------------------------------------------------------
   ! call roms initialize phase -- note: ocpl, roms & pop share ID, mpicom & infodata
   !--------------------------------------------------------------------------------------

   !--- init cdata_p --- contains ID,mpicom,infdodata,gGrid,gsMap 
   OCNID_r = OCNID_o  
   call seq_cdata_init(cdata_r, ID=OCNID_r, dom=gGrid_r, gsMap=gsMap_r, infodata=infodata_o, name='ROMS')
   call seq_cdata_setptrs(cdata_r, mpicom=mpicom_r) ! seq_cdata_init sets mpicom based on ID

   !--- initialize roms ---
   if (seq_comm_iamroot(OCNID_o)) write(o_logunit,F01) "roms_init_mct call" 
   call roms_init_mct( EClock, cdata_r, x2o_r, r2x_r)  ! , NLFilename )
   if (seq_comm_iamroot(OCNID_o)) write(o_logunit,F01) "roms_init_mct return" 
   call seq_timemgr_EClockGetData( EClock, curr_ymd=ymd, curr_tod=tod )
   if (seq_comm_iamroot(OCNID_o)) write(o_logunit,*) subname,'roms ymd, tod =',ymd,tod

   !--- restore ocpl dims in infodata, were changed to regional dims by roms_init_mct() ---
   call seq_infodata_GetData(infodata_o, ocn_nx = ni_r, ocn_ny = nj_r) ! get roms dims
   call seq_infodata_PutData(infodata_o, ocn_nx = ni_o, ocn_ny = nj_o) ! restore ocpl dims

   !--- init aVect for roms output on ocpl grid ----
   lsize_o = mct_aVect_lsize(o2x_o)
   call mct_aVect_init(r2x_o, r2x_r, lsize=lsize_o)
   call mct_aVect_zero(r2x_o)

   !--- sanity check on some data ---
   call seq_cdata_setptrs(cdata_r, ID=OCNID_r, mpicom=mpicom_r,name=name_r)
   if (seq_comm_iamroot(OCNID_o)) then       ! set logUnit to ocn.log.*
      write(o_logunit,'(3a    )') subName,'cdata_r name       = ' ,trim(name_r )
      write(o_logunit,'(2a,2i5)') subName,'cdata_r ID         = ' ,OCNID_r
      write(o_logunit,'(2a,2i5)') subName,'cdata_r mpicom     = ' ,mpicom_r
      write(o_logunit,'(2a,2i5)') subName,'ni_r,nj_r          = ' ,ni_r,nj_r
      write(o_logunit,'(2a,2i7)') subname,'x2o_r lsize, nflds = ' ,mct_aVect_lsize(x2o_r),mct_aVect_nRAttr(x2o_r)
      write(o_logunit,'(2a,2i7)') subname,'r2x_r lsize, nflds = ' ,mct_aVect_lsize(r2x_r),mct_aVect_nRAttr(r2x_r)
      write(o_logunit,'(2a,2i7)') subname,'r2x_o lsize, nflds = ' ,mct_aVect_lsize(r2x_o),mct_aVect_nRAttr(r2x_o)
      write(o_logunit,'(2a,2i7)') subname,'gsMap_r lsize      = ' , mct_gsMap_lsize(gsMap_r, mpicom_r)
   endif
   !--------------------------------------------------------------------------------------
   ! init additional data-types needed for ocpl's 3d global/regional ocean coupling 
   !--------------------------------------------------------------------------------------
   
   if (seq_comm_iamroot(OCNID_o)) write(o_logunit,F01) "call ocpl_pop_init" 
   call ocpl_pop_init( p2x_p, p2x_2d_p, p2x_3d_p)

   if (seq_comm_iamroot(OCNID_o)) write(o_logunit,F01) "call ocpl_roms_init" 
   call ocpl_roms_init()

   !--------------------------------------------------------------------------------------
   ! init all maps (surface, 3d, curtain)
   !--------------------------------------------------------------------------------------
   if (seq_comm_iamroot(OCNID_o)) write(o_logunit,*) subName, "initialize all maps..."
   call ocpl_map_init()

   !--------------------------------------------------------------------------------------
   ! map & merge pop & roms output to create ocpl output
   !--------------------------------------------------------------------------------------
   if (seq_comm_iamroot(OCNID_o)) write(o_logunit,'(2a)') subname,"map: r2x_r -> r2x_o"
   call mct_sMat_avMult(r2x_r, sMatp_r2o, r2x_o,vector=usevector)

   if (seq_comm_iamroot(OCNID_o)) write(o_logunit,'(2a)') subname,"map: p2x_p -> p2x_o"
   call mct_sMat_avMult(p2x_p, sMatp_p2o, p2x_o,vector=usevector)

   if (debug > 0) then
      k = mct_avect_indexra(r2x_o,'So_t')
      write(o_logunit,'(2a,2e12.4)') subname,'<DEBUG> min/max r2x_r SST = ',minval(r2x_r%rAttr(k,:)),maxval(r2x_r%rAttr(k,:))
      write(o_logunit,'(2a,2e12.4)') subname,'<DEBUG> min/max r2x_o SST = ',minval(r2x_o%rAttr(k,:)),maxval(r2x_o%rAttr(k,:))
      write(o_logunit,'(2a,2e12.4)') subname,'<DEBUG> min/max p2x_p SST = ',minval(p2x_p%rAttr(k,:)),maxval(p2x_p%rAttr(k,:))
      write(o_logunit,'(2a,2e12.4)') subname,'<DEBUG> min/max o2x_o SST = ',minval(p2x_o%rAttr(k,:)),maxval(p2x_o%rAttr(k,:))
   end if

   if (seq_comm_iamroot(OCNID_o)) write(o_logunit,'(a)') subname,"merge: r2x_o -> o2x_o"
   k = mct_avect_indexra(o2x_o,'So_t')
   lsize_o = mct_aVect_lsize( o2x_o )
   do n=1,lsize_o
      o2x_o%rAttr(:,n) = p2x_o%rAttr(:,n)  ! start with all pop data, replace some with roms data
      if (r2x_o%rAttr(k,n) > 1.0) then  ! has data mapped from roms (unmapped cells have sst = 0)
          o2x_o%rAttr(:,n) = r2x_o%rAttr(:,n)  ! merge all fields
      !   o2x_o%rAttr(k,n) = r2x_o%rAttr(k,n)  ! merge SST only
      end if
   end do

   if (debug > 0 ) then
      k = mct_avect_indexra(o2x_o,'So_t')
      write(o_logunit,'(2a,2e12.4)') subname,'<DEBUG> min/max o2x_o SST+ = ',minval(o2x_o%rAttr(k,:)),maxval(o2x_o%rAttr(k,:))
   end if

   !--------------------------------------------------------------------------------------
   ! done
   !--------------------------------------------------------------------------------------
   if (seq_comm_iamroot(OCNID_o)) write(o_logunit,F00) "EXIT"

   !--- Reset shr logging to original values ---
   call shr_file_setLogUnit (shrlogunit) ! restore log unit
   call shr_file_setLogLevel(shrloglev)  ! restore log level

 end subroutine ocn_init_mct

!=========================================================================================
!=========================================================================================
!BOP
!
! !IROUTINE: ocn_run_mct
!
! !INTERFACE:
  subroutine ocn_run_mct( EClock, cdata_o, x2o_o, o2x_o)
!
! !DESCRIPTION:
! Run POP for a coupling interval
!
! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)            , intent(inout) :: EClock
    type(seq_cdata)             , intent(inout) :: cdata_o
    type(mct_aVect)             , intent(inout) :: x2o_o
    type(mct_aVect)             , intent(inout) :: o2x_o

!
! !REVISION HISTORY:
! Author: Brian Kauffman
!EOP


   !----- local variables -----
   integer :: lbnum
   integer(IN) :: ymd,tod
   integer(IN) :: k       ! aVect field index 
   integer(IN) :: n       ! aVect cell index 
   integer(IN) :: lSize_o ! aVect local size of *_o grid

   character(*), parameter :: SubName = "(ocn_run_mct) "
   character(*), parameter :: F00 =  "( '(ocn_run_mct) ===== ',a,' ',70('=') )"
   character(*), parameter :: F01 =  "( '(ocn_run_mct) ----- ',a,' ',70('-') )"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
#if (defined _MEMTRACE)
    if(my_task == 0 ) then
       lbnum=1
       call memmon_dump_fort('memmon.out',SubName//':start::',lbnum) 
    endif
#endif

   !  call timer_start(timer_total) ! start up the main timer

   !--- reset shr logging to my log file ---
   call shr_file_getLogUnit (shrlogunit) ! save log unit
   call shr_file_getLogLevel(shrloglev)  ! save log level
   call shr_file_setLogUnit (o_logunit)     ! set shared log unit to ocpl's logUnit
   if (seq_comm_iamroot(OCNID_o)) write(o_logunit,F00) "ENTER" 

   call seq_timemgr_EClockGetData( EClock, curr_ymd=ymd, curr_tod=tod )
   if (seq_comm_iamroot(OCNID_o)) write(o_logunit,*) subname,'cpl target model date: ymd, tod =',ymd,tod

   !----------------------------------------------------------------------------
   ! export data from roms (for pop 3d restoring)
   !----------------------------------------------------------------------------
   if (seq_comm_iamroot(OCNID_o)) write(o_logunit,F01) "export 3d fields from roms (for pop 3d restoring)" 
   call ocpl_roms_export( )

   !----------------------------------------------------------------------------
   ! map: pop -> roms  (pop 3d restoring)
   !----------------------------------------------------------------------------
   if (seq_comm_iamroot(OCNID_o)) write(o_logunit,F01) "map: roms->pop (for pop 3d restoring)" 
   call ocpl_map_roms2pop()

   !----------------------------------------------------------------------------
   ! import data into pop (pop 3d restoring)
   !----------------------------------------------------------------------------
   if (seq_comm_iamroot(OCNID_o)) write(o_logunit,F01) "import ocean coupling fields into pop (3d restoring)" 
!  call ocpl_pop_import( )
   call ocpl_pop_import(p2x_p ) ! adds non-standard fields to validate/debug pop restoring

   !----------------------------------------------------------------------------
   ! map: ocpl -> pop
   !----------------------------------------------------------------------------
   if (seq_comm_iamroot(OCNID_o)) write(o_logunit,'(2a)') subname,"map: x2o_o -> x2o_p"
   call mct_sMat_avMult(x2o_o, sMatp_o2p, x2o_p,vector=usevector)

   !----------------------------------------------------------------------------
   ! run pop
   !----------------------------------------------------------------------------
   if (seq_comm_iamroot(OCNID_o)) write(o_logunit,F01) "call pop_run_mct"
   call pop_run_mct( EClock, cdata_p, x2o_p, p2x_p)

   !----------------------------------------------------------------------------
   ! export data from pop (roms lateral BCs)
   !----------------------------------------------------------------------------
   if (seq_comm_iamroot(OCNID_o)) write(o_logunit,F01) "export ocean coupling fields from pop (roms lateral BCs)" 
   call ocpl_pop_export( p2x_2d_p, p2x_3d_p)    ! to do?  remove args as aVects are in ocpl_data_mod

   !----------------------------------------------------------------------------
   ! map: pop -> roms  (roms lateral BCs)
   !----------------------------------------------------------------------------
   if (seq_comm_iamroot(OCNID_o)) write(o_logunit,F01) "map: pop->roms (roms lateral BCs)" 
   call ocpl_map_pop2roms()

   !----------------------------------------------------------------------------
   ! import data to roms (roms lateral BCs)
   !----------------------------------------------------------------------------
   if (seq_comm_iamroot(OCNID_o)) write(o_logunit,F01) "import ocean coupling fields into roms (lateral BCs)" 
   call ocpl_roms_import( )

   !----------------------------------------------------------------------------
   ! run roms
   !----------------------------------------------------------------------------
   if (seq_comm_iamroot(OCNID_o)) write(o_logunit,*) subname,"map: x2o_o -> x2o_r"
   call mct_sMat_avMult(x2o_o, sMatp_o2r, x2o_r,vector=usevector)

   if (debug > 0) then
      k = mct_avect_indexra(x2o_o,'Foxx_lwup')
      if (seq_comm_iamroot(OCNID_o)) write(o_logunit,'(2a,2i6)'   ) subname,'<DEBUG> Foxx_lwup index = ',k
      if (seq_comm_iamroot(OCNID_o)) write(o_logunit,'(2a,2e12.4)') subname,'<DEBUG> min/max x2o_o Foxx_lwup= ',minval(x2o_o%rAttr(k,:)),maxval(x2o_o%rAttr(k,:))
      if (seq_comm_iamroot(OCNID_o)) write(o_logunit,'(2a,2e12.4)') subname,'<DEBUG> min/max x2o_r Foxx_lwup= ',minval(x2o_r%rAttr(k,:)),maxval(x2o_r%rAttr(k,:))
   end if

   if (seq_comm_iamroot(OCNID_o)) write(o_logunit,F01) "call roms_run_mct"
   call roms_run_mct( EClock, cdata_r, x2o_r, r2x_r)


   !--------------------------------------------------------------------------------------
   ! map & merge pop & roms output to create ocpl output
   !--------------------------------------------------------------------------------------
   if (seq_comm_iamroot(OCNID_o)) write(o_logunit,'(2a)') subname,"map: r2x_r -> r2x_o"
   call mct_sMat_avMult(r2x_r, sMatp_r2o, r2x_o,vector=usevector)

   if (seq_comm_iamroot(OCNID_o)) write(o_logunit,'(2a)') subname,"map: p2x_p -> p2x_o"
   call mct_sMat_avMult(p2x_p, sMatp_p2o, p2x_o,vector=usevector)

   if (debug > 0) then
      k = mct_avect_indexra(r2x_o,'So_t')
      write(o_logunit,'(2a,2e12.4)') subname,'<DEBUG> min/max r2x_r SST = ',minval(r2x_r%rAttr(k,:)),maxval(r2x_r%rAttr(k,:))
      write(o_logunit,'(2a,2e12.4)') subname,'<DEBUG> min/max r2x_o SST = ',minval(r2x_o%rAttr(k,:)),maxval(r2x_o%rAttr(k,:))
      write(o_logunit,'(2a,2e12.4)') subname,'<DEBUG> min/max p2x_p SST = ',minval(p2x_p%rAttr(k,:)),maxval(p2x_p%rAttr(k,:))
      write(o_logunit,'(2a,2e12.4)') subname,'<DEBUG> min/max o2x_o SST = ',minval(p2x_o%rAttr(k,:)),maxval(p2x_o%rAttr(k,:))
   end if

   if (seq_comm_iamroot(OCNID_o)) write(o_logunit,'(a)') subname,"merge: r2x_o -> o2x_o"
   k = mct_avect_indexra(o2x_o,'So_t')
   lsize_o = mct_aVect_lsize( o2x_o )
   do n=1,lsize_o
      o2x_o%rAttr(:,n) = p2x_o%rAttr(:,n)  ! start with all pop data, replace some with roms data
      if (r2x_o%rAttr(k,n) > 1.0) then  ! has data mapped from roms (unmapped cells have sst = 0)
          o2x_o%rAttr(:,n) = r2x_o%rAttr(:,n)  ! merge all fields
      !   o2x_o%rAttr(k,n) = r2x_o%rAttr(k,n)  ! merge SST only
      end if
   end do

   if (debug > 0 ) then
      k = mct_avect_indexra(o2x_o,'So_t')
      write(o_logunit,'(2a,2e12.4)') subname,'<DEBUG> min/max o2x_o SST+ = ',minval(o2x_o%rAttr(k,:)),maxval(o2x_o%rAttr(k,:))
   end if

   !--------------------------------------------------------------------------------------
   ! done
   !--------------------------------------------------------------------------------------
!  call timer_stop(timer_total)

#if (defined _MEMTRACE)
    if(my_task == 0) then
       lbnum=1
       call memmon_dump_fort('memmon.out',SubName//':end::',lbnum) 
       call memmon_reset_addr()
    endif
#endif

   !--- Reset shr logging to original values ---
   if (seq_comm_iamroot(OCNID_o)) write(o_logunit,F00) "EXIT " 
   call shr_file_setLogUnit (shrlogunit) ! restore log unit
   call shr_file_setLogLevel(shrloglev)  ! restore log level

end subroutine ocn_run_mct

!================================================================================
!================================================================================
!BOP
!
! !IROUTINE: ocn_final_mct
!
! !INTERFACE:

subroutine ocn_final_mct( EClock, cdata_o, x2o_o, o2x_o)

!
! !DESCRIPTION:
! Finalize POP
!
! !USES:

! !ARGUMENTS:
    type(ESMF_Clock)            , intent(inout) :: EClock
    type(seq_cdata)             , intent(inout) :: cdata_o
    type(mct_aVect)             , intent(inout) :: x2o_o
    type(mct_aVect)             , intent(inout) :: o2x_o
!
! !REVISION HISTORY:
! Author: Brian Kauffman
!EOP
!BOC
   !-----  local variables -----

   character(*), parameter :: subName = "(ocn_final_mct) "
   character(*), parameter :: F00 =  "( '(ocn_final_mct) ===== ',a,' ',70('=') )"
   character(*), parameter :: F01 =  "( '(ocn_final_mct) ----- ',a,' ',70('-') )"

!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------

   !--- set shr logging to ocpl's unit ---
   call shr_file_getLogUnit (shrlogunit) ! save log unit
   call shr_file_getLogLevel(shrloglev)  ! save log level
   call shr_file_setLogUnit (o_logunit)  ! set shared log unit to ocpl's logUnit
   if (seq_comm_iamroot(OCNID_o)) write(o_logunit,F00) "ENTER" 

   if (seq_comm_iamroot(OCNID_o)) write(o_logunit,F01) "pop_final_mct call" 
   call pop_final_mct( EClock, cdata_o, x2o_o, o2x_o)

   if (seq_comm_iamroot(OCNID_o)) write(o_logunit,F01) "roms_final_mct call" 
   call roms_final_mct( EClock, cdata_r, x2o_r, r2x_r)

   if (seq_comm_iamroot(OCNID_o)) write(o_logunit,F00) "EXIT" 

   !--- restore shr logging to original values ---
   call shr_file_setLogUnit (shrlogunit) ! restore log unit
   call shr_file_setLogLevel(shrloglev)  ! restore log level

end subroutine ocn_final_mct

!================================================================================
!================================================================================

end module ocn_comp_mct
