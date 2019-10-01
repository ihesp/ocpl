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
   use io_types     , only :  pop_stdout => stdout  ! ocpl redirects stdout for pop
   use mod_iounits  , only : roms_stdout => stdout  ! ocpl redirects stdout for roms

   use mct_mod
   use shr_mct_mod
   use esmf

   use seq_flds_mod
   use seq_cdata_mod
   use seq_infodata_mod
   use seq_timemgr_mod
   use seq_comm_mct,      only : seq_comm_suffix, seq_comm_inst, seq_comm_name

   use shr_file_mod 
   use shr_cal_mod,       only : shr_cal_date2ymd
   use shr_sys_mod
   use shr_kind_mod, only: IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL

   use ocpl_data_mod
   use ocpl_pop_mod

! !PUBLIC MEMBER FUNCTIONS:

  implicit none
  private              ! default: all data is private
  SAVE                 ! save everything

  public :: ocn_init_mct
  public :: ocn_run_mct
  public :: ocn_final_mct


!
! ! PUBLIC DATA:
!
! !REVISION HISTORY:
! Author: Brian Kauffman
!
!EOP
! !PRIVATE MODULE FUNCTIONS:

   integer(IN) ::  stdout = 6             ! move to ocpl_data_mod
   integer(IN) ::  shrlogunit,shrloglev   ! to save & restore shared log units & levels

   integer (IN), private :: cpl_write_restart  ! flag id for write restart
   integer (IN), private :: cpl_write_history  ! flag id for write history
   integer (IN), private :: cpl_write_tavg     ! flag id for write tavg      
   integer (IN), private :: cpl_diag_global    ! flag id for computing diagnostics
   integer (IN), private :: cpl_diag_transp    ! flag id for computing diagnostics

   integer (IN)  :: nsend, nrecv

!===============================================================================
contains
!===============================================================================

!BOP ===========================================================================
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
! !INTERFACE: ------------------------------------------------------------------

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
   integer(IN) :: ni_r, nj_r
   integer(IN) :: rc
   integer(IN) :: k       ! aVect field index 
   integer(IN) :: n       ! aVect cell index 
   integer(IN) :: lSize_o ! aVect local size of *_o grid

   logical :: restart

   character(CL) :: start_type
   character(CS) :: case_name

   integer(IN) :: ymd,tod

   character(*), parameter :: subName = "(ocn_init_mct) "
   character(*), parameter :: F00 =  "( '(ocn_init_mct) ===== ',a,' ',70('=') )"
   character(*), parameter :: F01 =  "( '(ocn_init_mct) ----- ',a,' ',70('-') )"

!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------

   call shr_file_getLogUnit (shrlogunit) ! save log unit
   call shr_file_getLogLevel(shrloglev)  ! save log level
   stdout = shr_file_getUnit()           ! get/reserve unused unit number for ocpl's stdout
   call shr_file_setIO("ocn_modelio.nml",stdout)  ! assign ocpl stdout to log file
   pop_stdout  = stdout                  ! set pop    stdout = ocpl's stdout
   roms_stdout = stdout                  ! set roms   stdout = ocpl's stdout
   call shr_file_setLogUnit (stdout)     ! set shared stdout = ocpl's stdout

   write(*,F00) " ENTER"

   !--- Get data/pointers out of cdata ---
   call seq_cdata_setptrs(cdata_o, ID=OCNID_o, dom=dom_o, gsMap=gsMap_o, infodata=infodata_o, mpicom=mpicom_o,name=name_o)
   write(stdout,*) subName, 'cdata_o ID     : ' ,  OCNID_o
   write(stdout,*) subName, 'cdata_o mpicom : ' ,  mpicom_o
   write(stdout,*) subName, 'cdata_o name   : ' // trim(name_o )

   !--- Get data out of infodata ---
   call seq_infodata_GetData(infodata_o, case_name=case_name, start_type=start_type)
   write(stdout,*) subName, 'case_name  : ' // trim(case_name )
   write(stdout,*) subName, 'start_type : ' // trim(start_type)

   !----------------------------------------------------------------------------
   ! pop setup
   !----------------------------------------------------------------------------

   !--- initialize pop --
   write(stdout,F01) "pop_init_mct call" ; call shr_sys_flush(stdout)
   call pop_init_mct( EClock, cdata_o, x2o_o, o2x_o, NLFilename )
   write(stdout,F01) "pop_init_mct return" ; call shr_sys_flush(stdout)
   call seq_timemgr_EClockGetData( EClock, curr_ymd=ymd, curr_tod=tod )
   write(stdout,*) subname,'ymd, tod =',ymd,tod

   call seq_infodata_GetData(infodata_o, ocn_nx = ni_o, ocn_ny = nj_o)

   !----------------------------------------------------------------------------
   ! roms setup
   !----------------------------------------------------------------------------

   !--- set data, and assign pointers, inside cdata_r --- pop & roms share ID,mpicom,infdodata
   call seq_cdata_init(cdata_r, ID=OCNID_o, dom=dom_r, gsMap=gsMap_r, infodata=infodata_o, name='ROMS')

   !--- initialize roms ---
   write(stdout,F01) "roms_init_mct call" ; call shr_sys_flush(stdout)
   call roms_init_mct( EClock, cdata_r, x2o_r, r2x_r)  ! , NLFilename )
   write(stdout,F01) "roms_init_mct return" ; call shr_sys_flush(stdout)
   call seq_timemgr_EClockGetData( EClock, curr_ymd=ymd, curr_tod=tod )
   write(stdout,*) subname,'ymd, tod =',ymd,tod

   !--- put pop dims into infodata, were changed to regional dims by roms_comp_init() ---
   call seq_infodata_GetData(infodata_o, ocn_nx = ni_r, ocn_ny = nj_r)
   call seq_infodata_PutData(infodata_o, ocn_nx = ni_o, ocn_ny = nj_o)

   !----------------------------------------------------------------------------
   ! sanity check on some data
   !----------------------------------------------------------------------------
   call seq_cdata_setptrs(cdata_o, ID=OCNID_o,                    mpicom=mpicom_o,name=name_o)
   write(stdout,*) subName, 'cdata_o ID     : ' ,  OCNID_o
   write(stdout,*) subName, 'cdata_o mpicom : ' ,  mpicom_o
   write(stdout,*) subName, 'cdata_o name   : ' // trim(name_o )
   write(stdout,*) subName, 'ni_o,nj_o      : ' ,  ni_o,nj_o
   call seq_cdata_setptrs(cdata_r, ID=OCNID_r,                    mpicom=mpicom_r,name=name_r)
   write(stdout,*) subName, 'cdata_r ID     : ' ,  OCNID_r
   write(stdout,*) subName, 'cdata_r mpicom : ' ,  mpicom_r
   write(stdout,*) subName, 'cdata_r name   : ' // trim(name_r )
   write(stdout,*) subName, 'ni_r,nj_r      : ' ,  ni_r,nj_r

   !----------------------------------------------------------------------------
   ! initialize maps
   !----------------------------------------------------------------------------
   write(stdout,*) subName, "initialize r2o map..."
   call shr_mct_sMatPInitnc(sMatp_r2o,gsMap_r,gsMap_o, trim(r2o_mapfile),trim(r2o_maptype),mpicom_o)

   write(stdout,*) subName, "initialize o2r map..."
   call shr_mct_sMatPInitnc(sMatp_o2r,gsMap_o,gsMap_r, trim(o2r_mapfile),trim(o2r_maptype),mpicom_o)

   lsize_o = mct_aVect_lsize(o2x_o)
   call mct_aVect_init(r2x_o, r2x_r, lsize=lsize_o)
   call mct_aVect_zero(r2x_o)
   write(stdout,'(2a,2i6)'   ) subname,'<DEBUG> r2x_o lsize, nflds = ', mct_aVect_lsize (r2x_o),mct_aVect_nRAttr(r2x_o)

   write(stdout,*) subName, "done initializing maps"

   !----------------------------------------------------------------------------
   ! merge roms & pop output
   !----------------------------------------------------------------------------
   write(stdout,*) subname,"map: r2x_r -> r2x_o"
   call mct_sMat_avMult(r2x_r, sMatp_r2o, r2x_o,vector=usevector)
   write(stdout,*) subname,' merge roms & pop output'
   k = mct_avect_indexra(o2x_o,'So_t')
   write(stdout,'(2a,2i6)'   ) subname,'<DEBUG> So_t index = ',k
   write(stdout,'(2a,2e12.4)') subname,'<DEBUG> min/max o2x_o SST0= ',minval(o2x_o%rAttr(k,:)),maxval(o2x_o%rAttr(k,:))
   write(*     ,'(2a,2e12.4)') subname,'<DEBUG> min/max o2x_o SST0= ',minval(o2x_o%rAttr(k,:)),maxval(o2x_o%rAttr(k,:))
   lsize_o = mct_aVect_lsize( o2x_o )
   do n=1,lsize_o
      if (r2x_o%rAttr(k,n) > 1.0) then  ! has data mapped from roms (unmapped cells have sst = 0)
          o2x_o%rAttr(:,n) = r2x_o%rAttr(:,n)  ! merge all fields
      !   o2x_o%rAttr(k,n) = r2x_o%rAttr(k,n)  ! merge SST only
      end if
   end do
   write(stdout,'(2a,2e12.4)') subname,'<DEBUG> min/max r2x_r SST = ',minval(r2x_r%rAttr(k,:)),maxval(r2x_r%rAttr(k,:))
   write(stdout,'(2a,2e12.4)') subname,'<DEBUG> min/max r2x_o SST = ',minval(r2x_o%rAttr(k,:)),maxval(r2x_o%rAttr(k,:))
   write(stdout,'(2a,2e12.4)') subname,'<DEBUG> min/max o2x_o SST = ',minval(o2x_o%rAttr(k,:)),maxval(o2x_o%rAttr(k,:))
   write(*     ,'(2a,2e12.4)') subname,'<DEBUG> min/max r2x_r SST = ',minval(r2x_r%rAttr(k,:)),maxval(r2x_r%rAttr(k,:))
   write(*     ,'(2a,2e12.4)') subname,'<DEBUG> min/max r2x_o SST = ',minval(r2x_o%rAttr(k,:)),maxval(r2x_o%rAttr(k,:))
   write(*     ,'(2a,2e12.4)') subname,'<DEBUG> min/max o2x_o SST = ',minval(o2x_o%rAttr(k,:)),maxval(o2x_o%rAttr(k,:))

   !----------------------------------------------------------------------------
   write(stdout,*) subName, "init data for 3D coupling with pop" ; call shr_sys_flush(stdout)
   !----------------------------------------------------------------------------
   call ocpl_pop_init( o2x_o, p2x_2d_p, p2x_3d_p)

   write(stdout,F00) "EXIT" ; call shr_sys_flush(stdout)

   !--- Reset shr logging to original values ---
   call shr_file_setLogUnit (shrlogunit) ! restore log unit
   call shr_file_setLogLevel(shrloglev)  ! restore log level

 end subroutine ocn_init_mct

!================================================================================
!================================================================================
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

!-----------------------------------------------------------------------
!
!  start up the main timer
!
!-----------------------------------------------------------------------

!  call timer_start(timer_total)

   !--- reset shr logging to my log file ---
   call shr_file_getLogUnit (shrlogunit) ! save log unit
   call shr_file_getLogLevel(shrloglev)  ! save log level
   call shr_file_setLogUnit (stdout)     ! set shared log unit to ocpl's stdout
   write(stdout,F00) "ENTER" ; call shr_sys_flush(stdout)

   call seq_timemgr_EClockGetData( EClock, curr_ymd=ymd, curr_tod=tod )
   write(stdout,*) subname,'cpl target model date: ymd, tod =',ymd,tod

   !----------------------------------------------------------------------------
   ! run pop
   !----------------------------------------------------------------------------
   write(stdout,F01) "call pop_run_mct"
   call pop_run_mct( EClock, cdata_o, x2o_o, o2x_o)

   !----------------------------------------------------------------------------
   write(stdout,F01) "extract ocean coupling fields from pop (roms lateral BCs" ; call shr_sys_flush(stdout)
   !----------------------------------------------------------------------------
   call ocpl_pop_export( p2x_2d_p, p2x_3d_p)

   !----------------------------------------------------------------------------
   ! run roms
   !----------------------------------------------------------------------------
   write(stdout,*) subname,"map: x2o_o -> x2o_r"
   call mct_sMat_avMult(x2o_o, sMatp_o2r, x2o_r,vector=usevector)
   k = mct_avect_indexra(x2o_o,'Foxx_lwup')
   write(stdout,'(2a,2i6)'   ) subname,'<DEBUG> Foxx_lwup index = ',k
   write(stdout,'(2a,2e12.4)') subname,'<DEBUG> min/max x2o_o Foxx_lwup= ',minval(x2o_o%rAttr(k,:)),maxval(x2o_o%rAttr(k,:))
   write(stdout,'(2a,2e12.4)') subname,'<DEBUG> min/max x2o_r Foxx_lwup= ',minval(x2o_r%rAttr(k,:)),maxval(x2o_r%rAttr(k,:))
   write(*     ,'(2a,2e12.4)') subname,'<DEBUG> min/max x2o_o Foxx_lwup= ',minval(x2o_o%rAttr(k,:)),maxval(x2o_o%rAttr(k,:))
   write(*     ,'(2a,2e12.4)') subname,'<DEBUG> min/max x2o_r Foxx_lwup= ',minval(x2o_r%rAttr(k,:)),maxval(x2o_r%rAttr(k,:))

   write(stdout,F01) "call roms_run_mct"
   call roms_run_mct( EClock, cdata_r, x2o_r, r2x_r)

   !----------------------------------------------------------------------------
   ! merge roms output into pop output
   !----------------------------------------------------------------------------
   write(stdout,*) subname,"map: r2x_r -> r2x_o"
   call mct_sMat_avMult(r2x_r, sMatp_r2o, r2x_o,vector=usevector)
   write(stdout,*) subname,' merge roms & pop output'
   k = mct_avect_indexra(o2x_o,'So_t')
   write(stdout,'(2a,2i6)'   ) subname,'<DEBUG> So_t index = ',k
   write(stdout,'(2a,2e12.4)') subname,'<DEBUG> min/max o2x_o SST0= ',minval(o2x_o%rAttr(k,:)),maxval(o2x_o%rAttr(k,:))
   write(*     ,'(2a,2e12.4)') subname,'<DEBUG> min/max o2x_o SST0= ',minval(o2x_o%rAttr(k,:)),maxval(o2x_o%rAttr(k,:))
   lsize_o = mct_aVect_lsize( o2x_o )
   do n=1,lsize_o
      if (r2x_o%rAttr(k,n) > 1.0) then  ! has data mapped from roms (unmapped cells have sst = 0)
          o2x_o%rAttr(:,n) = r2x_o%rAttr(:,n)  ! merge all fields
      !   o2x_o%rAttr(k,n) = r2x_o%rAttr(k,n)  ! merge SST only
      end if
   end do
   write(stdout,'(2a,2e12.4)') subname,'<DEBUG> min/max r2x_r SST = ',minval(r2x_r%rAttr(k,:)),maxval(r2x_r%rAttr(k,:))
   write(stdout,'(2a,2e12.4)') subname,'<DEBUG> min/max r2x_o SST = ',minval(r2x_o%rAttr(k,:)),maxval(r2x_o%rAttr(k,:))
   write(stdout,'(2a,2e12.4)') subname,'<DEBUG> min/max o2x_o SST = ',minval(o2x_o%rAttr(k,:)),maxval(o2x_o%rAttr(k,:))
   write(*     ,'(2a,2e12.4)') subname,'<DEBUG> min/max r2x_r SST = ',minval(r2x_r%rAttr(k,:)),maxval(r2x_r%rAttr(k,:))
   write(*     ,'(2a,2e12.4)') subname,'<DEBUG> min/max r2x_o SST = ',minval(r2x_o%rAttr(k,:)),maxval(r2x_o%rAttr(k,:))
   write(*     ,'(2a,2e12.4)') subname,'<DEBUG> min/max o2x_o SST = ',minval(o2x_o%rAttr(k,:)),maxval(o2x_o%rAttr(k,:))

!  call timer_stop(timer_total)

#if (defined _MEMTRACE)
    if(my_task == 0) then
       lbnum=1
       call memmon_dump_fort('memmon.out',SubName//':end::',lbnum) 
       call memmon_reset_addr()
    endif
#endif

   !--- Reset shr logging to original values ---
   write(stdout,F00) "EXIT " ; call shr_sys_flush(stdout)
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
   call shr_file_setLogUnit (stdout)     ! set shared log unit to ocpl's stdout
   write(stdout,F00) "ENTER" ; call shr_sys_flush(stdout)

   write(stdout,F01) "pop_final_mct call" ; call shr_sys_flush(stdout)
   call pop_final_mct( EClock, cdata_o, x2o_o, o2x_o)
   write(stdout,F01) "pop_final_mct return" ; call shr_sys_flush(stdout)
   write(stdout,F01) "roms_final_mct call" ; call shr_sys_flush(stdout)
   call roms_final_mct( EClock, cdata_r, x2o_r, r2x_r)
   write(stdout,F01) "roms_final_mct return" ; call shr_sys_flush(stdout)

   write(stdout,F00) "EXIT" ; call shr_sys_flush(stdout)

   !--- restore shr logging to original values ---
   call shr_file_setLogUnit (shrlogunit) ! restore log unit
   call shr_file_setLogLevel(shrloglev)  ! restore log level

end subroutine ocn_final_mct

!================================================================================
!================================================================================

end module ocn_comp_mct
