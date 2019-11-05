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

   integer,parameter :: debug = 0   ! debug level

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

   if (do_Wcurtain) then
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

   use mod_param,only: BOUNDS,ROMS_imax => Lm,ROMS_jmax => Mm  ! roms internal data
   use mod_grid, only: ROMS_GRID => GRID   ! roms internal data
   use mod_parallel, only: MyRank          ! roms internal data
   use grid,     only: pop_depth => zw     ! pop  internal data
 
   integer(IN),parameter  :: nestID = 1    ! roms nest (grid/domain) #1 


!EOP
!BOC

!-----  local variables --------------------------------------------------------

!  use domain_size, only: km  ! # vertical levels, for 3d coupling

   integer(IN) :: i,j,ij              ! 1d & 2d array indicies
   integer(IN) :: k                   ! curtain index: N,E,S,W
   integer(IN) :: lsize               ! local tile size
   integer(IN) :: np,nr               ! level index for pop & roms
   integer(IN) :: np1,np2             ! lower,upper pop levels used in vertical interpolation
   real(R8)    :: w1 ,w2              ! interp weight assigned to lower,upper pop levels
   real(R8)    :: roms_sigma          ! depth of roms point being interpolated to
   logical     :: do_mapping          ! flags that this local tile reqires mapping
   logical     :: first_call = .true. ! flags 1st-time setup operations

   type(mct_aVect), pointer :: p2x_3d_pvert_rc(:,:) ! work aVect, data on roms curtain but pop vertical grid

   integer(IN) :: kfld                ! field index
   real(R8)    :: tmin,tmax           ! debug info

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
      if (k==k_Scurtain .and. do_Scurtain==.true.) do_mapping = .true.
      if (k==k_Ecurtain .and. do_Ecurtain==.true.) do_mapping = .true.
      if (k==k_Ncurtain .and. do_Ncurtain==.true.) do_mapping = .true.
      if (k==k_Wcurtain .and. do_Wcurtain==.true.) do_mapping = .true.

     
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
   ! 
   ! in the vertical column containing the target roms value Fr ...
   ! Fp1 is the the pop value immediately *below* the desired roms value
   ! Fp2 is the the pop value immediately *above* the desired roms value
   ! np1 is the pop level of Fp1
   ! np2 is the pop level of Fp2, pn2 = np1 - 1
   ! w1 is the weight assigned to Fp1 and is in the range [0,1]
   ! w2 is the weight assigned to Fp2, w2 = (1 - w1)
   ! Fr = w1*Fp1 + w2*Fp2
   !--------------------------------------------------------------------------------------
   if (first_call) then
      allocate(p2x_3d_pvert_rc(4,nlev_p))
   end if
   do k=1,4 ! four curtains: N,E,S,W

      do_mapping = .false.
      if (k==k_Scurtain .and. do_Scurtain==.true.) do_mapping = .true.
      if (k==k_Ecurtain .and. do_Ecurtain==.true.) do_mapping = .true.
      if (k==k_Ncurtain .and. do_Ncurtain==.true.) do_mapping = .true.
      if (k==k_Wcurtain .and. do_Wcurtain==.true.) do_mapping = .true.

      if (debug>0) write(*,*) subName,"3d map for curtain ",k,", do_mapping =",do_mapping

      !f (do_mapping) call mct_sMat_avMult(p2x_2d_p,sMatp_p2rc(k),p2x_2d_rc(k),rList=ocpl_fields_p2x_2d_fields)
      if (do_mapping) then

         if (debug>0) write(*,*) subName,"horizontal mapping"
         !--- horizontal mapping, still on pop vertical levels --------------------------
         do np=1,nlev_p
            if (first_call) then
               lsize = mct_aVect_lsize  (p2x_3d_rc(k,1))
               call mct_aVect_init(p2x_3d_pvert_rc(k,np),rlist=ocpl_fields_p2x_3d_fields,lsize=lsize)
               call mct_aVect_zero(p2x_3d_pvert_rc(k,np))
            end if
            call mct_sMat_avMult(p2x_3d_p(np),sMatp_p2rc(k),p2x_3d_pvert_rc(k,np))
         end do

         !--- vertical interpolation to roms layers -------------------------------------
         kfld = mct_aVect_indexRA(p2x_3d_rc(k,1),"So_temp" )
         lsize = mct_aVect_lsize (p2x_3d_rc(k,kfld))

         if (lsize < 1) then
            if (debug>0) write(*,*) subName,"NO vertical interp, lsize = 0"
         else
            if (debug>0) write(*,*) subName,"vertical interpolation..."
         do nr=1,nlev_r ! for each roms level
            call mct_aVect_zero(p2x_3d_rc(k,nr))
            if (debug>0) write(*,*) subName,"working on roms level = ",nr

            do ij=1,lsize ! for each grid cell
                if (debug>0) write(*,*) subName,"   working on cell ij = ",ij

                !--- determine depth of roms target value 
                if (k == k_Wcurtain) then        ! 
                   i = 1
                   j = ij + BOUNDS(nestID)%JstrR(MyRank) - 1
                else if (k == k_Ncurtain) then
                   i = ij + BOUNDS(nestID)%IstrR(MyRank) - 1
                   j = ROMS_jmax(nestID) + 1
                else if (k == k_Ecurtain) then
                   i = ROMS_imax(nestID) + 1
                   j = ij + BOUNDS(nestID)%JstrR(MyRank) - 1
                else if (k == k_Scurtain) then
                   i = ij + BOUNDS(nestID)%IstrR(MyRank) - 1
                   j = 1
                endif
                roms_sigma = -ROMS_GRID(nestID)%z_r(i,j,nr) * 100.0d0 ! convert: m -> cm

                if (debug>0) write(*,*) subName,"   determine weights"
                !--- determine interpolation weights of pop values surrounding roms target value
                w1  = 0.0_R8
                np1 = -1     ! an invalid value
                if (     roms_sigma >= pop_depth(nlev_p)) then ! roms deeper than pop 
                   w1  = 1.0_R8   ! all weight is given to lowest pop level
                   np1 = nlev_p
                else if (roms_sigma <= pop_depth(1)    ) then ! roms shallower than pop
                   w1  = 0.0_R8  ! all weight is given to shallowest pop level
                   np1 = 2
                else                 
                   do np = nlev_p,2,-1   ! find pop levels that surround roms level
                      if (pop_depth(np  ) >= roms_sigma .and. &
                          pop_depth(np-1) <  roms_sigma       ) then
                         w1 = ( pop_depth(np-1) - roms_sigma   ) &
                            / ( pop_depth(np-1) - pop_depth(np))
                         np1 = np
                         exit
                      endif
                   enddo
                endif
                w2  = 1.0_R8 - w1
                np2 = np1 + 1
                if (np1.eq.-1) write(*,*) 'FATAL ERROR in vertical interpolation'

                !--- apply interpolation weights
                if (debug>0) write(*,*) subName,"   apply weights"
                p2x_3d_rc(k,nr)%rAttr(:,ij) = w1*p2x_3d_pvert_rc(k,np1)%rAttr(:,ij) &
                                            + w2*p2x_3d_pvert_rc(k,np2)%rAttr(:,ij)

            end do ! grid cell
         end do ! level
         end if ! lsize > 0
      end if ! curtain is active
   end do  ! curtain: N,E,S,W

   first_call = .false.
   write(o_logunit,*) subname,"Exit" ;  call shr_sys_flush(o_logunit)

end subroutine ocpl_map_pop2roms

!=========================================================================================
!=========================================================================================
end module ocpl_map_mod
