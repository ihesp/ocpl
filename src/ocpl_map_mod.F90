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

   integer,parameter,private :: debug = 1   ! debug level

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

   if (do_Scurtain) then
      k = k_Scurtain
      call shr_mct_queryConfigFile(mpicom_o,configFile, &
           "pop2roms_Scurtain_file:",mapfile,"pop2roms_Scurtain_type:",maptype)
      write(o_logunit,'(4a)') subName, "file = ",trim(mapfile)
      call shr_mct_sMatPInitnc(sMatp_p2rc(k),gsMap_o,gsMap_rc(k),trim(mapfile),trim(maptype),mpicom_o)
   else
      write(o_logunit,*) subname,"do_Scurtain = false" ;  call shr_sys_flush(o_logunit)
   end if

   if (do_Ecurtain) then
      k = k_Ecurtain
      call shr_mct_queryConfigFile(mpicom_o,configFile, &
           "pop2roms_Ecurtain_file:",mapfile,"pop2roms_Ecurtain_type:",maptype)
      write(o_logunit,'(4a)') subName, "file = ",trim(mapfile)
      call shr_mct_sMatPInitnc(sMatp_p2rc(k),gsMap_o,gsMap_rc(k),trim(mapfile),trim(maptype),mpicom_o)
   else
      write(o_logunit,*) subname,"do_Ecurtain = false" ;  call shr_sys_flush(o_logunit)
   end if

   if (do_Ncurtain) then
      k = k_Ncurtain
      call shr_mct_queryConfigFile(mpicom_o,configFile, &
           "pop2roms_Ncurtain_file:",mapfile,"pop2roms_Ncurtain_type:",maptype)
      write(o_logunit,'(4a)') subName, "file = ",trim(mapfile)
      call shr_mct_sMatPInitnc(sMatp_p2rc(k),gsMap_o,gsMap_rc(k),trim(mapfile),trim(maptype),mpicom_o)
   else
      write(o_logunit,*) subname,"do_Ncurtain = false" ;  call shr_sys_flush(o_logunit)
   end if

   if (do_Wcurtain) then
      k = k_Wcurtain
      call shr_mct_queryConfigFile(mpicom_o,configFile, &
           "pop2roms_Wcurtain_file:",mapfile,"pop2roms_Wcurtain_type:",maptype)
      write(o_logunit,'(4a)') subName, "file = ",trim(mapfile)
      call shr_mct_sMatPInitnc(sMatp_p2rc(k),gsMap_o,gsMap_rc(k),trim(mapfile),trim(maptype),mpicom_o)
   else
      write(o_logunit,*) subname,"Wcurtain = false" ;  call shr_sys_flush(o_logunit)
   end if

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

   use mod_param   , only: BOUNDS,             & ! roms internal data
                           globalISize0 => Lm, &
                           globalJSize0 => Mm  
   use mod_grid    , only: ROMS_GRID => GRID     ! roms internal data
   use mod_parallel, only: MyRank                ! roms internal data
   use grid        ,     only: pop_depth => zw   ! pop  internal data
 
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
   real(R8)    :: depth_r             ! depth of roms point being interpolated to
   real(R8),allocatable :: depth_p(:) ! depths of pop data 
   real(R8)    :: zr,z1p,z2p          ! depth roms value and surrounding pop layers
   logical     :: do_mapping          ! flags that this local tile reqires mapping
   logical     :: first_call = .true. ! flags 1st-time setup operations

   type(mct_aVect), allocatable :: p2x_3d_pvert_rc(:,:) ! aVect on roms horiz grid but pop vertical grid

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

      if (debug>0) write(*,*) subName,"<DEBUG> 2d map for curtain ",k,", do_mapping =",do_mapping
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
   ! in the vertical column containing the target roms value location...
   ! np1 is the the pop layer immediately *above* the desired roms value
   ! np2 is the the pop layer immediately *below* the desired roms value, np2 = np + 1
   ! w1 is the weight assigned to the pop value at layer np1 and is in the range [0,1]
   ! w2 is the weight assigned to np2, w2 = (1 - w1)
   ! roms_value = w1*pop_value(np1)  + w2*pop_value(np2) 
   !--------------------------------------------------------------------------------------

   if (first_call) then ! one-time initializations 
      !--- pop depth array, for layer mid-point, in meters
      allocate(depth_p(nlev_p))
      depth_p(1) =     0.01_R8*(pop_depth(1)  + 0.0_R8       )/2.0_R8
      do np=2,nlev_p 
         depth_p(np) = 0.01_R8*(pop_depth(np)+pop_depth(np-1))/2.0_R8 
      end do
      !--- work array pop data aVect, on roms grid, but at pop vertical levels
      allocate(p2x_3d_pvert_rc(4,nlev_p))
      do k=1,4        ! four curtains: N,E,S,W
      do np=1,nlev_p  ! all pop layers
         lsize = mct_aVect_lsize  (p2x_3d_rc(k,1))
         call mct_aVect_init(p2x_3d_pvert_rc(k,np),rlist=ocpl_fields_p2x_3d_fields,lsize=lsize)
         call mct_aVect_zero(p2x_3d_pvert_rc(k,np))
      end do
      end do
   end if

   do k=1,4 ! four curtains: N,E,S,W

      do_mapping = .false.
      if (k==k_Scurtain .and. do_Scurtain==.true.) do_mapping = .true.
      if (k==k_Ecurtain .and. do_Ecurtain==.true.) do_mapping = .true.
      if (k==k_Ncurtain .and. do_Ncurtain==.true.) do_mapping = .true.
      if (k==k_Wcurtain .and. do_Wcurtain==.true.) do_mapping = .true.

      if (debug>0) write(*,*) subName,"<DEBUG> 3d map for curtain ",k,", do_mapping =",do_mapping

      if (do_mapping) then

         !--- horizontal mapping, still on pop vertical levels --------------------------
         if (debug==1) write(*,'(2a,i4)') subName,"<DEBUG> horiz mapping for k =",k
         do np=1,nlev_p
            if (debug>1) write(*,'(2a,2i4)') subName,"<DEBUG> horiz mapping for k,np =",k,np
            call mct_sMat_avMult(p2x_3d_p(np),sMatp_p2rc(k),p2x_3d_pvert_rc(k,np))
         end do

         !--- vertical interpolation to roms layers -------------------------------------
         kfld = mct_aVect_indexRA(p2x_3d_rc(k,1),"So_temp" )
         lsize = mct_aVect_lsize (p2x_3d_rc(k,1))

         if (lsize < 1 .and. debug>0) write(*,*) subName,"<DEBUG> lsize=0 => NO interp for k = ",k
         if (lsize > 0) then
         do nr=1,nlev_r ! for each roms level
            call mct_aVect_zero(p2x_3d_rc(k,nr))
            if (debug>2) write(*,'(2a,2i4)') subName,"<DEBUG> vert interp for k,nr = ",k,nr

            do ij=1,lsize ! for each grid cell

                !--- determine depth of roms target value 
                if (k == k_Wcurtain) then        ! 
                   i = 1
                   j = ij + BOUNDS(nestID)%JstrR(MyRank) - 1
                else if (k == k_Ncurtain) then
                   i = ij + BOUNDS(nestID)%IstrR(MyRank) - 1
                   j = globalJSize0(nestID) + 1
                else if (k == k_Ecurtain) then
                   i = globalISize0(nestID) + 1
                   j = ij + BOUNDS(nestID)%JstrR(MyRank) - 1
                else if (k == k_Scurtain) then
                   i = ij + BOUNDS(nestID)%IstrR(MyRank) - 1
                   j = 1
                endif
                depth_r = -ROMS_GRID(nestID)%z_r(i,j,nr) ! depth in meters, positive down

                !--- determine interpolation weights of pop values surrounding roms target value
                w1  = -1.0  ! an invalid value
                w2  = -1.0   
                np2 = -1    ! an invalid value
                np1 = -1     
                if (     depth_r >= depth_p(nlev_p)) then ! roms deeper than pop 
                   np1 = nlev_p - 1
                   np2 = nlev_p
                   w1  = 0.0_R8   
                   w2  = 1.0_R8   ! all weight to lowest pop level
                else if (depth_r <= depth_p(1)     ) then ! roms shallower than pop
                   np1 = 1
                   np2 = 2
                   w1  = 1.0_R8  ! all weight to shallowest pop level
                   w2  = 0.0_R8  
                else                         ! find pop levels that surround roms point
                   do np = 1,nlev_p-1
                      if (depth_p(np) <= depth_r .and. depth_r <= depth_p(np+1)) then
                         np1 = np        ! the shallower layer
                         np2 = np1 + 1   ! the deeper    layer
                         w2 = ( depth_p(np1) - depth_r     ) &
                            / ( depth_p(np1) - depth_p(np2))
                         w1  = 1.0_R8 - w2  
                         exit
                      endif
                   enddo
                endif
                if (np1.eq.-1) write(*,*) 'ERROR: failed to find bounding layers'

                !--- apply interpolation weights
                if (debug>0 .and. (nr==nlev_r .or. mod(nr-1,10)==0) .and. mod(ij,10)==0) then
                   zr  = depth_r
                   z1p = depth_p(np1) 
                   z2p = depth_p(np2) 
                   write(*,'(2a,5i4,2f5.2,3f9.2)') subName,&
                      "<DEBUG> k,nr,ij,np1,np2,w1,w2,z1,z,z2=",k,nr,ij,np1,np2,w1,w2,z1p,zr,z2p
                end if
                p2x_3d_rc(k,nr)%rAttr(:,ij) = w1*p2x_3d_pvert_rc(k,np1)%rAttr(:,ij) &
                                            + w2*p2x_3d_pvert_rc(k,np2)%rAttr(:,ij)

            end do ! do ij - over grid cell
         end do ! do nr ~ roms layer
         end if ! lsize > 0 
      end if ! curtain is active
   end do  ! do k=1,4  ~ curtains: N,E,S,W

   first_call = .false.
   write(o_logunit,*) subname,"Exit" ;  call shr_sys_flush(o_logunit)

end subroutine ocpl_map_pop2roms

!=========================================================================================
!=========================================================================================
end module ocpl_map_mod
