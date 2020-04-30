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
   private              ! all data private by default 
   SAVE                 ! save everything

! !PUBLIC MEMBER FUNCTIONS:

   public :: ocpl_map_init
   public :: ocpl_map_pop2roms
   public :: ocpl_map_roms2pop

! !PUBLIC DATA:

   !--- horizontal maps (2d) ---
   type(mct_sMatp),public :: sMatp_o2p    ! maps ocpl -> pop 
   type(mct_sMatp),public :: sMatp_p2o    ! maps pop  -> ocpl

   type(mct_sMatp),public :: sMatp_o2r    ! maps ocpl -> roms
   type(mct_sMatp),public :: sMatp_r2o    ! maps roms -> ocpl

   type(mct_sMatp),public :: sMatp_r2p    ! maps roms -> pop
   type(mct_sMatp),public :: sMatp_p2r    ! maps pop  -> roms

! !PRIVATE MODULE DATA:

!  integer(IN),parameter :: debug  = 0    ! debug level
   integer(IN)           :: debug  = 0    ! debug level
   integer(IN),parameter :: nestID = 1    ! roms nest (grid/domain) #1 

   type(mct_sMatp)       :: sMatp_p2rc(4) ! curtain amps: pop -> roms (4-curtains: S,E,N,W)

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

   write(o_logUnit,'(2a)') subname,"Enter" ;  call shr_sys_flush(o_logUnit)
   write(o_logUnit,'(2a,i2)') subname,"debug level = ",debug

   !--------------------------------------------------------------------------------------
   write(o_logUnit,*) subname,"Init pop<->roms surface and curtain maps"
   !--------------------------------------------------------------------------------------

   call shr_mct_queryConfigFile(mpicom_p,configFile,"roms2pop_file:",mapfile,"roms2pop_type:",maptype)
   write(o_logUnit,'(4a)') subName, "file = ",trim(mapfile)
   call shr_mct_sMatPInitnc(sMatp_r2p,gsMap_r,gsMap_p,trim(mapfile),trim(maptype),mpicom_p)

   call shr_mct_queryConfigFile(mpicom_p,configFile,"pop2roms_file:",mapfile,"pop2roms_type:",maptype)
   write(o_logUnit,'(4a)') subName, "file = ",trim(mapfile)
   call shr_mct_sMatPInitnc(sMatp_p2r,gsMap_p,gsMap_r,trim(mapfile),trim(maptype),mpicom_p)

   if (do_Scurtain) then
      k = k_Scurtain
      call shr_mct_queryConfigFile(mpicom_p,configFile, &
           "pop2roms_Scurtain_file:",mapfile,"pop2roms_Scurtain_type:",maptype)
      write(o_logUnit,'(4a)') subName, "file = ",trim(mapfile)
      call shr_mct_sMatPInitnc(sMatp_p2rc(k),gsMap_p,gsMap_rc(k),trim(mapfile),trim(maptype),mpicom_p)
   else
      write(o_logUnit,'(2a)') subname,"do_Scurtain = false" ;  call shr_sys_flush(o_logUnit)
   end if

   if (do_Ecurtain) then
      k = k_Ecurtain
      call shr_mct_queryConfigFile(mpicom_p,configFile, &
           "pop2roms_Ecurtain_file:",mapfile,"pop2roms_Ecurtain_type:",maptype)
      write(o_logUnit,'(4a)') subName, "file = ",trim(mapfile)
      call shr_mct_sMatPInitnc(sMatp_p2rc(k),gsMap_p,gsMap_rc(k),trim(mapfile),trim(maptype),mpicom_p)
   else
      write(o_logUnit,'(2a)') subname,"do_Ecurtain = false" ;  call shr_sys_flush(o_logUnit)
   end if

   if (do_Ncurtain) then
      k = k_Ncurtain
      call shr_mct_queryConfigFile(mpicom_p,configFile, &
           "pop2roms_Ncurtain_file:",mapfile,"pop2roms_Ncurtain_type:",maptype)
      write(o_logUnit,'(4a)') subName, "file = ",trim(mapfile)
      call shr_mct_sMatPInitnc(sMatp_p2rc(k),gsMap_p,gsMap_rc(k),trim(mapfile),trim(maptype),mpicom_p)
   else
      write(o_logUnit,'(2a)') subname,"do_Ncurtain = false" ;  call shr_sys_flush(o_logUnit)
   end if

   if (do_Wcurtain) then
      k = k_Wcurtain
      call shr_mct_queryConfigFile(mpicom_p,configFile, &
           "pop2roms_Wcurtain_file:",mapfile,"pop2roms_Wcurtain_type:",maptype)
      write(o_logUnit,'(4a)') subName, "file = ",trim(mapfile)
      call shr_mct_sMatPInitnc(sMatp_p2rc(k),gsMap_p,gsMap_rc(k),trim(mapfile),trim(maptype),mpicom_p)
   else
      write(o_logUnit,'(2a)') subname,"Wcurtain = false" ;  call shr_sys_flush(o_logUnit)
   end if

   !--------------------------------------------------------------------------------------
   write(o_logUnit,*) subname,"Init ocpl<->roms surface maps" 
   !--------------------------------------------------------------------------------------

   call shr_mct_queryConfigFile(mpicom_o,configFile,"ocpl2pop_file:",mapfile,"ocpl2pop_type:",maptype)
   write(o_logUnit,'(4a)') subName, "file = ",trim(mapfile)
   call shr_mct_sMatPInitnc(sMatp_o2p,gsMap_o,gsMap_p,trim(mapfile),trim(maptype),mpicom_o)

   call shr_mct_queryConfigFile(mpicom_o,configFile,"pop2ocpl_file:",mapfile,"pop2ocpl_type:",maptype)
   write(o_logUnit,'(4a)') subName, "file = ",trim(mapfile)
   call shr_mct_sMatPInitnc(sMatp_p2o,gsMap_p,gsMap_o,trim(mapfile),trim(maptype),mpicom_o)


   call shr_mct_queryConfigFile(mpicom_o,configFile,"roms2ocpl_file:",mapfile,"roms2ocpl_type:",maptype)
   write(o_logUnit,'(4a)') subName, "file = ",trim(mapfile)
   call shr_mct_sMatPInitnc(sMatp_r2o,gsMap_r,gsMap_o,trim(mapfile),trim(maptype),mpicom_o)

   call shr_mct_queryConfigFile(mpicom_o,configFile,"ocpl2roms_file:",mapfile,"ocpl2roms_type:",maptype)
   write(o_logUnit,'(4a)') subName, "file = ",trim(mapfile)
   call shr_mct_sMatPInitnc(sMatp_o2r,gsMap_o,gsMap_r,trim(mapfile),trim(maptype),mpicom_o)


   call shr_mct_queryConfigFile(mpicom_o,configFile,"roms2pop_file:",mapfile,"roms2pop_type:",maptype)
   write(o_logUnit,'(4a)') subName, "file = ",trim(mapfile)
   call shr_mct_sMatPInitnc(sMatp_r2p,gsMap_r,gsMap_p,trim(mapfile),trim(maptype),mpicom_o)

   call shr_mct_queryConfigFile(mpicom_o,configFile,"pop2roms_file:",mapfile,"pop2roms_type:",maptype)
   write(o_logUnit,'(4a)') subName, "file = ",trim(mapfile)
   call shr_mct_sMatPInitnc(sMatp_p2r,gsMap_p,gsMap_r,trim(mapfile),trim(maptype),mpicom_o)

   !--------------------------------------------------------------------------------------
   write(o_logUnit,*) subname,"Init ocpl<->pop  surface maps" 
   !--------------------------------------------------------------------------------------
   ! KLUDGE: temporarily assume ocpl has exactly pop domain & gsMap, "map" = copy 

   write(o_logUnit,'(2a)') subname,"Exit" ;  call shr_sys_flush(o_logUnit)

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
   use grid        , only: pop_depth => zt       ! pop  internal data
 
!  integer(IN),parameter  :: nestID = 1    ! roms nest (grid/domain) #1  


!EOP
!BOC

!-----  local variables --------------------------------------------------------

!  use domain_size, only: km  ! # vertical levels, for 3d coupling

   integer(IN) :: i,j,ij,n            ! 1d & 2d array indicies
   integer(IN) :: k                   ! curtain index: N,E,S,W
   integer(IN) :: lsize               ! local tile size
   integer(IN) :: np,nr               ! level index for pop & roms
   integer(IN) :: np1,np2             ! lower,upper pop levels used in vertical interpolation
   integer(IN) :: nFlds               ! number of fields in an aVect

   real(R8)    :: w1 ,w2              ! interp weight assigned to lower,upper pop levels
   real(R8)    :: F1 ,F2, Fr          ! pop data values and target roms data value
   real(R8)    :: depth_r             ! depth of roms point being interpolated to
   real(R8),allocatable :: depth_p(:) ! depths of pop data 
   real(R8)    :: zr,zp1,zp2          ! depth roms value and surrounding pop layers
   logical     :: do_mapping          ! flags that this local tile reqires mapping
   logical     :: first_call = .true. ! flags 1st-time setup operations

   type(mct_aVect),allocatable :: p2x_3d_pvert_rc(:,:) ! aVect on roms horiz grid but pop vertical grid

   integer(IN) :: kfld                ! field index
   real(R8)    :: xmin,xmax           ! debug info

   character(*), parameter :: subName = "(ocpl_map_pop2roms) "

   save

!-----------------------------------------------------------------------------------------
!  
!-----------------------------------------------------------------------------------------

   write(o_logUnit,'(2a)') subname,"Enter" ;  call shr_sys_flush(o_logUnit)
   if (debug>0) write(o_logUnit,'(2a,i2)') subname,"debug level = ",debug

   !--------------------------------------------------------------------------------------
   ! debug/sanity-check on data
   !--------------------------------------------------------------------------------------
   if (debug > 0) then
      write(o_logUnit,'(2a,i5)') subName,"pop lsize = ",mct_aVect_lsize(p2x_2d_p)
      do k=1,4 ! four curtains: N,E,S,W
         if (k==k_Scurtain) then
            lsize = mct_aVect_lsize(p2x_2d_rc(k_Scurtain))
            write(o_logUnit,'(2a,i2,L2,i4)') subname,"S curtain: k, do_mapping, lsize = ",k,do_Scurtain,lsize
         end if
         if (k==k_Ecurtain) then
            lsize = mct_aVect_lsize(p2x_2d_rc(k_Ecurtain))
            write(o_logUnit,'(2a,i2,L2,i4)') subname,"E curtain: k, do_mapping, lsize = ",k,do_Ecurtain,lsize
         end if
         if (k==k_Ncurtain) then
            lsize = mct_aVect_lsize(p2x_2d_rc(k_Ncurtain))
            write(o_logUnit,'(2a,i2,L2,i4)') subname,"N curtain: k, do_mapping, lsize = ",k,do_Ncurtain,lsize
         end if
         if (k==k_Wcurtain) then
            lsize = mct_aVect_lsize(p2x_2d_rc(k_Wcurtain))
            write(o_logUnit,'(2a,i2,L2,i4)') subname,"W curtain: k, do_mapping, lsize = ",k,do_Wcurtain,lsize
         end if
      end do
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
         lsize = mct_aVect_lsize(p2x_2d_rc(k))
         if ( lsize > 0 ) then
            kfld = mct_aVect_indexRA(p2x_2d_rc(k),"So_ssh" )
            xmin = minval(p2x_2d_rc(k)%rAttr(kfld,:))
            xmax = maxval(p2x_2d_rc(k)%rAttr(kfld,:))
            write(o_logUnit,'(2a,i2,2es11.3)') subname,"2d map: k, rc ssh min,max = ",k,xmin,xmax 
         else
            write(o_logUnit,'(2a,i3,a)')       subName,"2d map: k = ",k," lsize = 0 => no mapping"
         end if
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
      depth_p(:) = 0.01_R8*pop_depth(:) ! convert cm -> m
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
   first_call = .false.

   do k=1,4 ! four curtains: N,E,S,W

      do_mapping = .false.
      if (k==k_Scurtain .and. do_Scurtain==.true.) do_mapping = .true.
      if (k==k_Ecurtain .and. do_Ecurtain==.true.) do_mapping = .true.
      if (k==k_Ncurtain .and. do_Ncurtain==.true.) do_mapping = .true.
      if (k==k_Wcurtain .and. do_Wcurtain==.true.) do_mapping = .true.

      if (do_mapping) then

         !--- horizontal mapping, for each pop vertical level (global operation) --------
         do np=1,nlev_p
            call mct_sMat_avMult(p2x_3d_p(np),sMatp_p2rc(k),p2x_3d_pvert_rc(k,np))
         end do

         lsize = mct_aVect_lsize (p2x_3d_rc(k,1))
         if (debug>0 .and. lsize > 0) then
            kfld = mct_aVect_indexRA(p2x_3d_rc(k,1),"So_temp" )
            xmin = minval(p2x_3d_pvert_rc(k,1)%rAttr(kfld,:))
            xmax = maxval(p2x_3d_pvert_rc(k,1)%rAttr(kfld,:))
            write(o_logUnit,'(2a,i2,2es11.3)') subname,"3d map: k, pvert_rc SST min,max = ",k,xmin,xmax 
         end if

         !--- vertical interpolation to roms layers (local operation) -------------------
         if (lsize > 0) then
         do nr=1,nlev_r ! for each roms level
            call mct_aVect_zero(p2x_3d_rc(k,nr))

            do ij=1,lsize ! for each roms grid cell

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
                if (np1.eq.-1) write(o_logUnit,*) 'ERROR: failed to find bounding layers'

                !--- apply interpolation weights ---
                p2x_3d_rc(k,nr)%rAttr(:,ij) =  &
                       w1 * p2x_3d_pvert_rc(k,np1)%rAttr(:,ij) + &
                       w2 * p2x_3d_pvert_rc(k,np2)%rAttr(:,ij) 

                !--- document a sampling of interp specifics ---
                if (debug>1 .and. ((nr==nlev_r .or. mod(nr-1,10)==0) &
                            .and.  (ij==lsize  .or. mod(ij-1,10)==0))  ) then
                   zp1 = depth_p(np1) 
                   zp2 = depth_p(np2) 
                   zr  = depth_r
                   kfld = k_p2x_3d_So_temp
                   F1 = p2x_3d_pvert_rc(k,np1)%rAttr(kFld,ij) 
                   F2 = p2x_3d_pvert_rc(k,np2)%rAttr(kFld,ij)
                   Fr = p2x_3d_rc      (k,nr )%rAttr(kFld,ij)
                   write(o_logUnit,'(2a,i2,i3,2i4,2i3,2f5.2,3f7.1,3es11.2)') subName,&
                      "3d map: k,nr,i,j,np1,np2,w1,w2,z1,z,z2,F1,Fr,F2=", &
                               k,nr,i,j,np1,np2,w1,w2,zp1,zr,zp2,F1,Fr,F2
                end if

            end do ! do ij - loop over roms grid cells
         end do ! do nr ~ loop over roms layers
         end if ! lsize > 0 ~ non-zero tile size
         write(o_logUnit,'(2a,i3,a)') subName,"3d map: k = ",k," lsize = 0 => no mapping"
      end if ! do_mapping => curtain is active

   end do  ! do k=1,4  ~ over all four curtains: N,E,S,W

   write(o_logUnit,'(2a)') subname,"Exit" ;  call shr_sys_flush(o_logUnit)

end subroutine ocpl_map_pop2roms

!=========================================================================================
!BOP !====================================================================================
!
! !IROUTINE: ocpl_map_roms2pop
!
! !DESCRIPTION:
!     map data from roms to pop, for 3d restoring of pop to roms state 
!
! !INTERFACE:
!
! !REVISION HISTORY:
!    2020 Jan -- Brian Kauffman, initial version
!
! !INTERFACE: ----------------------------------------------------------------------------

subroutine ocpl_map_roms2pop()

! !INPUT/OUTPUT PARAMETERS:

   use grid,         only: z_p         => zw     ! pop depth array
   use mod_grid,     only: ROMS_GRID   => GRID
   use mod_param,    only: ROMS_BOUNDS => BOUNDS, & 
                           globalISize0 => Lm, &
                           globalJSize0 => Mm  
   use mod_parallel, only: MyRank

!EOP
!BOC

   !----- local variables -----
   integer(IN)        :: i,j,n                   ! spatial index wrt data(i,j) or data(n)
   integer(IN)        :: k,k_p,k_r,k0_r,k1_r     ! vertical index
   integer(IN)        :: iMin,iMax,jMin,jMax     ! range for roms(i,j), for local tile
   integer(IN)        :: LBi,LBj,UBi,UBj         ! range for roms(i,j), for global grid
   integer(IN)        :: localIsize,localJsize   ! roms array size
   integer(IN)        :: k_temp                  ! temp index for a field in an rAttr
   real(r8)           :: zMin_r,zMax_r           ! min/max roms column depth
   real(r8)           :: z0_r,z1_r               ! roms UB/LB column depths
   real(r8)           :: f0,f1                   ! vertical interp weights
   integer(IN)        :: lSize                   ! local aVect size
!  real(r8)           :: pmin, pmax, rmin, rmax  ! min/max values for debuging
   real(r8)           :: xmin, xmax              ! min/max values for debugging
   real(r8)           :: diff                    ! distance of index from edge of domain
   real(r8),parameter :: eps = 1.0e-12           ! epsilon for error checking
   logical            :: first_call = .true.

   character(*), parameter :: subName = "(ocpl_map_roms2pop) "

!-----------------------------------------------------------------------------------------
!  
!-----------------------------------------------------------------------------------------
 
   if (debug>0) write(o_logUnit,'(2a)') subname,"Enter" ;  call shr_sys_flush(o_logUnit)

   !--- for accessing roms arrays with (i,j) indexing ---
   iMin = ROMS_BOUNDS(nestID)%IstrR(MyRank)
   iMax = ROMS_BOUNDS(nestID)%IendR(MyRank)
   jMin = ROMS_BOUNDS(nestID)%JstrR(MyRank)
   jMax = ROMS_BOUNDS(nestID)%JendR(MyRank)
   localISize = iMax - iMin + 1
   localJSize = jMax - jMin + 1
   LBi  = 1
   LBj  = 1                              
   UBi  = globalISize0(nestID) + 1
   UBj  = globalJSize0(nestID) + 1

   !--- create weights that are 1.0 in roms interior, ramping down to 0.0 at roms boundary ---
   if (first_call) then
      write(o_logUnit,'(2a,4i7)') subname,"roms grid:  local: LBi,LBj,UBi,UBj=",iMin,jMin,iMax,jMax
      write(o_logUnit,'(2a,4i7)') subname,"roms grid, global: LBi,LBj,UBi,UBj=",LBi,LBj,UBi,UBj
      if ( mct_aVect_lSize(r2x_2d_r) == 0) then  ! local tile size = 0 (not data on this PE)
         write(o_logUnit,'(2a)') subname,"create roms domain edge weights (this tile lsize = 0)" 
         write(o_logUnit,'(2a)') subname,"roms edge weight min/max = N/A (lsize=0)"
      else
         write(o_logUnit,'(2a)') subname,"create roms domain edge weights" ;  call shr_sys_flush(o_logUnit)
         do n = 1, mct_aVect_lSize(r2x_2d_r) ! loop over roms aVect(n), need corresponding roms global i,j 
            i = mod(n-1,localISize) + iMin
            j = n/localISize        + jMin     ! requires/assumes fortran truncation

            !---- kludge: ramp to zero over 40 cells adjacent to boundary, need !a better algorithm
            diff = 1.0_r8
            if (abs(i-LBi)<40) diff = min(diff,abs(i-LBi)/40.0_r8)
            if (abs(i-UBi)<40) diff = min(diff,abs(i-UBi)/40.0_r8)
            if (abs(j-LBj)<40) diff = min(diff,abs(j-LBj)/40.0_r8)
            if (abs(j-UBj)<40) diff = min(diff,abs(j-UBj)/40.0_r8)
            r2x_2d_r%rAttr(k_r2x_2d_wgts,n) = diff
         end do
         if (debug>0 ) then
            xmin = minval(r2x_2d_r%rAttr(k_r2x_2d_wgts,:))
            xmax = maxval(r2x_2d_r%rAttr(k_r2x_2d_wgts,:))
            write(o_logUnit,'(2a,2f6.3)') subname,"roms edge weight min/max = ",xmin,xmax 
            call shr_sys_flush(o_logUnit)
         end if
      end if

      !----- map merge weights: roms -> pop -----
      call mct_sMat_avMult(r2x_2d_r, sMatp_r2p, r2x_2d_p,vector=usevector)

      if (debug>0 ) then
         if ( mct_aVect_lSize(r2x_2d_p) ==  0) then  ! local tile size = 0 (not data on this PE)
            write(o_logUnit,'(2a      )') subname,"pop  edge weight min/max = N/A (lsize=0)"
         else
            xmin = minval(r2x_2d_p%rAttr(k_r2x_2d_wgts,:))
            xmax = maxval(r2x_2d_p%rAttr(k_r2x_2d_wgts,:))
            write(o_logUnit,'(2a,2f6.3)') subname,"pop  edge weight min/max = ",xmin,xmax 
         end if
         call shr_sys_flush(o_logUnit)
      end if  

   end if
   first_call = .false.

   !--------------------------------------------------------------------------------------
   ! Step 1) vertical interpolation: roms levels -> pop levels
   ! Notes:
   ! o need roms(i,j) indicies to access internal roms cell depth array 
   ! o  "upper bound" means a roms level *deeper* than the target pop level
   ! o for pop,  level k+1 is deeper    than level k
   ! o for roms, level k+1 is shallower than level k
   ! o for pop,  z_p(k) units of depth are cm and are all positive numbers
   ! o for roms, z_r(k) units of depth are  m and are all negative numbers
   !--------------------------------------------------------------------------------------

!  alternate looping option
!  do j = jMin,jMax  ! j wrt roms(i,j)   ! loop over roms(i,j) and compute aVect(n)
!  do i = iMin,iMax  ! i wrt roms(i,j)
!  n = (j-jMin)*localISize + i - iMin + 1

   r2x_2d_r%rAttr(k_r2x_2d_reslev,:) = 1.0_r8 ! always restore at least one level (surface)

   if (debug>2) then
      write(o_logUnit,'(2a,5i5)') subName,"DEBUG: iMin,jMin,nLev_p,nLev_r,nLev_rp=",iMin,jMin,nLev_p,nLev_r,nLev_rp
   end if

   do k_p = 1,nLev_rp                    ! loop over pop levels (those involved in pop restoring)
   do n = 1,mct_aVect_lSize(r2x_3d_r(1)) ! loop over aVect(n), need corresponding roms global i,j 
      i = mod(n-1,localISize) + iMin
      j = n/localISize        + jMin     ! requires/assumes fortran truncation

      zMin_r = -ROMS_GRID(nestID)%z_r(i,j,nLev_r)*100.0_r8  ! shallowest roms level (convert m to cm)
      zMax_r = -ROMS_GRID(nestID)%z_r(i,j,     1)*100.0_r8  ! deepest    roms level (convert m to cm)

      if (debug>1 .and. mod(n-1,400)==0 ) then
      !  write(o_logUnit,'(2a,6i5)'  ) subName,"DEBUG: iMin,jMin,i,j,n,n2=",iMin,jMin,i,j,n,(j-jMin)*localISize + i - iMin + 1
         write(o_logUnit,'(2a,5i5)'  ) subName,"DEBUG: i,j,n,k_p   =",i,j,n,k_p
         write(o_logUnit,'(2a,3f9.1)') subName,"DEBUG: zMin_r,z_p(k_p),zMax_r=",zMin_r,z_p(k_p),zMax_r
      end if

      !----- if( pop is shallower than roms) then (restore down to this pop level) -----
      if (z_p(k_p) <= zMax_r) r2x_2d_r%rAttr(k_r2x_2d_reslev,n) = float(k_p) 

      if (z_p(k_p) >= zMax_r) then
         !----- target pop level deeper than max roms level => use max roms level -----
         r2x_3d_rp(k_p)%rAttr(:,n) = r2x_3d_r(nLev_r)%rAttr(:,n)
         if (debug>1 .and. mod(n-1,400)==0 ) then
            write(o_logUnit,'(2a)') subName,"  DEBUG: use max roms level"
         end if

      else if (z_p(k_p) <= zMin_r) then
         !----- target pop level shallower than min roms level => use min roms level -----
         r2x_3d_rp(k_p)%rAttr(:,n) = r2x_3d_r(1     )%rAttr(:,n)
         if (debug>1 .and. mod(n-1,400)==0 ) then
            write(o_logUnit,'(2a)') subName,"  DEBUG: use min roms level"
         end if
      else
         if (debug>1 .and. mod(n-1,400)==0 ) then
            write(o_logUnit,'(2a)') subName,"  DEBUG: vertical interpolation..."
         end if

         !----- search roms levels for UpperBound & LowerBound: z0_r < z_p <= z1_r -----
         z0_r = -ROMS_GRID(nestID)%z_r(i,j,nLev_r  )*100.0_r8  ! z0_r is the min roms level        (convert m to cm)
         do k_r=nLev_r,2,-1                                    ! start at surface, then go deeper
            z1_r = -ROMS_GRID(nestID)%z_r(i,j,k_r-1)*100.0_r8  ! z1_r is a candidate UB roms level (convert m to cm)
            if (z_p(k_p) <= z1_r) exit                         ! z1_r is the UB we're looking for, exit
            z0_r = z1_r                                        ! z1_r is NOT an UB, look deeper
         end do

         !----- compute UB & LB interp weights ----- z0_r < z_p <= z1_r 
         k0_r = k_r                                       ! LB roms index
         k1_r = k_r - 1                                   ! UB roms index
         f1   = (z_p(k_p) - z0_r)/(z1_r - z0_r)           ! fraction for UB field value
         f0   = 1.0_r8 - f1                               ! fraction for LB field value

         !----- do the linear vertical interpolation -----
         r2x_3d_rp(k_p)%rAttr(:,n) = f0* r2x_3d_r(k0_r)%rAttr(:,n) + f1* r2x_3d_r(k1_r)%rAttr(:,n)

         !----- check for errors ----
         if (f0 < -eps .or. f0 > (1.0_r8+eps) ) then
            write(o_logUnit,*) subName," ERROR: i,j,n,k_p,z_p(k_p) =", i,j,n, k_p,z_p(k_p)
            write(o_logUnit,*) subName," ERROR: k0_r,k1_r,z0_r,z1_r=", k0_r,k1_r, z0_r,z1_r
            write(o_logUnit,*) subName," ERROR: f0,f1              =", f0,f1
            call shr_sys_abort(subName//"ERROR: vertical interp")
         end if
         if (debug>1 .and. mod(n-1,400)==0 ) then
            write(o_logUnit,'(2a,5i5  )') subName,"  DEBUG: i,j,n,k_p,k1_r =", i,j,n,k_p,k1_r
            write(o_logUnit,'(2a,2f8.3)') subName,"  DEBUG: f0,f1 =", f0,f1
            write(o_logUnit,'(2a,3f8.1)') subName,"  DEBUG: z_LB_r,z_target,z_UB =",z0_r,z_p(k_p),z1_r
            write(o_logUnit,'(2a,3f8.1)') subName,"  DEBUG: T_LB_r,T_target,T_UB ="  &
            &                                , r2x_3d_r (k0_r)%rAttr(k_temp,n) &
            &                                , r2x_3d_rp( k_p)%rAttr(k_temp,n) &
            &                                , r2x_3d_r (k1_r)%rAttr(k_temp,n)
         end if
     end if

   end do ! n   = 1,mct_aVect_lsize(r2x_3d_r)
   end do ! k_p = 1,nLev_rp 


   !--------------------------------------------------------------------------------------
   ! Step 2) horizonal interpolation: roms -> pop  (all data at pop vertical levals)
   !--------------------------------------------------------------------------------------
   if (debug>0) write(o_logUnit,'(2a)') subname,"map 3d horizontal" ;  call shr_sys_flush(o_logUnit)

   call mct_sMat_avMult(r2x_2d_r, sMatp_r2p, r2x_2d_p,vector=usevector)

   do k_p = 1,nLev_rp  
      if (debug>0) write(o_logUnit,'(2a,i3)') subname,"map 3d horizontal, level=",k_p;  call shr_sys_flush(o_logUnit)
      call mct_sMat_avMult(r2x_3d_rp(k_p), sMatp_r2p, r2x_3d_p(k_p),vector=usevector)
   end do 

end subroutine ocpl_map_roms2pop

!=========================================================================================
!=========================================================================================
end module ocpl_map_mod
