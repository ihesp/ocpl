module ocpl_data_mod

!=========================================================================================
! MODULE: ocpl_data
!
! DESCRIPTION:
!    This module holds global data structures for the ocpl ocean driver 
!
! !REVISION HISTORY:
!    2018-Oct: first version by Brian Kauffman 
!=========================================================================================

   use mct_mod, only: mct_aVect, mct_gsmap, mct_ggrid
   use shr_kind_mod, only: IN=>SHR_KIND_IN,                  CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
   use seq_infodata_mod, only: seq_infodata_type
   use seq_cdata_mod, only : seq_cdata

   implicit none
   public                 ! default is data public, that's the point of this module

!  PUBLIC FUNCTIONS:

   ! none

!  PUBLIC DATA:

   character(CL) :: start_type
   character(CS) :: case_name
   integer(IN)   :: o_logUnit = 6 ! ocn   unit normally ocn.log
   integer(IN)   :: d_logUnit = 6 ! debug unit, ocn.log for PID==0, else stdout or /dev/nul 

   !--------------------------------------------------------------------------------------
   ! "_o" data -- ocean domain coupler interacts with, TENTITIVELY is the pop domain
   ! "_p" data -- pop/global    domain, coupler doesn't know about this
   ! "_r" data -- roms/regional domain, coupler doesn't know about this
   !--------------------------------------------------------------------------------------

   !--- pop  output: orig & re-gridded ---
   type(mct_aVect)         :: x2o_p          ! pop  input: x2o_o mapped to p grid 
   type(mct_aVect)         :: p2x_p          ! pop ouptut 
   type(mct_aVect)         :: p2x_o          ! pop ouptut: mapped to ocpl grid
   type(mct_aVect)         :: p2x_2d_p       ! pop output: 2D fields for pop/roms coupling
   type(mct_aVect),pointer :: p2x_3d_p(:)    ! pop output: 3D fields for pop/roms coupling (level)
   type(mct_aVect),pointer :: p2x_2d_rc(:)   ! pop output: 2D fields on roms curtain (curtain)
   type(mct_aVect),pointer :: p2x_3d_rc(:,:) ! pop output: 3D fields on roms curtain (curtain,level)

   integer(IN),parameter   :: k_Scurtain = 4 ! roms BC curtain index, South
   integer(IN),parameter   :: k_Ecurtain = 3 ! roms BC curtain index, East
   integer(IN),parameter   :: k_Ncurtain = 2 ! roms BC curtain index, North
   integer(IN),parameter   :: k_Wcurtain = 1 ! roms BC curtain index, West

   logical,parameter :: do_Scurtain = .true.  ! true => roms needs BCs for south curtain
   logical,parameter :: do_Ecurtain = .true.  ! true => roms needs BCs for  east curtain
   logical,parameter :: do_Ncurtain = .true.  ! true => roms needs BCs for north curtain
   logical,parameter :: do_Wcurtain = .true.  ! true => roms needs BCs for  west curtain

   !--- roms output: orig & re-gridded ---
   type(mct_aVect)         :: x2o_r          ! roms  input: x2o_o mapped to r grid 
   type(mct_aVect)         :: r2x_r          ! roms ouptut 
   type(mct_aVect)         :: r2x_o          ! roms ouptut: r2x_r mapped to o grid
   type(mct_aVect),pointer :: r2x_3d_r (:)   ! roms horiz grid, roms vert grid, (level)
   type(mct_aVect),pointer :: r2x_3d_rp(:)   ! roms horiz grid,  pop vert grid, (level)
   type(mct_aVect),pointer :: r2x_3d_p (:)   ! pop  horiz grid,  pop vert grid, (level)
   type(mct_aVect)         :: r2x_2d_r       ! additional 2d data needed wrt 3d data
   type(mct_aVect)         :: r2x_2d_p       ! additional 2d data needed wrt 3d data

   !--- roms & pop: number of vertical levels ---
   integer(IN)             :: nlev_p         ! pop : number of vertical levels
   integer(IN),parameter   :: nlev_rp = 34   ! number of pop restoring levels (34=> 527m)
   integer(IN)             :: nlev_r         ! roms: number of levels 

   !--- fundamental cesm coupler & mct data structures ---
   type(seq_cdata), pointer :: cdata_o       ! coupler instantiates this
   type(seq_cdata)          :: cdata_p       ! pop/regional
   type(seq_cdata)          :: cdata_r       ! roms/regional

   type(mct_gGrid), pointer :: gGrid_o      ! gGrid for ocpl (alloc'd in cpl)
   type(mct_gGrid), target  :: gGrid_r      ! gGrid for roms (alloc'd here)
   type(mct_gGrid), target  :: gGrid_p      ! gGrid for pop  (alloc'd here)

   type(mct_gsMap), pointer :: gsMap_o      ! gsMap for ocpl         (alloc'd in cpl)
   type(mct_gsMap), target  :: gsMap_p      ! gsMap for pop surface  (alloc'd here)
   type(mct_gsMap), target  :: gsMap_r      ! gsMap for roms surface (alloc'd here)
   type(mct_gsMap), pointer :: gsMap_rc(:)  ! gsMaps for 4 roms curtians (must alloc somewhere)

   integer(IN)              :: lsize_o      ! ocpl local 1d data size
   integer(IN)              :: lsize_p      ! pop  local 1d data size
   integer(IN)              :: lsize_r      ! roms local 1d data size

   integer(IN)              :: ni_o,nj_o    ! ocpl global 2d data size
   integer(IN)              :: ni_p,nj_p    ! pop  global 2d data size
   integer(IN)              :: ni_r,nj_r    ! roms global 2d data size

   type(seq_infodata_type), pointer :: infodata_o
   integer(IN)                      :: mpicom_o, OCNID_o ! all ocn comps share same mpi comm, ID
   character(16)                    :: name_o

!  type(seq_infodata_type), target  :: infodata_p   pop shares infodata_o with pop
   integer(IN)                      :: mpicom_p, OCNID_p
   character(16)                    :: name_p

!  type(seq_infodata_type), target  :: infodata_r   roms shares infodata_o with pop
   integer(IN)                      :: mpicom_r, OCNID_r
   character(16)                    :: name_r

#ifdef CPP_VECTOR
   logical :: usevector = .true.
#else
   logical :: usevector = .false.
#endif

   save ! save everything.

!================================================================================
end module ocpl_data_mod
