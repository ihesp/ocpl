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

   use mct_mod
!  use shr_kind_mod, only: IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
   use shr_kind_mod, only: IN=>SHR_KIND_IN,                  CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
   use shr_mct_mod
   use seq_infodata_mod
   use seq_cdata_mod

   implicit none
   public                 ! default is data public, that's the point of this module

!  PUBLIC FUNCTIONS:

   ! none

!  PUBLIC DATA:

   character(CL)   :: start_type
   character(CS)   :: case_name
   integer(IN)     :: o_logUnit = 6

   !--------------------------------------------------------------------------------------
   ! "_o" data -- ocean domain coupler interacts with, TENTITIVELY is the pop domain
   ! "_p" data -- pop/global    domain, coupler doesn't know about this
   ! "_r" data -- roms/regional domain, coupler doesn't know about this
   !--------------------------------------------------------------------------------------

   !--- pop  output: orig & re-gridded ---
!  type(mct_aVect)         :: x2o_p          ! pop  input: x2o_o mapped to p grid 
!  type(mct_aVect)         :: p2x_p          ! pop ouptut 
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
   logical,parameter :: do_Wcurtain = .false. ! true => roms needs BCs for  west curtain

   !--- roms output: orig & re-gridded ---
   type(mct_aVect)         :: x2o_r          ! roms  input: x2o_o mapped to r grid 
   type(mct_aVect)         :: r2x_r          ! roms ouptut 
   type(mct_aVect)         :: r2x_o          ! roms ouptut: r2x_r mapped to o grid
   type(mct_aVect),pointer :: r2x_3d_r (:)   ! roms horiz grid, roms vert grid, (level)
   type(mct_aVect),pointer :: r2x_3d_rp(:)   ! roms horiz grid,  pop vert grid, (level)
   type(mct_aVect),pointer :: r2x_3d_p (:)   ! pop  horiz grid,  pop vert grid, (level)

   !--- roms & pop: number of vertical levels ---
   integer(IN)             :: nlev_p         ! pop : number of vertical levels
   integer(IN)             :: nlev_pr        ! pop : number of vertical levels to be restored
   integer(IN)             :: nlev_r         ! roms: number of levels 

   !--- fundamental cesm coupler & mct data structures ---
!  type(seq_cdata), pointer :: cdata_o       ! coupler instantiates this
   type(seq_cdata)          :: cdata_r       ! roms/regional

   type(mct_gGrid), pointer :: dom_o 
   type(mct_gGrid), target  :: dom_r

   type(mct_gsMap), pointer :: gsMap_o      ! gsMap for pop surface
   type(mct_gsMap), target  :: gsMap_r      ! gsMap for roms surface
   type(mct_gsMap), pointer :: gsMap_rc(:)  ! gsMaps for four roms curtians

   type(seq_infodata_type), pointer :: infodata_o
   integer(IN)                      :: mpicom_o, OCNID_o ! all ocn comps share same mpi comm, ID
   character(16)                    :: name_o

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
