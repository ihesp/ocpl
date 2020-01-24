module ocpl_fields_mod

!===============================================================================
!BOP ===========================================================================
!
! !MODULE: ocpl_fields -- ocpl/pop 3d exchange fields
!
! !DESCRIPTION:
!     
! !REVISION HISTORY:
!     2019 Sep - B. Kauffman
!
! !INTERFACE: ------------------------------------------------------------------


! !USES:

   use shr_kind_mod, only : IN => SHR_KIND_IN

   public  ! except

! !PUBLIC TYPES:

  ! none

! !PUBLIC MEMBER FUNCTIONS:

! !PUBLIC DATA MEMBER:

   !----------------------------------------------------------------------------
   ! exchange fields: roms-> pop
   !----------------------------------------------------------------------------
   ! 2d data in support of 3d data exchange
   character(*), parameter :: ocpl_fields_r2x_2d_fields = &
         'reslev'      &    ! max restoring level in pop      DEF
      //':wgts'             ! merge weights ROMS+POP fields   DEF

   ! States
   character(*), parameter :: ocpl_fields_r2x_3d_states = &
         'So_temp'     &    ! temperature from ROMS           DEF
      //':So_salt'          ! salinity    from ROMS           DEF

   ! Fluxes
   character(*), parameter :: ocpl_fields_r2x_3d_fluxes = &
         '      '           !                                 DEF

   character(*), parameter :: ocpl_fields_r2x_3d_fields = &
      trim(ocpl_fields_r2x_3d_states)//":"//trim(ocpl_fields_r2x_3d_fluxes)

   !----------------------------------------------------------------------------
   ! exchange fields: pop -> roms  2d & 3d boundary data
   !----------------------------------------------------------------------------
   ! States
   character(*), parameter :: ocpl_fields_p2x_2d_states = &
         'So_ssh'      &    ! sea surface height              DEF
      //':So_ubar'     &    ! vertically integrated velocity  DEF
      //':So_vbar'          ! vertically integrated velocity  DEF

   ! Fluxes
   character(*), parameter :: ocpl_fields_p2x_2d_fluxes = &
         '      '           !                                 DEF

   character(*), parameter :: ocpl_fields_p2x_2d_fields = &
      trim(ocpl_fields_p2x_2d_states)//":"//trim(ocpl_fields_p2x_2d_fluxes)

   ! States
   character(*), parameter :: ocpl_fields_p2x_3d_states = &
         'So_temp'     &    ! ocean potential temperature     DEF
      //':So_salt'     &    ! ocean salinity                  DEF
      //':So_uvel'     &    ! ocean velocity in grid-x dir    DEF
      //':So_vvel'          ! ocean velocity in grid-y dir    DEF

   ! Fluxes
   character(*), parameter :: ocpl_fields_p2x_3d_fluxes = &
         '      '           !                                 DEF

   character(*), parameter :: ocpl_fields_p2x_3d_fields = &
      trim(ocpl_fields_p2x_3d_states)//":"//trim(ocpl_fields_p2x_3d_fluxes)


   integer(IN) :: k_p2x_2d_So_ssh
   integer(IN) :: k_p2x_2d_So_ubar
   integer(IN) :: k_p2x_2d_So_vbar

   integer(IN) :: k_p2x_3d_So_temp
   integer(IN) :: k_p2x_3d_So_salt
   integer(IN) :: k_p2x_3d_So_uvel
   integer(IN) :: k_p2x_3d_So_vvel

   integer(IN) :: k_r2x_3d_So_temp
   integer(IN) :: k_r2x_3d_So_salt

   integer(IN) :: k_r2x_2d_reslev
   integer(IN) :: k_r2x_2d_wgts

   save ! save everything

end module ocpl_fields_mod
