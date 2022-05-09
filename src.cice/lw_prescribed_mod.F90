!===================================================================
!BOP 
!
! !MODULE: lw - Prescribed Longwave Model
!
! !DESCRIPTION:
!
! 2018-Nov-28 Lantao Sun
! The modified prescribed lw model reads in target sea ice volume from a netCDF
! file. 
! The module then calculates the downward lw flux based on the formula in Sun,
! Alexander and Deser (2018) Journal of Climate, by using 
!       lwflx = rLfi*(vice - vice_target)/nudge_timescale
! where rLfi is the latent heat of fusion for sea ice times ice density, vice is the model
! instantaneous sea ice volume at each grid point, vice_target is the target sea
! ice volume, nudge_timescale is the relaxation time (default to be 5 days).

! The prescribed lw model reads in downward lw flux data from a netCDF
! file.  The prescribed downward lw flux is added to the flux
! sent to the ice model from the atmosphere, when it is run in a fully coupled mode.
! Regridding and data cycling capabilities are included.
! based on the prescibed ice model 
!
! 2020-May-15 - Alexandre Audette - modified to split bottom heat flux between bottom and lateral melt
! 2018-Nov-28 - Lantao Sun - modified to do instantaenous siv nudging 
! 2012-Aug-06 - Bob Tomas - test with the release 07 code  
!
! !REVISION HISTORY:
!  SVN:$Id: ice_prescribed_mod.F90 40 2006-12-01 19:09:30Z eclare $
!
! 2010-May-15 - Tony Craig and Mariana Vertenstein - updated to latest streams
! 2006-Aug-22 - D. Bailey, E. Hunke, modified to fit with CICE
! 2005-May-19 - J. Schramm - first version
! 2005-Apr-19 - B. Kauffman, J. Schramm, M. Vertenstein, NCAR - design
!
! !INTERFACE: ----------------------------------------------------------
 
module lw_prescribed_mod

! !USES:

   use shr_strdata_mod
   use shr_dmodel_mod
   use shr_string_mod
   use shr_ncread_mod
   use shr_sys_mod
   use shr_mct_mod
   use mct_mod
   use pio

   use ice_broadcast
   use ice_communicate, only : my_task, master_task, MPI_COMM_ICE
   use ice_kinds_mod
   use ice_fileunits
   use ice_exit,        only : abort_ice
   use ice_domain_size, only : nx_global, ny_global, ncat, nilyr, nslyr, max_blocks
   use ice_constants
   use ice_blocks,     only : nx_block, ny_block
   use ice_domain,     only : nblocks, distrb_info, blocks_ice
   use ice_grid,       only : TLAT,TLON,hm,tmask
   use ice_calendar,   only : idate, sec, calendar_type
   use ice_itd,        only : ilyr1, slyr1, hin_max
   use ice_read_write
! BT 06 August 2012 - this downward longwave flux will be updated to include a contribution from the values read from the file
!   use ice_flux,        only : flw

! Lantao Sun 2018-Nov-28 use ice volume to calculate the longwave flux
! Alexandre Audette 2020-May-15 use ice area to calculate longwave flux
   use ice_state, only : vice, aice

   implicit none
   save

   private ! except


! !PUBLIC TYPES:

! !PUBLIC MEMBER FUNCTIONS:

   public :: lw_prescribed_init      ! initialize input data stream
   public :: lw_prescribed_run       ! get time slices and time interp
!   public :: ice_prescribed_phys      ! set prescribed ice state and fluxes
! !PUBLIC DATA MEMBERS:

   logical(kind=log_kind), public :: prescribed_lw      ! true if prescribed lw 

!EOP

   logical(kind=log_kind)         :: prescribed_lw_fill ! true if data fill required
   integer(kind=int_kind)         :: stream_year_first   ! first year in stream to use
   integer(kind=int_kind)         :: stream_year_last    ! last year in stream to use
   integer(kind=int_kind)         :: model_year_align    ! align stream_year_first 
                                                         ! with this model year
   integer(SHR_KIND_IN),parameter :: nFilesMaximum = 100 ! max number of files

   character(len=char_len_long)   :: stream_fldVarName
   character(len=char_len_long)   :: stream_fldFileName(nFilesMaximum)
   character(len=char_len_long)   :: stream_domTvarName
   character(len=char_len_long)   :: stream_domXvarName
   character(len=char_len_long)   :: stream_domYvarName
   character(len=char_len_long)   :: stream_domAreaName
   character(len=char_len_long)   :: stream_domMaskName
   character(len=char_len_long)   :: stream_domFileName
   ! Alexandre Audette May-21-2020
   character(len=char_len_long)   :: sithick_fldVarName
   !character(len=char_len_long)   :: stream_fldFileName2(nFilesMaximum)

! Lantao Sun 2018-Nov-28 add logical paramter for NH and SH nudging seperately
   logical(kind=log_kind)         :: nh_nudge ! true if ice restore in the arctic
   logical(kind=log_kind)         :: sh_nudge ! true if ice restore in the AA
   
   type(shr_strdata_type)       :: sdat         ! prescribed data stream
   character(len=char_len_long) :: fldList      ! list of fields in data stream
   ! Alexandre Audette May-21-2020
   type(shr_strdata_type)       :: sdat2         ! prescribed data stream
!   real(kind=dbl_kind)          :: lw_pres(nx_block,ny_block,max_blocks) ! presecribed longwave 

! Alexandre Audette 2020-May-15 - read ice area and ice thickness
   real(kind=dbl_kind)          :: aice_target(nx_block,ny_block,max_blocks)
   real(kind=dbl_kind)          :: hice_target(nx_block,ny_block,max_blocks)

! Lantao Sun 2018-Nov-28 read ice volume for nudging target
   real(kind=dbl_kind)          :: vice_target(nx_block,ny_block,max_blocks) !

! Alexandre Audette 2020-May-15 - read ice area and ice thickness
!   real(kind=dlb_kind)         :: aice_target(nx_block,ny_block,max_blocks)
!   real(kind=dlb_kind)          :: hice_target(nx_block,ny_block,max_blocks)

! delta_flx added to the bottom of the ice (fbot)
   real(kind=dbl_kind), public  :: fbot_delta(nx_block,ny_block,max_blocks) !

! Alexandre Audette 2020-May-15 - fraction of ice that needs to be gone to reach target added to the side of the ice (rside) 
   real(kind=dbl_kind), public  :: rside_delta(nx_block,ny_block,max_blocks) ! 
   real(kind=dbl_kind) :: hice(nx_block,ny_block,max_blocks) ! ice thickness from volume and area
! Lantao Sun 2018-Nov-28 add nudging time scale and latent heat of fusion
   real (kind=dbl_kind)         :: bot_nudge_timescale
   real (kind=dbl_kind)         :: lat_nudge_timescale
   real (kind=dbl_kind), parameter :: &
       rLfi = Lfresh*rhoi   ! latent heat of fusion ice               (J/m^3)

! BT 6 Aug 2012 Don't need this
!    real (kind=dbl_kind), parameter :: &
!       cp_sno = 0.0_dbl_kind & ! specific heat of snow                (J/kg/K)
!    ,  rLfi = Lfresh*rhoi & ! latent heat of fusion ice               (J/m^3)
!    ,  rLfs = Lfresh*rhos & ! latent heat of fusion snow              (J/m^3)
!    ,  rLvi = Lvap*rhoi   & ! latent heat of vapor*rhoice             (J/m^3)
!    ,  rLvs = Lvap*rhos   & ! latent heat of vapor*rhosno             (J/m^3)
!    ,  rcpi = cp_ice*rhoi & ! heat capacity of fresh ice              (J/m^3)
!    ,  rcps = cp_sno*rhos & ! heat capacity of snow                   (J/m^3)
!    ,  rcpidepressT = rcpi*depressT & ! param for finding T(z) from q (J/m^3)
!    ,  rLfidepressT = rLfi*depressT ! param for heat capacity   (J deg/m^3)
!         ! heat capacity of sea ice, rhoi*C=rcpi+rLfidepressT*salinity/T^2
!
!=======================================================================
contains
!===============================================================================
!BOP
!
! !IROUTINE: lw_prescribed_init -  prescribed lw initialization
!
! !INTERFACE: 
 subroutine lw_prescribed_init(compid, gsmap, dom)
   use shr_pio_mod, only : shr_pio_getiotype, shr_pio_getiosys
! !DESCRIPTION:
!    Prescribed lw initialization - needed to 
!    work with new shr_strdata module derived type 
!
! !REVISION HISTORY:
!    2009-Oct-12 - M. Vertenstein
!
! !INPUT/OUTPUT PARAMETERS:
!
   implicit none
   include 'mpif.h'
   integer(kind=int_kind), intent(in) :: compid
   type(mct_gsMap) :: gsmap
   type(mct_gGrid) :: dom

!EOP
   !----- Local ------
   integer(kind=int_kind) :: nml_error ! namelist i/o error flag
   integer(kind=int_kind) :: n, nFile, nFile2, ierr   
   character(len=8)       :: fillalgo
   character(*),parameter :: subName = "('lw_prescribed_init2')"
   character(*),parameter :: F00 = "('(lw_prescribed_init2) ',4a)"

   namelist /lw_prescribed_nml/  &
        prescribed_lw,      &
        prescribed_lw_fill, &
	     stream_year_first ,  &
        stream_year_last  ,  &
        model_year_align,    &
        stream_fldVarName ,  &
        stream_fldFileName,  &
        stream_domTvarName,  &
        stream_domXvarName,  &
        stream_domYvarName,  &
        stream_domAreaName,  &
        stream_domMaskName,  &
        stream_domFileName,  &
        nh_nudge,            & !Lantao Sun Nov-28-2018 add three parameters
        sh_nudge,            &
        bot_nudge_timescale,     &
        sithick_fldVarName,  & !Alexandre Audette May-21-2020 added two parameters
        lat_nudge_timescale

   ! default values for namelist
   prescribed_lw         = .false.           ! if true, prescribe lw
   stream_year_first      = 1                ! first year in  pice stream to use
   stream_year_last       = 1                ! last  year in  pice stream to use
   model_year_align       = 1                ! align stream_year_first with this model year
   stream_fldVarName      = 'aice_target'
   sithick_fldVarName     = 'hice_target'
   stream_fldFileName(:)  = ' '
   stream_domTvarName     = 'time'
   stream_domXvarName     = 'lon'
   stream_domYvarName     = 'lat'
   stream_domAreaName     = 'area'
   stream_domMaskName     = 'mask'
   stream_domFileName     =  ' '
   prescribed_lw_fill    = .false.          ! true if plw data fill required

   !Lantao Sun Nov-28-2018
   nh_nudge               = .false.
   sh_nudge               = .false.
   bot_nudge_timescale    = 5.0 ! default to be 5 days
   lat_nudge_timescale    = 5.0 ! default to be 5 days
  
   ! read from input file
   call get_fileunit(nu_nml)
   if (my_task == master_task) then
      open (nu_nml, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      do while (nml_error > 0)
         read(nu_nml, nml=lw_prescribed_nml,iostat=nml_error)
         if (nml_error > 0) read(nu_nml,*)  ! for Nagware compiler
      end do
      if (nml_error == 0) close(nu_nml)
   endif
   call release_fileunit(nu_nml)
   call broadcast_scalar(nml_error,master_task)
   if (nml_error /= 0) then
      call abort_ice ('ice: Namelist read error in lw_prescribed_mod')
   endif

   call broadcast_scalar(prescribed_lw,master_task)

   ! *** If not prescribed lw then return ***
   if (.not. prescribed_lw) RETURN

   call broadcast_scalar(prescribed_lw_fill,master_task)
   call broadcast_scalar(stream_year_first,master_task)
   call broadcast_scalar(stream_year_last,master_task)
   call broadcast_scalar(model_year_align,master_task)
   call broadcast_scalar(stream_fldVarName,master_task)
   call broadcast_scalar(stream_domTvarName,master_task)
   call broadcast_scalar(stream_domXvarName,master_task)
   call broadcast_scalar(stream_domYvarName,master_task)
   call broadcast_scalar(stream_domAreaName,master_task)
   call broadcast_scalar(stream_domMaskName,master_task)
   call broadcast_scalar(stream_domFileName,master_task)

   ! Lantao Sun Nov-28-2018
   call broadcast_scalar(nh_nudge,master_task)
   call broadcast_scalar(sh_nudge,master_task)
   call broadcast_scalar(bot_nudge_timescale,master_task)

   ! Alexandre Audette May-21-2020
   call broadcast_scalar(sithick_fldVarName,master_task)
   !call broadcast_scalar(stream_fldFileName2,master_task)


   call mpi_bcast(stream_fldFileName, len(stream_fldFileName(1))*NFilesMaximum, &
        MPI_CHARACTER, 0, MPI_COMM_ICE, ierr)

   nFile = 0
   do n=1,nFilesMaximum
      if (stream_fldFileName(n) /= ' ') nFile = nFile + 1
   end do

   


   if (my_task == master_task) then
      write(nu_diag,*) ' '
      write(nu_diag,*) 'This is the prescribed lw coverage option.'
      write(nu_diag,*) '  stream_year_first  = ',stream_year_first  
      write(nu_diag,*) '  stream_year_last   = ',stream_year_last   
      write(nu_diag,*) '  model_year_align   = ',model_year_align   
      write(nu_diag,*) '  stream_fldVarName  = ',trim(stream_fldVarName)
      write(nu_diag,*) '  sithick_fldVarName = ',trim(sithick_fldVarName) ! Alexandre Audette
      do n = 1,nFile
         write(nu_diag,*) '  stream_fldFileName = ',trim(stream_fldFileName(n)),n
      end do
      write(nu_diag,*) '  stream_domTvarName = ',trim(stream_domTvarName)
      write(nu_diag,*) '  stream_domXvarName = ',trim(stream_domXvarName)
      write(nu_diag,*) '  stream_domYvarName = ',trim(stream_domYvarName)
      write(nu_diag,*) '  stream_domFileName = ',trim(stream_domFileName)
      write(nu_diag,*) '  prescribed_lw_fill= ',prescribed_lw_fill

      !Lantao Sun Nov-28-2018
      write(nu_diag,*) '  nh_nudge = ',nh_nudge
      write(nu_diag,*) '  sh_nudge = ',sh_nudge
      write(nu_diag,*) '  bot_nudge_timescale = ',bot_nudge_timescale

      write(nu_diag,*) ' '
   endif

   ! Read shr_strdata_nml namelist
   if (prescribed_lw_fill) then
      fillalgo='nn'
   else
      fillalgo='none'
   endif

   call shr_strdata_create(sdat,name="prescribed_lw", &
        mpicom=MPI_COMM_ICE, compid=compid, &
        gsmap=gsmap, ggrid=dom,          &
        nxg=nx_global,nyg=ny_global,     &
        yearFirst=stream_year_first,     &
        yearLast=stream_year_last,       &
        yearAlign=model_year_align,      &
        offset=0,                        &
        domFilePath='',                  &
        domFileName=trim(stream_domFileName), &
        domTvarName=stream_domTvarName,  &
        domXvarName=stream_domXvarName,  &
        domYvarName=stream_domYvarName,  &
        domAreaName=stream_domAreaName,  &
        domMaskName=stream_domMaskName,  &
        filePath='',                     &
        filename=stream_fldFileName(1:nFile), &
        fldListFile=stream_fldVarName,   &
        fldListModel=stream_fldVarName,  &
        pio_subsystem=shr_pio_getiosys('ICE'),&
        pio_iotype=shr_pio_getiotype('ICE'),&
        fillalgo = trim(fillalgo),       &
        calendar = trim(calendar_type))

   ! Alexandre Audette May-21-2020
   ! Read in the sea ice thickness field.
   call shr_strdata_create(sdat2,name="prescribed_lw_sithick", &
        mpicom=MPI_COMM_ICE, compid=compid, &
        gsmap=gsmap, ggrid=dom,          &
        nxg=nx_global,nyg=ny_global,     &
        yearFirst=stream_year_first,     &
        yearLast=stream_year_last,       &
        yearAlign=model_year_align,      &
        offset=0,                        &
        domFilePath='',                  &
        domFileName=trim(stream_domFileName), &
        domTvarName=stream_domTvarName,  &
        domXvarName=stream_domXvarName,  &
        domYvarName=stream_domYvarName,  &
        domAreaName=stream_domAreaName,  &
        domMaskName=stream_domMaskName,  &
        filePath='',                     &
        filename=stream_fldFileName(1:nFile), &
        fldListFile=sithick_fldVarName,   &
        fldListModel=sithick_fldVarName,  &
        pio_subsystem=shr_pio_getiosys('ICE'),&
        pio_iotype=shr_pio_getiotype('ICE'),&
        fillalgo = trim(fillalgo),       &
        calendar = trim(calendar_type))


   if (my_task == master_task) then
      call shr_strdata_print(sdat, 'SPRECICE data')
      call shr_strdata_print(sdat2,'SPRESICE data')
   endif

   !-----------------------------------------------------------------
   ! For one ice category, set hin_max(1) to something big
   !-----------------------------------------------------------------
!   if (ncat == 1) then
!      hin_max(1) = 999._dbl_kind
!   end if
end subroutine lw_prescribed_init
  
!=======================================================================
!BOP ===================================================================
!
! !IROUTINE: lw_prescribed_run -- Update presecribed longwave 
!
! !DESCRIPTION:
!
!  Finds two time slices bounding current model time, remaps if necessary
!
! !REVISION HISTORY:
!     2005-May-19 - J. Schramm - first version
!     2009-Oct-15 - M. Vertenstein - update to new data model changes
!
! !INTERFACE: -----------------------------------------------------------

subroutine lw_prescribed_run(mDateIn, secIn)

! !USES:
   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(kind=int_kind), intent(in) :: mDateIn  ! Current model date (yyyymmdd)
   integer(kind=int_kind), intent(in) :: secIn    ! Elapsed seconds on model date

!EOP

   integer(kind=int_kind) :: i,j,n,iblk       ! loop indices and counter
   integer(kind=int_kind) :: ilo,ihi,jlo,jhi  ! beginning and end of physical domain
   type (block)           :: this_block
   real(kind=dbl_kind)    :: lw_max         ! maximun ice concentration
   logical, save          :: first_time = .true.
   character(*),parameter :: subName = "('lw_prescribed_run')"
   character(*),parameter :: F00 = "('(lw_prescribed_run) ',a,2g20.13)"
 
   !------------------------------------------------------------------------
   ! Interpolate to new prescribed lw
   !------------------------------------------------------------------------
! BT 31 Jul 2012 - I believe the last two arguments in the presecibed ice code are optional
! For now, I will try omitting them in the prescribed lw code and see what happens
! I think the last one is for timers, I don't know what the second from last one is for
!
! BT 1 Aug 2012 - I get a compiler error if I omit the last two arguments:
! line 327.9: 1513-061 (S) Actual argument attributes do not match those specified by an accessible explicit interface.
! so I'll put them back in for now and later determine if I need to create new ones for this imposed longwave
! code
!   call shr_strdata_advance(sdat,mDateIn,SecIn)
   call shr_strdata_advance(sdat,mDateIn,SecIn,MPI_COMM_ICE,'cice_pice')
   call shr_strdata_advance(sdat2,mDateIn,SecIn,MPI_COMM_ICE,'cice_pice')
   
!   lw_pres(:,:,:) = c0  ! This initializes ghost cells as well 

   n=0
   do iblk = 1, nblocks
      this_block = get_block(blocks_ice(iblk),iblk)         
      ilo = this_block%ilo
      ihi = this_block%ihi
      jlo = this_block%jlo
      jhi = this_block%jhi
      
      do j = jlo, jhi
      do i = ilo, ihi
         n = n+1
! BT 31 Jul 2012 - retained sdat%avs from prescribed ice model; may need to change avs to something else
! 3 Aug 2012 trying to update the flw here
!         lw_pres(i,j,iblk) = sdat%avs(1)%rAttr(1,n)
!         flw(i,j,iblk) = flw(i,j,iblk) + lw_pres(i,j,iblk)

! Lantao Sun Nov-28-2018 - Alexandre Audette 2020-May-15 - Modified to add lateral melt
!         vice_target(i,j,iblk) = sdat%avs(1)%rAttr(1,n)
         ! Alexandre Audette May-21-2020
         aice_target(i,j,iblk) = sdat%avs(1)%rAttr(1,n)
         hice_target(i,j,iblk) = sdat2%avs(1)%rAttr(1,n)
         fbot_delta(i,j,iblk) = 0.0_dbl_kind
         rside_delta(i,j,iblk)= 0.0_dbl_kind

               
         
         if( (nh_nudge .eq. .true. .and. TLAT(i,j,iblk)*rad_to_deg >  40.0_dbl_kind) .or. &
             (sh_nudge .eq. .true. .and. TLAT(i,j,iblk)*rad_to_deg < -40.0_dbl_kind) ) then
           !fbot_delta(i,j,iblk) = rLfi*(vice(i,j,iblk)-vice_target(i,j,iblk))/(nudge_timescale*secday)

          
           if (aice(i,j,iblk) .gt. 0.0_dbl_kind) then
            hice(i,j,iblk) = vice(i,j,iblk)/aice(i,j,iblk)
           else
            hice(i,j,iblk) = 0.0_dbl_kind
           end if
           
           rside_delta(i,j,iblk) = (aice(i,j,iblk) - aice_target(i,j,iblk))/(lat_nudge_timescale*secday)
           fbot_delta(i,j,iblk) = rLfi*aice(i,j,iblk)*(hice(i,j,iblk)-hice_target(i,j,iblk))/(bot_nudge_timescale*secday)
         else 
           fbot_delta(i,j,iblk) = 0.0_dbl_kind
           rside_delta(i,j,iblk) = 0.0_dbl_kind ! Alexandre Audette May-21-2020
         end if

      end do
      end do
   end do

! BT 31 Jul 2012 - don't need this check; maybe some other error checking but remove this one for now
   !--------------------------------------------------------------------
   ! Check to see that ice concentration is in fraction, not percent
   !--------------------------------------------------------------------
!   if (first_time) then
!      aice_max = maxval(ice_cov)
!
!      if (aice_max > c10) then
!         write(nu_diag,F00) "ERROR: Ice conc data must be in fraction, aice_max= ",&
!              aice_max
!         call abort_ice(subName)
!      end if
!      first_time = .false.
!   end if
!
   !-----------------------------------------------------------------
   ! Set prescribed ice state and fluxes
   !-----------------------------------------------------------------
!
!   call ice_prescribed_phys()
!
end subroutine lw_prescribed_run

!==============================================================================

end module lw_prescribed_mod

!==============================================================================
