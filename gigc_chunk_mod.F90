#include "MAPL_Generic.h"
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gigc_chunk_mod
!
! !DESCRIPTION: Module GC\_CHUNK\_MOD is the module that contains the init,
!  run, and finalize methods for the ESMF interface to GEOS-Chem.
!\\
!\\
! !INTERFACE: 
!      
MODULE GIGC_Chunk_Mod
!
! !USES:
!      
  USE MAPL_MOD
  USE ESMF
  USE ErrCode_Mod
  USE Precision_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: GIGC_Chunk_Init
  PUBLIC :: GIGC_Chunk_Run
  PUBLIC :: GIGC_Chunk_Final
!
! !PRIVATE MEMBER FUNCTIONS:
!
#if defined( MODEL_GEOS )
  PRIVATE :: SET_OZONOPAUSE
#endif

  INTEGER  ::  MemDebugLevel
!
! !REVISION HISTORY:
!  22 Jun 2009 - R. Yantosca & P. Le Sager - Chunkized & cleaned up.
!  09 Oct 2012 - R. Yantosca - Now pass am_I_Root to all routines
!  09 Oct 2012 - R. Yantosca - Added comments, cosmetic changes
!  16 Oct 2012 - R. Yantosca - Renamed GC_STATE argument to State_Chm
!  22 Oct 2012 - R. Yantosca - Renamed to gigc_chunk_mod.F90
!  01 Nov 2012 - R. Yantosca - Now pass Input Options object to routines
!  15 Mar 2013 - R. Yantosca - Add routine GIGC_Cap_Tropopause_Prs
!  08 Mar 2018 - E. Lundgren - Move gigc_initialization_mod contents to 
!                              gigc_chunk_init now that LOC is much reduced
!  14 Dec 2018 - E. Lundgren - Combine offline and GEOS-5 code and simplify
!EOP
!------------------------------------------------------------------------------
!BOC

CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gigc_chunk_init
!
! !DESCRIPTION: Subroutine GIGC\_CHUNK\_INIT is the ESMF init method for
!  GEOS-Chem.  This routine calls routines within core GEOS-Chem to allocate 
!  arrays and read input files.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GIGC_Chunk_Init( am_I_Root,                                  &
                              nymdB,         nhmsB,      nymdE,           &
                              nhmsE,         tsChem,     tsDyn,           &
                              lonCtr,        latCtr,     myPET,           &
#if !defined( MODEL_GEOS )
                              GC,            EXPORT,                      &
#endif
                              Input_Opt,     State_Chm,  State_Diag,      &
                              State_Grid,    State_Met,  HcoConfig,       &
                              HistoryConfig, RC )
!
! !USES:
!
    USE Chemistry_Mod,           ONLY : Init_Chemistry
    USE Emissions_Mod,           ONLY : Emissions_Init
    USE GC_Environment_Mod
    USE GC_Grid_Mod,             ONLY : SetGridFromCtr
    USE GIGC_HistoryExports_Mod, ONLY : HistoryConfigObj
    USE HCO_Types_Mod,           ONLY : ConfigObj
    USE Input_Mod,               ONLY : Read_Input_File
    USE Input_Opt_Mod,           ONLY : OptInput, Set_Input_Opt
    USE Linoz_Mod,               ONLY : Linoz_Read
    USE PBL_Mix_Mod,             ONLY : Init_PBL_Mix
    USE PhysConstants,           ONLY : PI_180
    USE Pressure_Mod,            ONLY : Init_Pressure
    USE Roundoff_Mod,            ONLY : RoundOff
    USE State_Chm_Mod,           ONLY : ChmState, Ind_
    USE State_Diag_Mod,          ONLY : DgnState
    USE State_Grid_Mod,          ONLY : GrdState, Init_State_Grid
    USE State_Met_Mod,           ONLY : MetState
    USE Strat_Chem_Mod,          ONLY : Init_Strat_Chem
#if defined( MODEL_GEOS )
    USE Tendencies_Mod,          ONLY : TEND_INIT
#endif
    USE Time_Mod,                ONLY : Set_Timesteps
    USE UCX_MOD,                 ONLY : INIT_UCX
    USE UnitConv_Mod,            ONLY : Convert_Spc_Units
#ifdef ADJOINT
    USE Charpak_Mod,             ONLY : To_UpperCase
#endif
!
! !INPUT PARAMETERS:
!
    LOGICAL,            INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
    INTEGER,            INTENT(IN)    :: myPET       ! Local PET
    INTEGER,            INTENT(IN)    :: nymdB       ! YYYYMMDD @ start of run
    INTEGER,            INTENT(IN)    :: nhmsB       ! hhmmss   @ start of run
    INTEGER,            INTENT(IN)    :: nymdE       ! YYYYMMDD @ end of run
    INTEGER,            INTENT(IN)    :: nhmsE       ! hhmmss   @ end of run
    REAL,               INTENT(IN)    :: tsChem      ! Chemistry timestep [s]
    REAL,               INTENT(IN)    :: tsDyn       ! Chemistry timestep [s]
    REAL(ESMF_KIND_R4), INTENT(IN)    :: lonCtr(:,:) ! Lon centers [radians]
    REAL(ESMF_KIND_R4), INTENT(IN)    :: latCtr(:,:) ! Lat centers [radians]
!
! !INPUT/OUTPUT PARAMETERS:
!
#if !defined( MODEL_GEOS )
    TYPE(ESMF_State),    INTENT(INOUT), TARGET :: EXPORT ! Export state object
    TYPE(ESMF_GridComp), INTENT(INOUT)         :: GC     ! Ref to this GridComp
#endif
    TYPE(OptInput),      INTENT(INOUT) :: Input_Opt      ! Input Options object
    TYPE(ChmState),      INTENT(INOUT) :: State_Chm      ! Chem State object 
    TYPE(DgnState),      INTENT(INOUT) :: State_Diag     ! Diag State object
    TYPE(GrdState),      INTENT(INOUT) :: State_Grid     ! Grid State object
    TYPE(MetState),      INTENT(INOUT) :: State_Met      ! Met State object
    TYPE(ConfigObj),     POINTER       :: HcoConfig      ! HEMCO config obj 
    TYPE(HistoryConfigObj), POINTER    :: HistoryConfig  ! History config obj 
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT)   :: RC             ! Success or failure?
!
! !REMARKS:
!  Need to add better error checking
!
! !REVISION HISTORY: 
!  18 Jul 2011 - M. Long     - Initial Version
!  28 Mar 2012 - M. Long     - Rewrite per structure of BCC init interface
!  09 Oct 2012 - R. Yantosca - Added comments, cosmetic changes
!  16 Oct 2012 - R. Yantosca - Renamed GC_STATE argument to State_Chm
!  16 Oct 2012 - R. Yantosca - Renamed GC_MET argument to State_Met
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_chm_mod.F90
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_met_mod.F90
!  22 Oct 2012 - R. Yantosca - Renamed to GIGC_Chunk_Init
!  01 Nov 2012 - R. Yantosca - Now reference gigc_input_opt_mod.F90
!  01 Nov 2012 - R. Yantosca - Reordered arguments for clarit
!  28 Nov 2012 - M. Long     - Now pass lonCtr, latCtr, latEdg as arguments
!                              to routine GIGC_Init_Simulation
!  03 Dec 2012 - R. Yantosca - Now call Init_CMN_SIZE (in CMN_SIZE_mod.F)
!                              instead of GIGC_Init_Dimensions to initialize
!                              the size parameters.
!  03 Dec 2012 - R. Yantosca - Rename NI, NJ, NL to IM, JM, LM for clarity
!  03 Dec 2012 - R. Yantosca - Now pass I_LO, J_LO, I_HI, J_HI, IM_WORLD, 
!                              JM_WORLD, LM_WORLD via the arg list
!  05 Dec 2012 - R. Yantosca - Remove latEdg argument
!  06 Dec 2012 - R. Yantosca - Add nymdB, nhmsB, nymdB, nhmsB arguments,
!                              and remove nymd, nhms
!  06 Mar 2018 - E. Lundgren - Remove Set_Initial_MixRatios
!  08 Mar 2018 - E. Lundgren - Move gigc_initialized_mod code here and move
!                              dlat/dlon calculation to gc_init_grid;
!                              GC timesteps are now seconds;
!                              Call set_input_opt to initialize input_opt vars;
!                              Add error handling using MAPL Assert_;
!                              Rename Initialize_Geos_Grid to GC_Init_Grid;
!                              Now call GC_Allocate_All after input.geos read;
!                              Restructure grid init based on gcbe v11-02e;
!                              Remove all unused code and simplify comments
!  14 Dec 2018 - E. Lundgren - Combine offline and GEOS-5 code and simplify
!  14 Jan 2019 - E. Lundgren - Read input.geos on all threads; remove broadcast
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                        :: I, J, L, STATUS
    CHARACTER(LEN=ESMF_MAXSTR)     :: Iam
    TYPE(ESMF_Config)              :: CF            ! Grid comp config object

#ifdef ADJOINT
    ! Adoint variables
    ! Local Finite Difference variables
    REAL(fp)                       :: FD_LAT, FD_LON
    INTEGER                        :: FD_STEP
    CHARACTER(LEN=ESMF_MAXSTR)     :: FD_SPEC
    REAL(fp)                       :: d, dmin
    INTEGER                        :: imin, jmin, NFD, LFD
    INTEGER                        :: IFD, JFD
    CHARACTER(LEN=ESMF_MAXSTR)     :: FD_TYPE

    ! At present, we are unable to load cube-sphere files through ExtData
    ! so we will define the cost function region thusly in GCHP.rc
    INTEGER                        :: CF_IMIN, CF_IMAX
    INTEGER                        :: CF_JMIN, CF_JMAX
    INTEGER                        :: CF_LMIN, CF_LMAX

    ! Need to get gloabl grid information for some FD spot tests
    TYPE(ESMF_Grid)                :: grid           ! ESMF Grid object
    INTEGER                        :: IL_PET, IU_PET ! Global lon bounds on this PET
    INTEGER                        :: JL_PET, JU_PET ! Global lat bounds on this PET

    ! Model phase: fwd, TLM, ADJOINT
    CHARACTER(LEN=ESMF_MAXSTR)     :: ModelPhase
#endif

    !=======================================================================
    ! GIGC_CHUNK_INIT begins here 
    !=======================================================================

    ! Error trap
    Iam = 'GIGC_CHUNK_INIT (gigc_chunk_mod.F90)'

    ! Assume success
    RC = GC_SUCCESS

#if defined( MODEL_GEOS )
    ! ckeller, 01/16/17
    ! ewl to do: revisit if we need these here
    Input_Opt%MAX_DIAG      = 1 
    Input_Opt%MAX_FAM       = 250
    Input_Opt%LINOZ_NLAT    = 18
    Input_Opt%LINOZ_NMONTHS = 12
    Input_Opt%LINOZ_NFIELDS = 7
    Input_Opt%RootCPU       = am_I_Root
#endif

    ! Get memory debug level
    call ESMF_GridCompGet ( GC, config=CF, RC=STATUS )
    _VERIFY(STATUS)
    call ESMF_ConfigGetAttribute(CF, MemDebugLevel, &
                                 Label="MEMORY_DEBUG_LEVEL:" , RC=STATUS)
    _VERIFY(STATUS)

    ! Initialize Input_Opt fields to zeros or equivalent
    CALL Set_Input_Opt( am_I_Root, Input_Opt, RC )
    _ASSERT(RC==GC_SUCCESS, 'informative message here')

    ! Read input.geos at very beginning of simulation on every thread
    CALL Read_Input_File( am_I_Root, Input_Opt, State_Grid, RC )
    _ASSERT(RC==GC_SUCCESS, 'informative message here')

    ! Initialize GEOS-Chem horizontal grid structure
    CALL GC_Init_Grid( am_I_Root, Input_Opt, State_Grid, RC )
    _ASSERT(RC==GC_SUCCESS, 'informative message here')

    ! Set maximum number of levels in the chemistry grid
    IF ( Input_Opt%LUCX ) THEN
       State_Grid%MaxChemLev  = State_Grid%MaxStratLev
    ELSE
       State_Grid%MaxChemLev  = State_Grid%MaxTropLev
    ENDIF

    ! In the ESMF/MPI environment, we can get the total overhead ozone
    ! either from the met fields (GIGCsa) or from the Import State (GEOS-5)
    Input_Opt%USE_O3_FROM_MET = .TRUE.
    
    ! Read LINOZ climatology
    IF ( Input_Opt%LLINOZ ) THEN
       CALL Linoz_Read( am_I_Root, Input_Opt, RC ) 
       _ASSERT(RC==GC_SUCCESS, 'informative message here')
    ENDIF

    ! Allocate all lat/lon arrays
    CALL GC_Allocate_All( am_I_Root, Input_Opt, State_Grid, RC )
    _ASSERT(RC==GC_SUCCESS, 'informative message here')

    ! Set grid based on passed mid-points
    CALL SetGridFromCtr( am_I_Root, State_Grid, lonCtr, latCtr, RC )
    _ASSERT(RC==GC_SUCCESS, 'informative message here')

    ! Update Input_Opt with timing fields
    Input_Opt%NYMDb   = nymdB
    Input_Opt%NHMSb   = nhmsB
    Input_Opt%NYMDe   = nymdE
    Input_Opt%NHMSe   = nhmsE
    Input_Opt%TS_CHEM = INT( tsChem )   ! Chemistry timestep [sec]
    Input_Opt%TS_EMIS = INT( tsChem )   ! Chemistry timestep [sec]
    Input_Opt%TS_DYN  = INT( tsDyn  )   ! Dynamic   timestep [sec]
    Input_Opt%TS_CONV = INT( tsDyn  )   ! Dynamic   timestep [sec]
    Input_Opt%myCPU   = myPET

    ! Set GEOS-Chem timesteps on all CPUs
    CALL SET_TIMESTEPS( am_I_Root  = am_I_Root,                          &
                        Chemistry  = Input_Opt%TS_CHEM,                  &
                        Convection = Input_Opt%TS_CONV,                  &
                        Dynamics   = Input_Opt%TS_DYN,                   &
                        Emission   = Input_Opt%TS_EMIS,                  &
                        Radiation  = Input_Opt%TS_RAD,                   &
                        Unit_Conv  = MAX( Input_Opt%TS_DYN,              &
                                          Input_Opt%TS_CONV ),           &
                        Diagnos    = Input_Opt%TS_DIAG         )

    ! Initialize derived-type objects for met, chem, and diag
    CALL GC_Init_StateObj( am_I_Root, HistoryConfig%DiagList, Input_Opt, &
                           State_Chm, State_Diag, State_Grid, State_Met, RC )
    _ASSERT(RC==GC_SUCCESS, 'informative message here')

#ifdef ADJOINT
    ! Are we running the adjoint?
    call ESMF_ConfigGetAttribute(CF, ModelPhase,            &
                                 Label="MODEL_PHASE:" ,         &
                                 Default="FORWARD",  RC=STATUS)
    _VERIFY(STATUS)
    call WRITE_PARALLEL('Checking if this is adjoint. Model phase = "' // trim(ModelPhase) // '"')
    input_opt%IS_ADJOINT = .FALSE.
    if (TRIM(ModelPhase) .eq. 'ADJOINT') THEN
       call WRITE_PARALLEL('Yes! Setting IS_ADJOINT to true.')
       input_opt%IS_ADJOINT = .TRUE.
    endif

    call ESMF_ConfigGetAttribute(CF, FD_TYPE, &
         Label="FD_TYPE:" , Default='NONE', RC=STATUS)
    _VERIFY(STATUS)

    Input_Opt%IS_FD_GLOBAL = TRIM(To_UpperCase(FD_TYPE(1:4))) == 'GLOB'
    Input_Opt%IS_FD_SPOT   = TRIM(To_UpperCase(FD_TYPE(1:4))) == 'SPOT'
    IF (MAPL_Am_I_Root()) THEN
       WRITE(*,1091) TRIM(FD_TYPE), Input_Opt%IS_FD_GLOBAL, Input_Opt%IS_FD_SPOT
    ENDIF
1091   FORMAT('FD_TYPE = ', a6, ', FD_GLOB = ', L1, ', FD_SPOT = ', L1)



    call ESMF_ConfigGetAttribute(CF, FD_STEP, &
         Label="FD_STEP:" , Default=-1, RC=STATUS)
    _VERIFY(STATUS)

    IF (Input_Opt%IS_FD_GLOBAL .or. Input_Opt%IS_FD_SPOT)  THEN
       _ASSERT(FD_STEP /= -1, 'FD_GLOB or FD_SPOT require FD_STEP')
    ENDIF


    if (.not. FD_STEP == -1 .or. input_opt%IS_ADJOINT) THEN
       Input_Opt%FD_STEP = FD_STEP

       call ESMF_ConfigGetAttribute(CF, FD_SPEC, &
            Label="FD_SPEC:", default="", RC=STATUS)
       _VERIFY(STATUS)
       IF (TRIM(FD_SPEC ) == "") THEN
          NFD = -1
       ELSE
          NFD = Ind_(FD_SPEC)
       ENDIF

       call ESMF_ConfigGetAttribute(CF, IFD, &
            Label="IFD:", default=-1, RC=STATUS)
       _VERIFY(STATUS)

       call ESMF_ConfigGetAttribute(CF, JFD, &
            Label="JFD:", default=-1, RC=STATUS)
       _VERIFY(STATUS)

       ! Get the ESMF grid attached to this gridded component
       CALL ESMF_GridCompGet( GC, grid=Grid, __RC__ )

       ! Get the upper and lower bounds of on each PET using MAPL
       CALL MAPL_GridGetInterior( Grid, IL_PET, IU_PET, JL_PET, JU_PET )

       ! See if we specified IFD and JFD in GCHP.rc
       IF ( IFD > 0 .and. JFD > 0 ) THEN

          if (IL_PET .le. IFD .and. IFD .le. IU_PET .and. &
               JL_PET .le. JFD .and. JFD .le. JU_PET) THEN
             Input_Opt%IS_FD_SPOT_THIS_PET = .true.
             Input_opt%IFD = IFD - IL_PET + 1
             Input_Opt%JFD = JFD - JL_PET + 1

             ! set these for debug printing
             DMIN = 0.0
             IMIN = Input_Opt%IFD
             JMIN = Input_Opt%JFD
          ENDIF

       ELSE

          call ESMF_ConfigGetAttribute(CF, FD_LAT, &
               Label="FD_LAT:", default=-999.0d0, RC=STATUS)
          _VERIFY(STATUS)

          call ESMF_ConfigGetAttribute(CF, FD_LON, &
               Label="FD_LON:", default=-999.0d0, RC=STATUS)
          _VERIFY(STATUS)

          _ASSERT( FD_LAT .ne. -999.0d0 .and. FD_LON .ne. -999.0d0, 'FD_SPOT requires either IFD and JFD or FD_LAT and FD_LON be set in GCHP.rc')


          dmin = 99999.9
          imin = -1
          jmin = -1
          ! try to find lat lon grid cell closest to 44.65, -63.58 (Halifax, NS)
          DO I = 1, state_grid%nx
             DO J = 1, state_grid%ny
                d = sqrt((state_grid%XMID(I,J) - FD_LON)**2 + &
                     (state_grid%YMID(I,J) - FD_LAT)**2)
                if (d < dmin) then
                   dmin = d
                   imin = i
                   jmin = j
                endif
             enddo
          enddo
          ! this is terrible. We need a better way to figure out if we're really in
          ! a grid cell, bbut I don't know how to do that. For now we're just hardcoding
          ! to the value for C24 and hoping for no points near cubed-sphere face
          ! boundaries
          if (dmin < 3.2) then
             ! getting the global grid offset is possible, see Chem_GridCompMod.F90:Extract_
             Input_Opt%IS_FD_SPOT_THIS_PET = .true.
             Input_Opt%IFD = IMIN
             Input_Opt%JFD = JMIN

          end if
       ENDIF

       Input_Opt%NFD = NFD

       call ESMF_ConfigGetAttribute(CF, LFD, &
            Label="LFD:", RC=STATUS)
       _VERIFY(STATUS)

       Input_Opt%LFD = LFD

       ! Read in cost function region

       call ESMF_ConfigGetAttribute(CF, CF_IMIN, &
            Label="CF_IMIN:", default=-1, RC=STATUS)
       _VERIFY(STATUS)

       CF_IMIN = CF_IMIN - IL_PET + 1

       call ESMF_ConfigGetAttribute(CF, CF_IMAX, &
            Label="CF_IMAX:", default=-1, RC=STATUS)
       _VERIFY(STATUS)

       CF_IMAX = CF_IMAX - IL_PET + 1

       call ESMF_ConfigGetAttribute(CF, CF_JMIN, &
            Label="CF_JMIN:", default=-1, RC=STATUS)
       _VERIFY(STATUS)

       CF_JMIN = CF_JMIN - JL_PET + 1

       call ESMF_ConfigGetAttribute(CF, CF_JMAX, &
            Label="CF_JMAX:", default=-1, RC=STATUS)
       _VERIFY(STATUS)

       CF_JMAX = CF_JMAX - JL_PET + 1

       call ESMF_ConfigGetAttribute(CF, CF_LMIN, &
            Label="CF_LMIN:", default=-1, RC=STATUS)
       _VERIFY(STATUS)

       call ESMF_ConfigGetAttribute(CF, CF_LMAX, &
            Label="CF_LMAX:", default=-1, RC=STATUS)
       _VERIFY(STATUS)

       IF (CF_IMIN < 1 .OR. CF_IMIN > State_Grid%NX .OR. &
            CF_IMAX < 1 .OR. CF_IMAX > State_Grid%NX .OR. &
            CF_JMIN < 1 .OR. CF_JMIN > State_Grid%NY .OR. &
            CF_JMAX < 1 .OR. CF_JMAX > State_Grid%NY) THEN
       WRITE(*,1028) Input_Opt%MyCPU,   &
            Input_Opt%CF_IMIN, Input_Opt%CF_IMAX, &
            Input_Opt%CF_JMIN, Input_Opt%CF_JMAX, &
            Input_Opt%CF_LMIN, Input_Opt%CF_LMAX
1028   FORMAT('Pre-CF on Pet ', i3, ' I = (', i3, ', ', i3, ') &
             J = ( ', i3, ', ', i3, ') &
             L = (', i3, ', ', i3, ')')

          CF_IMIN = -1
          CF_IMAX = -1
          CF_JMIN = -1
          CF_JMAX = -1
          CF_LMIN = -1
          CF_LMAX = -1
       ENDIF

       _ASSERT(CF_IMIN * CF_IMAX > 0, 'Please define both max and min for CF_I')
       _ASSERT(CF_JMIN * CF_JMAX > 0, 'Please define both max and min for CF_J')
       _ASSERT(CF_LMIN * CF_LMAX > 0, 'Please define both max and min for CF_L')

       _ASSERT(CF_LMIN * CF_IMIN > 0, 'If CF_I: is defined, please define CF_L')
       _ASSERT(CF_JMIN * CF_IMIN > 0, 'If CF_I: is defined, please define CF_J')
       
       ! At this point, they should all be set or all be negative (probably -1)
       IF (CF_IMIN > 0) THEN
          Input_Opt%CF_IMIN = CF_IMIN
          Input_Opt%CF_IMAX = CF_IMAX
          Input_Opt%CF_JMIN = CF_JMIN
          Input_Opt%CF_JMAX = CF_JMAX
          Input_Opt%CF_LMIN = CF_LMIN
          Input_Opt%CF_LMAX = CF_LMAX
       ELSEIF (Input_Opt%IS_FD_SPOT_THIS_PET) THEN
          Input_Opt%CF_IMIN = Input_Opt%IFD
          Input_Opt%CF_IMAX = Input_Opt%IFD
          Input_Opt%CF_JMIN = Input_Opt%JFD
          Input_Opt%CF_JMAX = Input_Opt%JFD
          Input_Opt%CF_LMIN = Input_Opt%LFD
          Input_Opt%CF_LMAX = Input_Opt%LFD
       WRITE(*,1027) Input_Opt%MyCPU,   &
            Input_Opt%CF_IMIN, Input_Opt%CF_IMAX, &
            Input_Opt%CF_JMIN, Input_Opt%CF_JMAX, &
            Input_Opt%CF_LMIN, Input_Opt%CF_LMAX
1027   FORMAT('CF on Pet ', i3, ' I = (', i3, ', ', i3, ') &
             J = ( ', i3, ', ', i3, ') &
             L = (', i3, ', ', i3, ')')

       ELSE
          Input_Opt%CF_IMIN = -1
          Input_Opt%CF_IMAX = -1
          Input_Opt%CF_JMIN = -1
          Input_Opt%CF_JMAX = -1
          Input_Opt%CF_LMIN = -1
          Input_Opt%CF_LMAX = -1
       ENDIF

       IF ( Input_Opt%IS_FD_SPOT_THIS_PET ) THEN
          write (*,1011) myPet, dmin, imin, jmin, &
               state_grid%YMID(IMIN,JMIN), state_grid%XMID(IMIN,JMIN)

#ifdef DEBUG
          ! Get the ESMF grid attached to this gridded component
          CALL ESMF_GridCompGet( GC, grid=Grid, __RC__ )

          ! Get the upper and lower bounds of on each PET using MAPL
          CALL MAPL_GridGetInterior( Grid, IL_PET, IU_PET, JL_PET, JU_PET )
          WRITE(*,1013) IL_PET, IU_PET
          WRITE(*,1014) JL_PET, JU_PET
#endif
          
         ENDIF

1011   FORMAT('Found FD_SPOT on PET ', i5, ' ', f7.2, &
            ' degrees from cell ', i3, ', ', i3, ' (', f7.2, ', ', f7.2, ')')
1012   FORMAT('Did not find FD_SPOT on PET ', i5, ' ', f7.2,&
            ' degrees from cell ', i3, ', ', i3, ' (', f7.2, ', ', f7.2, ')')
1013   FORMAT('   XminOffset = ', i3, '     XmaxOffset = ', i3)
1014   FORMAT('   YminOffset = ', i3, '     YmaxOffset = ', i3)
1015   FORMAT('   GlobalXMid(', i3, ', ', i3, ') = (', f7.2, ', ' f7.2, ')')
1016   FORMAT('       SPC(', a10, ', FD_SPOT) = ', e22.10)
1019   FORMAT('   SPC_ADJ(', a10, ', FD_SPOT) = ', e22.10)
    ENDIF
#endif

    ! Initialize other GEOS-Chem modules
    CALL GC_Init_Extra( am_I_Root, HistoryConfig%DiagList, Input_Opt,    &
                        State_Chm, State_Diag, State_Grid, RC ) 
    _ASSERT(RC==GC_SUCCESS, 'informative message here')

    ! Set initial State_Chm%Species units to units expected in transport
# if defined( MODEL_GEOS )
    State_Chm%Spc_Units = 'kg/kg total'
#else
    State_Chm%Spc_Units = 'kg/kg dry'
#endif

    ! Initialize the GEOS-Chem pressure module (set Ap & Bp)
    CALL Init_Pressure( am_I_Root, State_Grid, RC )

    ! Initialize the PBL mixing module
    CALL INIT_PBL_MIX( am_I_Root, State_Grid, RC )
    _ASSERT(RC==GC_SUCCESS, 'informative message here')

    ! Initialize chemistry mechanism
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM .OR. Input_Opt%ITS_AN_AEROSOL_SIM ) THEN
       CALL INIT_CHEMISTRY ( am_I_Root,  Input_Opt,  State_Chm, &
                             State_Diag, State_Grid, RC )
       _ASSERT(RC==GC_SUCCESS, 'informative message here')
    ENDIF

    ! Initialize HEMCO
    CALL EMISSIONS_INIT ( am_I_Root,  Input_Opt, State_Chm, &
                          State_Grid, State_Met, RC, &
                          HcoConfig=HcoConfig )
    _ASSERT(RC==GC_SUCCESS, 'informative message here')

    ! Stratosphere - can't be initialized without HEMCO because of STATE_PSC
    IF ( Input_Opt%LUCX ) THEN

       ! Initialize stratospheric routines
       CALL INIT_UCX( am_I_Root, Input_Opt, State_Chm, State_Diag, State_Grid )

    ENDIF

#if defined( MODEL_GEOS )
    ! Keep commented out line with LLSTRAT as a GEOS-5 option reminder
    !IF ( Input_Opt%LSCHEM .AND. Input_Opt%LLSTRAT < value_LM ) THEN
#endif
     IF ( Input_Opt%LSCHEM ) THEN
       CALL INIT_STRAT_CHEM( am_I_Root, Input_Opt,  State_Chm, & 
                             State_Met, State_Grid, RC )
       IF (RC /= GC_SUCCESS) RETURN
    ENDIF

    !-------------------------------------------------------------------------
    ! Diagnostics and tendencies 
    !-------------------------------------------------------------------------

#if defined( MODEL_GEOS )
    ! The GEOS-Chem diagnostics list, stored in HistoryConfig, is initialized 
    ! during GIGC_INIT_SIMULATION, and corresponding arrays in State_Diag are 
    ! allocated accordingly when initializing State_Diag. Here, we thus 
    ! only need to initialize the tendencies, which have not been initialized
    ! yet (ckeller, 11/29/17). 
    CALL Tend_Init ( am_I_Root, Input_Opt, State_Met, State_Chm, RC ) 
    _ASSERT(RC==GC_SUCCESS, 'informative message here')
#endif

#if !defined( MODEL_GEOS )
    ! GCHP only: Convert species units to internal state units (v/v dry)
    CALL Convert_Spc_Units( am_I_Root,  Input_Opt, State_Chm, &
                            State_Grid, State_Met, 'v/v dry', RC )
    _ASSERT(RC==GC_SUCCESS, 'informative message here')
#endif

    ! Return success
    RC = GC_Success

  END SUBROUTINE GIGC_Chunk_Init
!EOC

!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gigc_chunk_run
!
! !DESCRIPTION: Subroutine GIGC\_CHUNK\_RUN is the ESMF run method for
!  GEOS-Chem.
!
! !INTERFACE:
!
  SUBROUTINE GIGC_Chunk_Run( am_I_Root,  GC,                                 &
                             nymd,       nhms,       year,       month,      &
                             day,        dayOfYr,    hour,       minute,     &
                             second,     utc,        hElapsed,   Input_Opt,  &
                             State_Chm,  State_Diag, State_Grid, State_Met,  &
                             Phase,      IsChemTime,                         &
#if defined( MODEL_GEOS )
                             FrstRewind, &
#endif
#if defined( ADJOINT )
                             IsStarttime, &
#endif
                             RC )
!
! !USES:
!
    ! GEOS-Chem state objects 
    USE HCO_Interface_Mod,  ONLY : HcoState
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState

    ! GEOS-Chem components
    USE Chemistry_Mod,      ONLY : Do_Chemistry, Recompute_OD
    USE Convection_Mod,     ONLY : Do_Convection
    USE DryDep_Mod,         ONLY : Do_DryDep
    USE Emissions_Mod,      ONLY : Emissions_Run
    USE Mixing_Mod,         ONLY : Do_Tend, Do_Mixing
    USE WetScav_Mod,        ONLY : Setup_WetScav, Do_WetDep

    ! Specialized subroutines
    USE Dao_Mod,            ONLY : AirQnt, Set_Dry_Surface_Pressure
    USE Dao_Mod,            ONLY : GIGC_Cap_Tropopause_Prs
    USE Set_Global_CH4_Mod, ONLY : Set_CH4
    USE MODIS_LAI_Mod,      ONLY : Compute_XLAI
    USE PBL_Mix_Mod,        ONLY : Compute_PBL_Height
    USE Pressure_Mod,       ONLY : Set_Floating_Pressures
    USE TOMS_Mod,           ONLY : Compute_Overhead_O3
    USE UCX_Mod,            ONLY : Set_H2O_Trac

    ! Utilities
    USE ErrCode_Mod
    USE HCO_Error_Mod
    USE HCO_Interface_Mod,  ONLY : SetHcoTime
    USE MAPL_MemUtilsMod
    USE Pressure_Mod,       ONLY : Accept_External_Pedge
    USE State_Chm_Mod,      ONLY : IND_
    USE Time_Mod,           ONLY : Accept_External_Date_Time
    USE UnitConv_Mod,       ONLY : Convert_Spc_Units, Print_Global_Species_Kg

    ! Diagnostics
    USE Diagnostics_Mod,    ONLY : Set_Diagnostics_EndofTimestep, Set_SpcAdj_Diagnostic
    USE Aerosol_Mod,        ONLY : Set_AerMass_Diagnostic

#if defined( MODEL_GEOS )
    USE CMN_SIZE_MOD,       ONLY : IIPAR, JJPAR, LLPAR
    USE DAO_MOD,            ONLY : GET_COSINE_SZA
    USE DIAG_MOD,           ONLY : AD21
    USE HCOI_GC_MAIN_MOD,   ONLY : HCOI_GC_WriteDiagn
#endif
    USE Species_Mod,   ONLY : Species

!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on root CPU?
    INTEGER,        INTENT(IN)    :: nymd        ! YYYY/MM/DD @ current time
    INTEGER,        INTENT(IN)    :: nhms        ! hh:mm:ss   @ current time
    INTEGER,        INTENT(IN)    :: year        ! UTC year 
    INTEGER,        INTENT(IN)    :: month       ! UTC month
    INTEGER,        INTENT(IN)    :: day         ! UTC day
    INTEGER,        INTENT(IN)    :: dayOfYr     ! UTC day of year
    INTEGER,        INTENT(IN)    :: hour        ! UTC hour
    INTEGER,        INTENT(IN)    :: minute      ! UTC minute
    INTEGER,        INTENT(IN)    :: second      ! UTC second
    REAL*4,         INTENT(IN)    :: utc         ! UTC time [hrs]
    REAL*4,         INTENT(IN)    :: hElapsed    ! Elapsed hours
    INTEGER,        INTENT(IN)    :: Phase       ! Run phase (-1, 1 or 2)
    LOGICAL,        INTENT(IN)    :: IsChemTime  ! Time for chemistry? 
#if defined( MODEL_GEOS )
    LOGICAL,        INTENT(IN)    :: FrstRewind  ! Is it the first rewind? 
#endif
#if defined ( ADJOINT )
    LOGICAL,        INTENT(IN)    :: IsStarttime ! Have we reached the start time
                                                 ! in an adjoint run
#endif 
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ESMF_GridComp), INTENT(INOUT) :: GC          ! Ref to this GridComp
    TYPE(OptInput),      INTENT(INOUT) :: Input_Opt   ! Input Options obj
    TYPE(ChmState),      INTENT(INOUT) :: State_Chm   ! Chemistry State obj
    TYPE(DgnState),      INTENT(INOUT) :: State_Diag  ! Diagnostics State obj
    TYPE(GrdState),      INTENT(INOUT) :: State_Grid  ! Grid State obj
    TYPE(MetState),      INTENT(INOUT) :: State_Met   ! Meteorology State obj
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT)   :: RC          ! Return code
!
! !REMARKS:
!
! !REVISION HISTORY:
!  18 Jul 2011 - M. Long     - Initial Version
!  09 Oct 2012 - R. Yantosca - Added extra comments & cosmetic changes
!  16 Oct 2012 - R. Yantosca - Renamed GC_STATE argument to State_Chm
!  16 Oct 2012 - R. Yantosca - Renamed GC_MET argument to State_Met
!  17 Oct 2012 - R. Yantosca - Need to call AIRQNT before chemistry
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_chm_mod.F90
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_met_mod.F90
!  22 Oct 2012 - R. Yantosca - Renamed to GIGC_Chunk_Run
!  25 Oct 2012 - R. Yantosca - Now pass RC to GIGC_DO_CHEM
!  01 Nov 2012 - R. Yantosca - Now reference gigc_input_opt_mod.F90
!  08 Nov 2012 - R. Yantosca - Now pass Input_Opt to GIGC_Do_Chem
!  13 Nov 2012 - M. Long     - Added Dry Deposition method
!  29 Nov 2012 - R. Yantosca - Now block off calls to GIGC_DO_DRYDEP and
!                              GIGC_DO_CHEM w/ the appropriate logical flags
!  04 Dec 2012 - R. Yantosca - Now convert units of State_Chm%TRACERS here
!                              instead of in lower-level routines
!  07 Dec 2012 - R. Yantosca - Now call Accept_Date_Time_From_ESMF to pass the
!                              date & time from ESMF to GeosUtil/time_mod.F
!  07 Dec 2012 - R. Yantosca - Now pass UTC via Accept_Date_Time_From_ESMF;
!                              this ensures proper localtime computation
!  11 Dec 2012 - R. Yantosca - Now call DO_DRYDEP directly; no longer call
!                              GIGC_DO_DRYDEP, this is moved to obsolete dir.
!  11 Dec 2012 - R. Yantosca - Now call routine ACCEPT_EXTERNAL_PEDGE to pass
!                              the pressure edges from ESMF to GEOS-Chem
!  15 Mar 2013 - R. Yantosca - Now call GIGC_CAP_TROPOPAUSE_PRS to cap the
!                              State_Met%TROPP field to 200 hPa polewards
!                              of 60S and 60N.  We do this in the std G-C.
!  05 Jun 2013 - R. Yantosca - Remove obsolete code
!  22 Sep 2014 - C. Keller   - Added run phase argument
!  14 Oct 2014 - C. Keller   - Various updates to include drydep and emissions
!                              to tracer arrays, etc.
!  26 Nov 2014 - C. Keller   - Added IsChemTime variable.
!  19 Oct 2016 - R. Yantosca - Now call Set_Init_Cond_Strat_Chem after the
!                              1st call to AIRQNT to save initial conditions
!  01 Dec 2016 - E. Lundgren - Calculate LAI using new routine for GCHP
!  13 Feb 2018 - E. Lundgren - Call Recompute_OD at end of chem dt for aer diags
!  14 Dec 2018 - E. Lundgren - Combine offline and GEOS-5 code and simplify
!EOP
!------------------------------------------------------------------------------
!BOC
    TYPE(MAPL_MetaComp), POINTER   :: STATE
    TYPE(ESMF_VM)                  :: VM            ! ESMF VM object
    REAL*8                         :: DT
    CHARACTER(LEN=ESMF_MAXSTR)     :: Iam, OrigUnit
    INTEGER                        :: STATUS, HCO_PHASE
#if defined( MODEL_GEOS )
    INTEGER                        :: N, I, J, L
#endif

    ! Local logicals to turn on/off individual components
    ! The parts to be executed are based on the input options,
    ! the time step and the phase.
    LOGICAL                        :: DoConv 
    LOGICAL                        :: DoDryDep
    LOGICAL                        :: DoEmis
    LOGICAL                        :: DoTend 
    LOGICAL                        :: DoTurb 
    LOGICAL                        :: DoChem
    LOGICAL                        :: DoWetDep

    ! First call?
    LOGICAL, SAVE                  :: FIRST = .TRUE.

    ! # of times this routine has been called. Only temporary for printing 
    ! processes on the first 10 calls.
    INTEGER, SAVE                  :: NCALLS = 0

    ! Strat. H2O settings 
    LOGICAL                        :: SetStratH2O 
#if defined( MODEL_GEOS )
    LOGICAL, SAVE                  :: LSETH2O_orig
#endif

    ! Debug variables
    INTEGER, parameter             :: I_DBG = 6, J_DBG = 5, L_DBG=1
#ifdef ADJOINT
    ! Adjoint Finitie Difference Variables
    INTEGER                        :: IFD, JFD, LFD, NFD
    INTEGER                        :: I, J, L
    REAL*8                         :: CFN
#endif
#if !defined( MODEL_GEOS )
    INTEGER                        :: N
#endif
#if defined ( ADJOINT )
    CHARACTER(len=ESMF_MAXSTR)     :: FD_SPEC, TRACNAME
    TYPE(Species),       POINTER   :: ThisSpc
#endif

    !=======================================================================
    ! GIGC_CHUNK_RUN begins here 
    !=======================================================================

    ! Error trap
    Iam = 'GIGC_CHUNK_RUN (gigc_chunk_mod.F90)'

    ! Assume success
    RC = GC_SUCCESS

    ! Get state object (needed for timers)
    CALL MAPL_GetObjectFromGC(GC, STATE, __RC__)

    ! Get the VM for optional memory prints (level >= 2)
    !-----------------------------------
    if ( MemDebugLevel > 0 ) THEN
       call ESMF_VmGetCurrent(VM, RC=STATUS)
       _VERIFY(STATUS)
    endif

    !=======================================================================
    ! Define processes to be covered in this phase
    !
    ! In the standard GEOS-Chem, the following operator sequence is used:
    ! 1. DryDep (kg)
    ! 2. Emissions (kg)
    ! 3. Turbulence (v/v)
    ! 4. Convection (v/v)
    ! 5. Chemistry (kg)
    ! 6. Wetdep (kg)
    !
    ! The GEOS-5 operator sequence is:
    ! 1. Gravity wave drag
    ! 2. Moist (convection)
    ! 3. Chemistry 1 (drydep and emissions)
    ! 4. Surface 1
    ! 5. Turbulence 1
    ! 6. Surface 2
    ! 7. Turbulence 2
    ! 8. Chemistry 2 (chemistry and wet deposition)
    ! 9. Radiation 
    !
    ! Here, we use the following operator sequence:
    ! 
    ! 1.  Convection (v/v) --> Phase 1
    ! 2.  DryDep (kg)      --> Phase 1
    ! 3.  Emissions (kg)   --> Phase 1
    ! 4a. Tendencies (v/v) --> Phase 1
    ! -------------------------------
    ! 4b. Turbulence (v/v) --> Phase 2 
    ! 5.  Chemistry (kg)   --> Phase 2
    ! 6.  WetDep (kg)      --> Phase 2     
    ! 
    ! Any of the listed processes is only executed if the corresponding switch
    ! in the input.geos file is enabled. If the physics component already
    ! covers convection or turbulence, they should not be applied here!
    ! The tendencies are only applied if turbulence is not done within
    ! GEOS-Chem (ckeller, 10/14/14).
    ! 
    ! The standard number of phases in GCHP is 1, set in GCHP.rc, which
    ! results in Phase -1 in gigc_chunk_run. This results in executing
    ! all GEOS-Chem components in a single run rather than splitting up
    ! across two runs as is done in GEOS-5. (ewl, 10/26/18)
    !=======================================================================

    ! By default, do processes as defined in input.geos. DoTend defined below. 
    DoConv   = Input_Opt%LCONV                    ! dynamic time step
    DoDryDep = Input_Opt%LDRYD .AND. IsChemTime   ! chemistry time step
    DoEmis   = Input_Opt%LEMIS .AND. IsChemTime   ! chemistry time step
#if defined( MODEL_GEOS )
    DoTurb   = Input_Opt%LTURB .AND. IsChemTime   ! dynamic time step
#else
    DoTurb   = Input_Opt%LTURB                    ! dynamic time step
#endif
    DoChem   = Input_Opt%LCHEM .AND. IsChemTime   ! chemistry time step
    DoWetDep = Input_Opt%LWETD                    ! dynamic time step 

    ! If Phase is not -1, only do selected processes for given phases: 
    ! Phase 1: disable turbulence, chemistry and wet deposition.
    IF ( Phase == 1 ) THEN
       DoTurb   = .FALSE.
       DoChem   = .FALSE.
       DoWetDep = .FALSE.

    ! Phase 2: disable convection, drydep and emissions. 
    ELSEIF ( Phase == 2 ) THEN
       DoConv   = .FALSE.
       DoDryDep = .FALSE.
       DoEmis   = .FALSE. 
    ENDIF

    ! Check if tendencies need be applied. The drydep and emission calls
    ! only calculates the emission / drydep rates, but do not apply the
    ! tendencies to the tracer array yet. If turbulence is done as part of
    ! GEOS-5, we need to make sure that these tendencies are applied to the
    ! tracer array. If turbulence is explicitly covered by GEOS-Chem,
    ! however, the tendencies become automatically applied within the PBL
    ! mixing routines (DO_MIXING), so we should never apply the tendencies
    ! in this case.
    DoTend = ( DoEmis .OR. DoDryDep ) .AND. .NOT. Input_Opt%LTURB

    ! testing only
    IF ( am_I_Root .and. NCALLS < 10 ) THEN 
       write(*,*) 'GEOS-Chem phase ', Phase, ':'
       write(*,*) 'DoConv   : ', DoConv
       write(*,*) 'DoDryDep : ', DoDryDep
       write(*,*) 'DoEmis   : ', DoEmis
       write(*,*) 'DoTend   : ', DoTend
       write(*,*) 'DoTurb   : ', DoTurb
       write(*,*) 'DoChem   : ', DoChem
       write(*,*) 'DoWetDep : ', DoWetDep
       write(*,*) ' '
    ENDIF

    !-------------------------------------------------------------------------
    ! Pre-Run assignments
    !-------------------------------------------------------------------------

    ! Eventually initialize/reset wetdep
    IF ( DoConv .OR. DoChem .OR. DoWetDep ) THEN
       CALL SETUP_WETSCAV( am_I_Root,  Input_Opt, State_Chm, &
                           State_Grid, State_Met, RC )
    ENDIF

    ! Pass time values obtained from the ESMF environment to GEOS-Chem
    CALL Accept_External_Date_Time( am_I_Root      = am_I_Root,  &
                                    value_NYMD     = nymd,       &  
                                    value_NHMS     = nhms,       &  
                                    value_YEAR     = year,       &  
                                    value_MONTH    = month,      &  
                                    value_DAY      = day,        &  
                                    value_DAYOFYR  = dayOfYr,    &  
                                    value_HOUR     = hour,       &  
                                    value_MINUTE   = minute,     &  
                                    value_HELAPSED = hElapsed,   & 
                                    value_UTC      = utc,        &
                                    RC             = RC         )

    ! Set HEMCO time
    CALL SetHcoTime ( am_I_Root, DoEmis, RC )

    ! Calculate MODIS leaf area indexes needed for dry deposition
    CALL Compute_XLAI( am_I_Root, Input_Opt, State_Grid, State_Met, RC )

    ! Set the pressure at level edges [hPa] from the ESMF environment
    CALL Accept_External_Pedge    ( am_I_Root      = am_I_Root,  &
                                    State_Met      = State_Met,  &
                                    RC             = RC         )

#if defined( MODEL_GEOS )
    ! Eventually set tropopause pressure according to ozone values
    ! (use ozonopause) 
    CALL SET_OZONOPAUSE ( am_I_Root,  Input_Opt, State_Chm, &
                          State_Grid, State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN 
#endif

    ! Set dry surface pressure (PS1_DRY) from State_Met%PS1_WET
    CALL SET_DRY_SURFACE_PRESSURE( State_Grid, State_Met, 1 )

    ! Set dry surface pressure (PS2_DRY) from State_Met%PS2_WET
    CALL SET_DRY_SURFACE_PRESSURE( State_Grid, State_Met, 2 )

    ! Initialize surface pressures to match the post-advection pressures
    State_Met%PSC2_WET = State_Met%PS1_WET
    State_Met%PSC2_DRY = State_Met%PS1_DRY
    CALL SET_FLOATING_PRESSURES( am_I_Root, State_Grid, State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Define airmass and related quantities. Do not scale mixing ratio
    ! since mass conservation across timesteps is handled in FV3 advection.
    ! Beware that this means tracer mass will not be conserved across timesteps
    ! if advection is turned off.
    CALL AirQnt( am_I_Root, Input_opt, State_Chm, State_Grid, &
                 State_Met, RC, .FALSE. )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Cap the polar tropopause pressures at 200 hPa, in order to avoid
    ! tropospheric chemistry from happening too high up (cf. J. Logan)
    CALL GIGC_Cap_Tropopause_Prs  ( am_I_Root      = am_I_Root,  &
                                    Input_Opt      = Input_Opt,  &
                                    State_Grid     = State_Grid, &
                                    State_Met      = State_Met,  &
                                    RC             = RC         )

    ! Call PBL quantities. Those are always needed
    CALL COMPUTE_PBL_HEIGHT( am_I_Root, State_Grid, State_Met, RC )

    !=======================================================================
    ! Convert State_Chm%Species units to kg/kg dry
    !=======================================================================
#if defined( MODEL_GEOS )
    ! Adjust total mixing ratio to account for change in specific
    ! humidity and dry air mass since end of last timestep. Always skip
    ! first step and do not apply if transport is turned on. May need to
    ! further adjust conditional if this does not catch all cases (ewl, 11/8/18)
    IF ( (.NOT. FIRST) .AND. ( .NOT. Input_Opt%LTRAN ) ) THEN
       DO N = 1, State_Chm%nSpecies
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX
         State_Chm%Species(I,J,L,N) = State_Chm%Species(I,J,L,N)          &
                 / ( 1e0_fp - ( State_Met%SPHU_PREV(I,J,L) * 1e-3_fp ) )  &
                 * ( 1e0_fp - ( State_Met%SPHU(I,J,L) * 1e-3_fp ) )       &
                 * State_Met%DP_DRY_PREV(I,J,L) / State_Met%DELP_DRY(I,J,L)    
       ENDDO
       ENDDO
       ENDDO
       ENDDO
    ENDIF
#endif
#ifdef ADJOINT
    if (.not. first) &
         CALL Print_Global_Species_Kg( am_I_Root, I_DBG, J_DBG, L_DBG, &
                                       'CO2', Input_Opt, State_Chm,   &
                                       State_Grid, State_Met, trim(Iam) // &
                                       ' before first unit conversion', RC)
    CALL GIGC_PRINT_MET( am_I_root, I_DBG, J_DBG, L_DBG, Input_Opt,&
         State_Grid, State_Met, trim(Iam) // ' before first unit conversion.', RC)
    
    IF (first .and. Input_Opt%IS_FD_SPOT_THIS_PET .and.  Input_Opt%IS_FD_SPOT) THEN
       FD_SPEC = transfer(state_chm%SpcData(Input_Opt%NFD)%Info%Name, FD_SPEC)
       IFD = Input_Opt%IFD
       JFD = Input_Opt%JFD
       LFD = Input_Opt%LFD
       NFD = Input_Opt%NFD
       WRITE (*, 1017) TRIM(FD_SPEC), state_chm%species(IFD, JFD, LFD, NFD)
       IF (Input_Opt%IS_ADJOINT) THEN
          WRITE(*,*) ' Computing final cost function'
          CFN = 0d0
          state_chm%SpeciesAdj(:,:,:,NFD) = 0d0
          DO L = 1,State_Grid%NZ
          DO J = 1,State_Grid%NY
          DO I = 1,State_Grid%NX
             if (State_chm%CostFuncMask(I,J,L) > 0d0) THEN
                WRITE (*, 1047) I, J, L, state_chm%species(I, J, L, NFD)
                state_chm%SpeciesAdj(I,J,L, NFD) = 1.0d0
                CFN = CFN + state_chm%species(I,J,L,NFD)
             endif
          ENDDO
          ENDDO
          ENDDO
          WRITE(*,'(a7, e22.10)') ' CFN = ', CFN
1047      FORMAT('  SPC(', i2, ', ', i2, ', ', i2, ') = ', e22.10)
       ELSE
          IF (Input_Opt%FD_STEP .eq. 0) THEN
             WRITE(*, *) '    Not perturbing'
          ELSEIF (Input_Opt%FD_STEP .eq. 1) THEN
             WRITE(*, *) '    Perturbing +0.1'
             state_chm%species(IFD, JFD, LFD, NFD) = state_chm%species(IFD, JFD, LFD, NFD) * 1.1d0
          ELSEIF (Input_Opt%FD_STEP .eq. 2) THEN
             WRITE(*, *) '    Perturbing -0.1'
             state_chm%species(IFD, JFD, LFD, NFD) = state_chm%species(IFD, JFD, LFD, NFD) * 0.9d0
          ELSE
             WRITE(*, *) '    FD_STEP = ', Input_Opt%FD_STEP, ' NOT SUPPORTED!'
          ENDIF
          WRITE (*, 1017) TRIM(FD_SPEC), state_chm%species(IFD, JFD, LFD, NFD)
       ENDIF
    ENDIF

    IF (first .and. Input_Opt%IS_FD_GLOBAL) THEN
       FD_SPEC = transfer(state_chm%SpcData(Input_Opt%NFD)%Info%Name, FD_SPEC)
       NFD = Input_Opt%NFD
       IF (Input_Opt%IS_FD_SPOT_THIS_PET) THEN
          IFD = Input_Opt%IFD
          JFD = Input_Opt%JFD
          LFD = Input_Opt%LFD
          WRITE (*, 1017) TRIM(FD_SPEC), state_chm%species(IFD, JFD, LFD, NFD)
          IF (Input_Opt%Is_Adjoint) &
               WRITE (*, 1018) TRIM(FD_SPEC), state_chm%SpeciesAdj(IFD, JFD, LFD, NFD)
       ENDIF
       IF (.not. Input_Opt%IS_ADJOINT) THEN
          IF (Input_Opt%FD_STEP .eq. 0) THEN
             WRITE(*, *) '    Not perturbing'
          ELSEIF (Input_Opt%FD_STEP .eq. 1) THEN
             WRITE(*, *) '    Perturbing +0.1'
             state_chm%species(:, :, :, NFD) = state_chm%species(:, :, :, NFD) * 1.1d0
          ELSEIF (Input_Opt%FD_STEP .eq. 2) THEN
             WRITE(*, *) '    Perturbing -0.1'
             state_chm%species(:, :, :, NFD) = state_chm%species(:, :, :, NFD) * 0.9d0
          ELSE
             WRITE(*, *) '    FD_STEP = ', Input_Opt%FD_STEP, ' NOT SUPPORTED!'
          ENDIF
          IF (Input_Opt%IS_FD_SPOT_THIS_PET) &
               WRITE (*, 1017) TRIM(FD_SPEC), state_chm%species(IFD, JFD, LFD, NFD)
       ELSE
          IF (NFD > 0) THEN
             state_chm%SpeciesAdj(:,:,:,NFD) = 1d0
          ENDIF
       ENDIF
    ENDIF

1017 FORMAT('       SPC(', a10, ', FD_SPOT) = ', e22.10)
1018   FORMAT('   SPC_ADJ(', a10, ', FD_SPOT) = ', e22.10)
    IF (Input_Opt%IS_FD_SPOT_THIS_PET ) THEN
       FD_SPEC = transfer(state_chm%SpcData(Input_Opt%NFD)%Info%Name, FD_SPEC)
       NFD = Input_Opt%NFD
       IFD = Input_Opt%IFD
       JFD = Input_Opt%JFD
       LFD = Input_Opt%LFD
       WRITE(*,1017) TRIM(FD_SPEC), state_chm%species(IFD, JFD, LFD, NFD)
       IF (Input_Opt%Is_Adjoint) &
            WRITE (*, 1018) TRIM(FD_SPEC), state_chm%SpeciesAdj(IFD, JFD, LFD, NFD)
    ENDIF
#endif
    
    ! Convert to dry mixing ratio
    CALL Convert_Spc_Units ( am_I_Root, Input_Opt, State_Chm, State_Grid, &
                             State_Met, 'kg/kg dry', RC, OrigUnit=OrigUnit )

    ! SDE 05/28/13: Set H2O to STT if relevant
    IF ( IND_('H2O','A') > 0 ) THEN
       SetStratH2O = .FALSE.
#if defined( MODEL_GEOS )
       !=======================================================================
       ! Tropospheric H2O is always prescribed (using GEOS Q). For strat H2O
       ! there are three options, controlled by toggles 'set initial global MR'
       ! in input.geos and 'Prescribe_strat_H2O' in GEOSCHEMchem_GridComp.rc: 
       ! (A) never prescribe strat H2O -> both toggles off
       ! (B) prescribe strat H2O on init time step -> toggle in input.goes on
       ! (C) always prescribe strat H2O -> toggle in GEOSCHEMchem_GridComp.rc on
       !=======================================================================
       IF ( FIRST ) THEN
          LSETH2O_orig = Input_Opt%LSETH2O
       ENDIF
       IF ( FIRST .OR. FrstRewind ) THEN
          Input_Opt%LSETH2O = LSETH2O_orig
       ENDIF
       IF ( Input_Opt%LSETH2O .OR. .NOT. Input_Opt%LUCX .OR. &
            Input_Opt%AlwaysSetH2O ) THEN
#else
       IF ( Input_Opt%LSETH2O .OR. .NOT. Input_Opt%LUCX ) THEN
#endif
          SetStratH2O = .TRUE.
       ENDIF
       CALL SET_H2O_TRAC( am_I_Root, SetStratH2O, Input_Opt, & 
                          State_Chm, State_Grid,  State_Met, RC )
       _ASSERT(RC==GC_SUCCESS, 'informative message here')

      ! Only force strat once if using UCX
       IF (Input_Opt%LSETH2O) Input_Opt%LSETH2O = .FALSE.
    ENDIF

#if defined( MODEL_GEOS )
    ! GEOS-5 only: needed in GCHP?
    ! Compute the cosine of the solar zenith angle array:
    !    State_Met%SUNCOS     => COS(SZA) at the current time
    !    State_Met%SUNCOSmid  => COS(SZA) at the midpt of the chem timestep
    !    COS(SZA) at the midpt of the chem timestep 5hrs ago is now
    !    calculated elsewhere, in the HEMCO PARANOx extension
    CALL GET_COSINE_SZA( am_I_Root, Input_Opt, State_Grid, State_Met, RC )
    _ASSERT(RC==GC_SUCCESS, 'informative message here')
#endif

    !=======================================================================
    ! EMISSIONS. Pass HEMCO Phase 1 which only updates the HEMCO clock
    ! and the HEMCO data list. Should be called every time to make sure 
    ! that the HEMCO clock and the HEMCO data list are up to date.
    !=======================================================================
    HCO_PHASE = 1
    CALL EMISSIONS_RUN( am_I_Root,  Input_Opt, State_Chm, State_Diag, &
                        State_Grid, State_Met, DoEmis, HCO_PHASE, RC  )
    _ASSERT(RC==GC_SUCCESS, 'informative message here')

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!                                PHASE 1 or -1                           !!!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !=======================================================================
    ! 1. Convection
    ! 
    ! Call GEOS-Chem internal convection routines if convection is enabled
    ! in input.geos. This should only be done if convection is not covered
    ! by another gridded component and/or the GC species are not made
    ! friendly to this component!!
    !=======================================================================
    IF ( DoConv ) THEN
       if(am_I_Root.and.NCALLS<10) write(*,*) ' --- Do convection now'
       CALL MAPL_TimerOn( STATE, 'GC_CONV' )

       CALL DO_CONVECTION ( am_I_Root,  Input_Opt, State_Chm, State_Diag, &
                            State_Grid, State_Met, RC )
       _ASSERT(RC==GC_SUCCESS, 'informative message here')
 
       CALL MAPL_TimerOff( STATE, 'GC_CONV' )
       if(am_I_Root.and.NCALLS<10) write(*,*) ' --- Convection done!'
    ENDIF   

    !=======================================================================
    ! 2. Dry deposition
    !
    ! Calculates the deposition rates in [s-1].
    !=======================================================================
    IF ( DoDryDep ) THEN
       if(am_I_Root.and.NCALLS<10) THEN
          write(*,*) ' --- Do drydep now'
          write(*,*) '     Use FULL PBL: ', Input_Opt%PBL_DRYDEP
       endif
       CALL MAPL_TimerOn( STATE, 'GC_DRYDEP' )
    
       ! Do dry deposition
       CALL Do_DryDep ( am_I_Root,  Input_Opt, State_Chm, State_Diag, &
                        State_Grid, State_Met, RC ) 
       _ASSERT(RC==GC_SUCCESS, 'informative message here')

       CALL MAPL_TimerOff( STATE, 'GC_DRYDEP' )
       if(am_I_Root.and.NCALLS<10) write(*,*) ' --- Drydep done!'
    ENDIF

    !=======================================================================
    ! 3. Emissions (HEMCO)
    !
    ! HEMCO must be called on first time step to make sure that the HEMCO
    ! data lists are all properly set up. 
    !=======================================================================
    IF ( DoEmis ) THEN
       ! Optional memory prints (level >= 3)
       if ( MemDebugLevel > 0 ) THEN
          call ESMF_VMBarrier(VM, RC=STATUS)
          _VERIFY(STATUS)
          call MAPL_MemUtilsWrite(VM, &
                  'gigc_chunk_run, before Emissions_Run', RC=STATUS )
          _VERIFY(STATUS)
       endif

       if(am_I_Root.and.NCALLS<10) write(*,*) ' --- Do emissions now'
       CALL MAPL_TimerOn( STATE, 'GC_EMIS' )

       ! Do emissions. Pass HEMCO Phase 2 which performs the emissions 
       ! calculations.
       HCO_PHASE = 2
       CALL EMISSIONS_RUN ( am_I_Root,  Input_Opt, State_Chm, State_Diag, &
                            State_Grid, State_Met, DoEmis, HCO_PHASE, RC )
       _ASSERT(RC==GC_SUCCESS, 'informative message here')

       CALL MAPL_TimerOff( STATE, 'GC_EMIS' )
       if(am_I_Root.and.NCALLS<10) write(*,*) ' --- Emissions done!'

       ! Optional memory prints (level >= 3)
       if ( MemDebugLevel > 0 ) THEN
          call ESMF_VMBarrier(VM, RC=STATUS)
          _VERIFY(STATUS)
          call MAPL_MemUtilsWrite(VM,&
                  'gigc_chunk_run, after  Emissions_Run', RC=STATUS )
          _VERIFY(STATUS)
       endif

    ENDIF

    !=======================================================================
    ! If physics covers turbulence, simply add the emission and dry 
    ! deposition fluxes calculated above to the tracer array, without caring
    ! about the vertical distribution. The tracer tendencies are only added
    ! to the tracers array after emissions, drydep. So we need to use the
    ! emissions time step here.
    !=======================================================================
    IF ( DoTend ) THEN 
       if(am_I_Root.and.NCALLS<10) write(*,*)   &
                           ' --- Add emissions and drydep to tracers'
       CALL MAPL_TimerOn( STATE, 'GC_FLUXES' )

       ! Get emission time step [s]. 
       _ASSERT(ASSOCIATED(HcoState), 'informative message here')
       DT = HcoState%TS_EMIS 

       ! Apply tendencies over entire PBL. Use emission time step.
       CALL DO_TEND ( am_I_Root,  Input_Opt, State_Chm, State_Diag, &
                      State_Grid, State_Met, .FALSE., RC, DT=DT )
       _ASSERT(RC==GC_SUCCESS, 'informative message here')

       ! testing only
       if(am_I_Root.and.NCALLS<10) write(*,*)   &
                                 '     Tendency time step [s]: ', DT 

       CALL MAPL_TimerOff( STATE, 'GC_FLUXES' )
       if(am_I_Root.and.NCALLS<10) write(*,*)   &
                                 ' --- Fluxes applied to tracers!' 
    ENDIF ! Tendencies 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!                              PHASE 2 or -1                             !!!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !=======================================================================
    ! 4. Turbulence
    !
    ! Call GEOS-Chem internal turbulence routines if turbulence is enabled
    ! in input.geos. This should only be done if turbulence is not covered
    ! by another gridded component and/or the GC species are not made
    ! friendly to this component!!
    !=======================================================================
    IF ( DoTurb ) THEN
       if(am_I_Root.and.NCALLS<10) write(*,*) ' --- Do turbulence now'
       CALL MAPL_TimerOn( STATE, 'GC_TURB' )

       ! Do mixing and apply tendencies. This will use the dynamic time step,
       ! which is fine since this call will be executed on every time step. 
       CALL DO_MIXING ( am_I_Root,  Input_Opt, State_Chm, State_Diag, &
                        State_Grid, State_Met, RC )
       _ASSERT(RC==GC_SUCCESS, 'informative message here')

       CALL MAPL_TimerOff( STATE, 'GC_TURB' )
       if(am_I_Root.and.NCALLS<10) write(*,*) ' --- Turbulence done!'
    ENDIF

    ! Set tropospheric CH4 concentrations and fill species array with
    ! current values. 
#if defined( MODEL_GEOS )
!    IF ( DoTurb .OR. DoTend ) THEN
!      IF ( Input_Opt%LCH4SBC ) THEN
    IF ( Phase /= 2 ) THEN
#else
    IF ( Phase /= 2 .AND. Input_Opt%ITS_A_FULLCHEM_SIM  &
         .AND. IND_('CH4','A') > 0 ) THEN
#endif
       CALL SET_CH4 ( am_I_Root,  Input_Opt, State_Chm, State_Diag, &
                      State_Grid, State_Met, RC )
       _ASSERT(RC==GC_SUCCESS, 'informative message here')
    ENDIF

    !=======================================================================
    ! 5. Chemistry
    !=======================================================================
    IF ( DoChem ) THEN
       if(am_I_Root.and.NCALLS<10) write(*,*) ' --- Do chemistry now'
       CALL MAPL_TimerOn( STATE, 'GC_CHEM' )

       ! Calculate TOMS O3 overhead. For now, always use it from the
       ! Met field. State_Met%TO3 is imported from PCHEM (ckeller, 10/21/2014).
       CALL COMPUTE_OVERHEAD_O3( am_I_Root, State_Grid, DAY, &
                                 .TRUE., State_Met%TO3 )

#if !defined( MODEL_GEOS )
       ! Set H2O to species value if H2O is advected
       IF ( IND_('H2O','A') > 0 ) THEN
          CALL SET_H2O_TRAC( am_I_Root, (.not. Input_Opt%LUCX), Input_Opt, &
                             State_Chm, State_Grid, State_Met, RC )
       ENDIF
#endif

       ! Optional memory prints (level >= 3)
       if ( MemDebugLevel > 0 ) THEN
          call ESMF_VMBarrier(VM, RC=STATUS)
          _VERIFY(STATUS)
          call MAPL_MemUtilsWrite(VM, &
                  'gigc_chunk_run:, before Do_Chemistry', RC=STATUS )
          _VERIFY(STATUS)
       endif

       ! Do chemistry
       CALL Do_Chemistry( am_I_Root,  Input_Opt, State_Chm, State_Diag, &
                          State_Grid, State_Met, RC ) 
       _ASSERT(RC==GC_SUCCESS, 'informative message here')

       CALL MAPL_TimerOff( STATE, 'GC_CHEM' )
       if(am_I_Root.and.NCALLS<10) write(*,*) ' --- Chemistry done!'

       ! Optional memory prints (level >= 3)
       if ( MemDebugLevel > 0 ) THEN
          call ESMF_VMBarrier(VM, RC=STATUS)
          _VERIFY(STATUS)
          call MAPL_MemUtilsWrite(VM, &
                  'gigc_chunk_run, after  Do_Chemistry', RC=STATUS )
          _VERIFY(STATUS)
       endif

    ENDIF

    !=======================================================================
    ! 6. Wet deposition
    !=======================================================================
    IF ( DoWetDep ) THEN
       if(am_I_Root.and.NCALLS<10) write(*,*) ' --- Do wetdep now'
       CALL MAPL_TimerOn( STATE, 'GC_WETDEP' )

       ! Do wet deposition
       CALL DO_WETDEP( am_I_Root,  Input_Opt, State_Chm, State_Diag, &
                       State_Grid, State_Met, RC )
       _ASSERT(RC==GC_SUCCESS, 'informative message here')

       CALL MAPL_TimerOff( STATE, 'GC_WETDEP' )
       if(am_I_Root.and.NCALLS<10) write(*,*) ' --- Wetdep done!'
    ENDIF

    !=======================================================================
    ! Diagnostics 
    !=======================================================================

    !==============================================================
    !      ***** U P D A T E  O P T I C A L  D E P T H *****          
    !==============================================================
    ! Recalculate the optical depth at the wavelength(s) specified
    ! in the Radiation Menu. This must be done before the call to any
    ! diagnostic and only on a chemistry timestep.
    ! (skim, 02/05/11)
#if defined( MODEL_GEOS )
    ! RECOMPUTE_OD also contains a call to AEROSOL_CONC, which updates
    ! the PM25 diagnostics.
    !IF ( DoChem .AND. ND21 > 0 ) THEN
    !IF ( DoChem ) THEN
    ! Recompute_OD populates the AD21 diagnostics. Make sure to flush
    ! them first (ckeller, 8/9/17)
    ! Will need to come up with an alternative since all AD arrays
    ! will be removed with binary diagnostics (ewl, 12/1/18)
    IF ( Input_Opt%ND21 > 0 ) THEN
          AD21(:,:,:,:) = 0.0
    ENDIF
#else
    IF ( DoChem ) THEN
#endif
       CALL RECOMPUTE_OD ( am_I_Root,  Input_Opt, State_Chm, State_Diag, &
                           State_Grid, State_Met, RC )
       _ASSERT(RC==GC_SUCCESS, 'informative message here')
#if defined( MODEL_GEOS )
       !ENDIF
       !ENDIF
#else
    ENDIF
#endif

    if(am_I_Root.and.NCALLS<10) write(*,*) ' --- Do diagnostics now'
    CALL MAPL_TimerOn( STATE, 'GC_DIAGN' )

    ! Set certain diagnostics dependent on state at end of step. This
    ! includes species concentration and dry deposition flux.
    CALL Set_Diagnostics_EndofTimestep( am_I_Root,  Input_Opt,  &
                                        State_Chm,  State_Diag, &
                                        State_Grid, State_Met, RC )
    _ASSERT(RC==GC_SUCCESS, 'informative message here')

    ! Archive aerosol mass and PM2.5 diagnostics
    IF ( State_Diag%Archive_AerMass ) THEN
       CALL Set_AerMass_Diagnostic( am_I_Root,  Input_Opt,  State_Chm, &
                                    State_Diag, State_Grid, State_Met, RC )
       _ASSERT(RC==GC_SUCCESS, 'informative message here')
    ENDIF

    CALL MAPL_TimerOff( STATE, 'GC_DIAGN' )
    if(am_I_Root.and.NCALLS<10) write(*,*) ' --- Diagnostics done!'

    !=======================================================================
    ! Convert State_Chm%Species units
    !=======================================================================
    CALL Convert_Spc_Units ( am_I_Root, Input_Opt, State_Chm, State_Grid, &
                             State_Met, OrigUnit, RC )

#if defined( MODEL_GEOS )
    ! Save specific humidity and dry air mass for total mixing ratio 
    ! adjustment in next timestep, if needed (ewl, 11/8/18)
    State_Met%SPHU_PREV = State_Met%SPHU
#endif

#ifdef ADJOINT
       if (Input_Opt%IS_FD_SPOT_THIS_PET .and. Input_opt%IFD > 0) THEN
       DO N = 1, State_Chm%nSpecies
          ThisSpc => State_Chm%SpcData(N)%Info
          write(*,*) 'SpcAdj(', TRIM(thisSpc%Name), ') = ',  &
               State_Chm%SpeciesAdj(Input_Opt%IFD,Input_Opt%JFD,Input_Opt%LFD,N)
       ENDDO
       ENDIF
    !=======================================================================
    ! If this is an adjoint run, we need to check for the final (first)
    ! timestep and multiply the scaling factor adjoint by the initial concs
    !=======================================================================
    IF (Input_Opt%IS_ADJOINT .and. IsStarttime) THEN
       if (am_I_Root) WRITE(*,*) '   Adjoint multiplying SF_ADJ by ICS'
       DO N = 1, State_Chm%nSpecies
          ThisSpc => State_Chm%SpcData(N)%Info

          ! Find the non-adjoint variable or this
          TRACNAME = ThisSpc%Name

          State_Chm%SpeciesAdj(:,:,:,N) = State_Chm%SpeciesAdj(:,:,:,N) * State_Chm%Species(:,:,:,N)
          if (Input_Opt%IS_FD_SPOT_THIS_PET .and. Input_Opt%IFD > 0) THEN
             write(*,*) 'After conversion ',  &
                  State_Chm%SpeciesAdj(Input_Opt%IFD,Input_Opt%JFD,Input_Opt%LFD,N)
          ENDIF
       ENDDO

       CALL Set_SpcAdj_Diagnostic( am_I_Root, 'SpeciesAdj',                 &
                                   State_Diag%SpeciesAdj,                   &
                                   Input_Opt,  State_Chm,                   &
                                   State_Grid, State_Met,  RC              )
    ENDIF
#endif


    !=======================================================================
    ! Clean up
    !=======================================================================

    ! testing only
    IF ( PHASE /= 1 .AND. NCALLS < 10 ) NCALLS = NCALLS + 1 

    ! First call is done
    FIRST = .FALSE.

    ! Return success
    RC = GC_SUCCESS

  END SUBROUTINE GIGC_Chunk_Run
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gigc_chunk_final
!
! !DESCRIPTION: Subroutine GIGC\_CHUNK\_FINAL is the ESMF finalize method for
!  GEOS-Chem.  This routine deallocates pointers and arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GIGC_Chunk_Final( am_I_Root,  Input_Opt, State_Chm, State_Diag, &
                               State_Grid, State_Met, RC )
!
! !USES:
!
    USE Input_Opt_Mod,         ONLY : OptInput
#if !defined( MODEL_GEOS )
    USE Input_Opt_Mod,         ONLY : Cleanup_Input_Opt
#endif
    USE State_Chm_Mod,         ONLY : ChmState, Cleanup_State_Chm
    USE State_Diag_Mod,        ONLY : DgnState, Cleanup_State_Diag
    USE State_Grid_Mod,        ONLY : GrdState, Cleanup_State_Grid
    USE State_Met_Mod,         ONLY : MetState, Cleanup_State_Met
    USE HCOI_GC_MAIN_MOD,      ONLY : HCOI_GC_FINAL
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root     ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt     ! Input Options object
    TYPE(ChmState), INTENT(INOUT) :: State_Chm     ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag    ! Diagnostics State object
    TYPE(GrdState), INTENT(INOUT) :: State_Grid    ! Grid State object
    TYPE(MetState), INTENT(INOUT) :: State_Met     ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC            ! Success or failure
!
! !REVISION HISTORY: 
!  18 Jul 2011 - M. Long     - Initial Version
!  09 Oct 2012 - R. Yantosca - Added comments & cosmetic changes
!  16 Oct 2012 - R. Yantosca - Renamed GC_STATE argument to State_Chm
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_chm_mod.F90
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_met_mod.F90
!  22 Oct 2012 - R. Yantosca - Renamed to GIGC_Chunk_Final
!  01 Nov 2012 - R. Yantosca - Now reference gigc_input_opt_mod.F90
!  19 Sep 2017 - E. Lundgren - Move gigc_finalize content to within this routine
!  14 Dec 2018 - E. Lundgren - Combine offline and GEOS-5 code and simplify
!EOP
!------------------------------------------------------------------------------
!BOC
! LOCAL VARIABLES:
!
#ifdef ADJOINT
    ! Finite difference test variables
    INTEGER                        :: IFD, JFD, LFD, NFD
    INTEGER                        :: I, J, L
    REAL*8                         :: CFN
    CHARACTER(len=ESMF_MAXSTR)     :: FD_SPEC
#endif

    ! Assume succes
    RC = GC_SUCCESS

#ifdef ADJOINT
    IF (Input_Opt%IS_FD_SPOT_THIS_PET .and. .not. Input_Opt%IS_FD_GLOBAL) THEN
       FD_SPEC = transfer(state_chm%SpcData(Input_Opt%NFD)%Info%Name, FD_SPEC)
       IFD = Input_Opt%IFD
       JFD = Input_Opt%JFD
       LFD = Input_Opt%LFD
       NFD = Input_Opt%NFD
       ! print out the cost function
       WRITE(*,*) ' Computing final cost function'
       CFN = 0d0
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX
          if (State_Chm%CostFuncMask(I,J,L) > 0d0) THEN
             WRITE (*, 1047) I, J, L, state_chm%species(I, J, L, NFD)
             CFN = CFN + state_chm%Species(I, J, L, NFD)
          endif
       ENDDO
        ENDDO
       ENDDO
       WRITE(*,'(a7, e22.10)') ' CFN = ', CFN
1047   FORMAT('  SPC(', i2, ', ', i2, ', ', i2, ') = ', e22.10)
    ENDIF
#endif

    ! Finalize HEMCO
    CALL HCOI_GC_FINAL( am_I_Root, .FALSE., RC )
    IF ( am_I_Root ) THEN
       IF ( RC == GC_SUCCESS ) THEN
          write(*,'(a)') 'HEMCO::Finalize... OK.'
       ELSE
          write(*,'(a)') 'HEMCO::Finalize... FAILURE.'
       ENDIF
    ENDIF

    ! Deallocate fields of the Chemistry State object
    CALL Cleanup_State_Chm( am_I_Root, State_Chm, RC )
    IF ( am_I_Root ) THEN
       IF ( RC == GC_SUCCESS ) THEN
          write(*,'(a)') 'Chem::State_Chm Finalize... OK.'
       ELSE
          write(*,'(a)') 'Chem::State_Chm Finalize... FAILURE.'
       ENDIF
    ENDIF

    ! Deallocate fields of the Diagnostics State object
    CALL Cleanup_State_Diag( am_I_Root, State_Diag, RC )
    IF ( am_I_Root ) THEN
       IF ( RC == GC_SUCCESS ) THEN
          write(*,'(a)') 'Chem::State_Diag Finalize... OK.'
       ELSE
          write(*,'(a)') 'Chem::State_Diag Finalize... FAILURE.'
       ENDIF
    ENDIF

    ! Deallocate fields of the Grid State object
    CALL Cleanup_State_Grid( am_I_Root, State_Grid, RC )
    IF ( am_I_Root ) THEN
       IF ( RC == GC_SUCCESS ) THEN
          write(*,'(a)') 'Chem::State_Grid Finalize... OK.'
       ELSE
          write(*,'(a)') 'Chem::State_Grid Finalize... FAILURE.'
       ENDIF
    ENDIF

    ! Deallocate fields of the Meteorology State object
    CALL Cleanup_State_Met( am_I_Root, State_Met, RC )
    IF ( am_I_Root ) THEN
       IF ( RC == GC_SUCCESS ) THEN
          write(*,'(a)') 'Chem::State_Met Finalize... OK.'
       ELSE
          write(*,'(a)') 'Chem::State_Met Finalize... FAILURE.'
       ENDIF
    ENDIF

#if !defined( MODEL_GEOS )
    ! Deallocate fields of the Input Options object
    ! The call to Cleanup_Input_Opt causes a memory leak error. Comment
    ! for now (ckeller, 11/29/16).
    ! Does this still cause a memory leak? (ewl, 12/14/18)
     CALL Cleanup_Input_Opt( am_I_Root, Input_Opt, RC )
    IF ( am_I_Root ) THEN
       IF ( RC == GC_SUCCESS ) THEN
          write(*,'(a)') 'Chem::Input_Opt Finalize... OK.'
       ELSE
          write(*,'(a)') 'Chem::Input_Opt Finalize... FAILURE.'
       ENDIF
    ENDIF
#endif

  END SUBROUTINE GIGC_Chunk_Final
!EOC
#if defined( MODEL_GEOS )
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_ozonopause
!
! !DESCRIPTION: 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SET_OZONOPAUSE ( am_I_Root,  Input_Opt, State_Chm, &
                              State_Grid, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState, Ind_
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
    USE PRESSURE_MOD,   ONLY : GET_PCENTER
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root     ! Are we on the root CPU?
    TYPE(OptInput), INTENT(IN)    :: Input_Opt     ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid    ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm     ! Chem state object
    TYPE(MetState), INTENT(INOUT) :: State_Met     ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC            ! Success or failure
!
! !REMARKS:
! 
! !REVISION HISTORY: 
!  10 Nov 2015 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: I, J, L, N, IDTO3
    LOGICAL            :: IsInStrat
    REAL(fp)           :: O3val
    INTEGER, PARAMETER :: NLEVEL = 3

    ! Assume success
    RC = GC_SUCCESS

    ! Species index
    IDTO3 = Ind_('O3')

    ! Return here if ozonopause parameter has invalid value
    IF ( Input_Opt%OZONOPAUSE <= 0.0_fp .OR. IDTO3 < 0 ) RETURN

    ! Reset tropopause pressures
    State_Met%TROPP(:,:) = 0.0_fp

    ! ozonopause value (ppb --> v/v) 
    O3val = Input_Opt%OZONOPAUSE * 1.0e-9_fp

    ! Loop over grid boxes on this PET
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       IsInStrat = .FALSE.

       ! Find first level where ozone concentration exceeds threshold 
       DO L = 1, State_Grid%NZ
          IF ( State_Chm%Species(I,J,L,IDTO3) >= O3val ) THEN
             ! Sanity check: level above should have higher ozone and  
             ! pressure should be above 500hPa
             IF ( L > ( LM-NLEVEL+1 ) ) THEN
                IsInStrat = .TRUE.
             ELSEIF ( GET_PCENTER(I,J,L) >= 500.0_fp ) THEN
                IsInStrat = .FALSE.
             ELSE
                IsInStrat = .TRUE.
                DO N = 1, NLEVEL
                   IF ( State_Chm%Species(I,J,L+N,IDTO3) < State_Chm%Species(I,J,L,IDTO3) ) THEN
                      IsInStrat = .FALSE.
                      EXIT
                   ENDIF
                ENDDO
             ENDIF
          ENDIF

          ! Set TROPP to value in middle of this grid box
          IF ( IsInStrat ) THEN
                State_Met%TROPP(I,J) = GET_PCENTER(I,J,L) 
             EXIT ! End L loop
          ENDIF
       ENDDO
    ENDDO
    ENDDO

  END SUBROUTINE SET_OZONOPAUSE
!EOC
#endif

!BOP
  SUBROUTINE GIGC_PRINT_MET(am_I_Root, I, J, L,         &
                            Input_Opt, State_Grid, State_Met, LOC, RC )

!
! !USES:
!
    USE State_Met_Mod,        ONLY : MetState
    USE Input_Opt_Mod,        ONLY : OptInput
    USE State_Grid_Mod,       ONLY : GrdState

!
! !INPUT PARAMETERS: 
!
      LOGICAL,          INTENT(IN)    :: am_I_Root ! Are we on root CPU?
      INTEGER,          INTENT(IN)    :: I         ! Grid cell lat index
      INTEGER,          INTENT(IN)    :: J         ! Grid cell lon index
      INTEGER,          INTENT(IN)    :: L         ! Grid cell lev index
      CHARACTER(LEN=*), INTENT(IN)    :: LOC       ! Call location string
      TYPE(OptInput),   INTENT(IN)    :: Input_Opt ! Input Options object
      TYPE(GrdState),   INTENT(IN)    :: State_Grid! Grid State object
      TYPE(MetState),   INTENT(IN)    :: State_Met ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS: 
!

!
! !OUTPUT PARAMETERS:
!
      INTEGER,          INTENT(OUT)   :: RC        ! Success or failure?! 
! !REMARKS:
!
! !REVISION HISTORY: 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!     
      CHARACTER(LEN=255) :: ErrorMsg, ThisLoc


      !=========================================================================
      ! GIGC_PRINT_MET begins here!
      !=========================================================================

      ErrorMsg  = ''
      ThisLoc   = ' -> at Print_Global_Species_Kg (in module ' // &
                  'GeosUtil/unitconv_mod.F)'

      ! Echo info
      IF ( am_I_Root ) THEN
         WRITE( 6, 100 ) TRIM( LOC )
         WRITE( 6, 113 ) State_Grid%YMid(I,J), State_Grid%XMid(I,J)
      ENDIF
100   FORMAT( /, '%%%%% GIGC_PRINT_MET at ', a )
113   FORMAT( 'Lat: ', f5.1, '   Lon: ', f5.1 )

      ! Write formatted output
      IF ( am_I_Root ) THEN
         ! 2-D Fields
         WRITE( 6, 114 ) 'PBLH',     State_Met%PBLH(I,J),     I, J
         WRITE( 6, 114 ) 'PSC2_WET', State_Met%PSC2_WET(I,J), I, J
         WRITE( 6, 114 ) 'PSC2_DRY', State_Met%PSC2_DRY(I,J), I, J
         WRITE( 6, 114 ) 'PS1_WET',  State_Met%PS1_WET(I,J), I, J
         WRITE( 6, 114 ) 'PS1_DRY',  State_Met%PS1_DRY(I,J), I, J
         WRITE( 6, 114 ) 'PS2_WET',  State_Met%PS2_WET(I,J), I, J
         WRITE( 6, 114 ) 'PS2_DRY',  State_Met%PS2_DRY(I,J), I, J
         WRITE( 6, 114 ) 'TS',       State_Met%TS(I,J),       I, J
         WRITE( 6, 114 ) 'U10M',     State_Met%U10M(I,J),     I, J
         ! 3-D Fields
         WRITE( 6, 115 ) 'CLDF',     State_Met%CLDF(I,J,L),      I, J, L
         WRITE( 6, 115 ) 'OMEGA',    State_Met%OMEGA(I,J,L),     I, J, L
         WRITE( 6, 115 ) 'PEDGE',    State_Met%PEDGE(I,J,L),     I, J, L
         WRITE( 6, 115 ) 'T',        State_Met%T(I,J,L),         I, J, L
         WRITE( 6, 115 ) 'U',        State_Met%U(I,J,L),         I, J, L
         WRITE( 6, 115 ) 'V',        State_Met%V(I,J,L),         I, J, L
         WRITE( 6, 115 ) 'AD',       State_Met%AD(I,J,L),        I, J, L
         WRITE( 6, 115 ) 'PREVSPHU', State_Met%SPHU_PREV(I,J,L), I, J, L
         WRITE( 6, 115 ) 'SPHU',     State_Met%SPHU(I,J,L),      I, J, L
         ! terminator
         WRITE( 6, 120 )
      ENDIF
 114  FORMAT( 'Grid cell  for ', a8, ' = ', es24.16, ', I,J  = ',2I4 )
 115  FORMAT( 'Grid cell  for ', a8, ' = ', es24.16, ', I,J,L= ',3I4 )
 120  FORMAT( / )


  END SUBROUTINE GIGC_PRINT_MET

END MODULE GIGC_Chunk_Mod
