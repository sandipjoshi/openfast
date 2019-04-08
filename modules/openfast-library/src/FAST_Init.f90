!**********************************************************************************************************************************
! FAST_Solver.f90, FAST_Subs.f90, FAST_Lin.f90, and FAST_Mods.f90 make up the FAST glue code in the FAST Modularization Framework.
! FAST_Prog.f90, FAST_Library.f90, FAST_Prog.c are different drivers for this code.
!..................................................................................................................................
! LICENSING
! Copyright (C) 2013-2016  National Renewable Energy Laboratory
!
!    This file is part of FAST.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!**********************************************************************************************************************************
MODULE FAST_Initialization

   USE FAST_ModTypes
   USE FAST_IO
   USE FAST_VTK, ONLY: SetVTKParameters_B4HD, SetVTKParameters
   USE FAST_Linear
   USE FAST_Solver,       ONLY: initmodulemappings
   USE AeroDyn,           ONLY: AD_Init
   USE AeroDyn14,         ONLY: AD14_Init
   USE InflowWind,        ONLY: InflowWind_Init
   USE ElastoDyn,         ONLY: ED_Init
   USE BeamDyn,           ONLY: BD_Init
   USE FEAMooring,        ONLY: FEAM_Init
   USE MoorDyn,           ONLY: MD_Init
   USE MAP,               ONLY: MAP_Init
   USE OrcaFlexInterface, ONLY: Orca_Init
   USE HydroDyn,          ONLY: HydroDyn_Init
   USE IceDyn,            ONLY: IceD_Init
   USE IceFloe,           ONLY: IceFloe_Init
   USE ServoDyn,          ONLY: SrvD_Init, Cmpl4SFun, Cmpl4LV
   USE SubDyn,            ONLY: SD_Init
   USE OpenFOAM,          ONLY: Init_OpFM
   USE SuperController,   ONLY: Init_SC
   Use ExtPtfm_MCKF,      ONLY: ExtPtfm_Init
     
   IMPLICIT NONE

CONTAINS
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! INITIALIZATION ROUTINES
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> a wrapper routine to call FAST_Initialize a the full-turbine simulation level (makes easier to write top-level driver)
SUBROUTINE FAST_InitializeAll_T( t_initial, TurbID, Turbine, ErrStat, ErrMsg, InFile, ExternInitData )

   REAL(DbKi),                        INTENT(IN   ) :: t_initial      !< initial time
   INTEGER(IntKi),                    INTENT(IN   ) :: TurbID         !< turbine Identifier (1-NumTurbines)
   TYPE(FAST_TurbineType),            INTENT(INOUT) :: Turbine        !< all data for one instance of a turbine
   INTEGER(IntKi),                    INTENT(  OUT) :: ErrStat        !< Error status of the operation
   CHARACTER(*),                      INTENT(  OUT) :: ErrMsg         !< Error message if ErrStat /= ErrID_None
   CHARACTER(*),             OPTIONAL,INTENT(IN   ) :: InFile         !< A CHARACTER string containing the name of the primary FAST input file (if not present, we'll get it from the command line)   
   TYPE(FAST_ExternInitType),OPTIONAL,INTENT(IN   ) :: ExternInitData !< Initialization input data from an external source (Simulink)
   
   Turbine%TurbID = TurbID  
   
   
   IF (PRESENT(InFile)) THEN
      IF (PRESENT(ExternInitData)) THEN
         CALL FAST_InitializeAll( t_initial, Turbine%p_FAST, Turbine%y_FAST, Turbine%m_FAST, &
                     Turbine%ED, Turbine%BD, Turbine%SrvD, Turbine%AD14, Turbine%AD, Turbine%IfW, Turbine%OpFM, Turbine%SC,&
                     Turbine%HD, Turbine%SD, Turbine%ExtPtfm, Turbine%MAP, Turbine%FEAM, Turbine%MD, Turbine%Orca, &
                     Turbine%IceF, Turbine%IceD, Turbine%MeshMapData, ErrStat, ErrMsg, InFile, ExternInitData )
      ELSE         
         CALL FAST_InitializeAll( t_initial, Turbine%p_FAST, Turbine%y_FAST, Turbine%m_FAST, &
                     Turbine%ED, Turbine%BD, Turbine%SrvD, Turbine%AD14, Turbine%AD, Turbine%IfW, Turbine%OpFM, Turbine%SC, &
                     Turbine%HD, Turbine%SD, Turbine%ExtPtfm, Turbine%MAP, Turbine%FEAM, Turbine%MD, Turbine%Orca, &
                     Turbine%IceF, Turbine%IceD, Turbine%MeshMapData, ErrStat, ErrMsg, InFile  )
      END IF
   ELSE
      CALL FAST_InitializeAll( t_initial, Turbine%p_FAST, Turbine%y_FAST, Turbine%m_FAST, &
                     Turbine%ED, Turbine%BD, Turbine%SrvD, Turbine%AD14, Turbine%AD, Turbine%IfW, Turbine%OpFM, Turbine%SC, &
                     Turbine%HD, Turbine%SD, Turbine%ExtPtfm, Turbine%MAP, Turbine%FEAM, Turbine%MD, Turbine%Orca, &
                     Turbine%IceF, Turbine%IceD, Turbine%MeshMapData, ErrStat, ErrMsg )
   END IF
   
         
END SUBROUTINE FAST_InitializeAll_T
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to call Init routine for each module. This routine sets all of the init input data for each module.
SUBROUTINE FAST_InitializeAll( t_initial, p_FAST, y_FAST, m_FAST, ED, BD, SrvD, AD14, AD, IfW, OpFM, SC, HD, SD, ExtPtfm, &
                               MAPp, FEAM, MD, Orca, IceF, IceD, MeshMapData, ErrStat, ErrMsg, InFile, ExternInitData )

   use ElastoDyn_Parameters, only: Method_RK4

   REAL(DbKi),               INTENT(IN   ) :: t_initial           !< initial time
   TYPE(FAST_ParameterType), INTENT(INOUT) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code
   TYPE(FAST_MiscVarType),   INTENT(INOUT) :: m_FAST              !< Miscellaneous variables
     
   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  !< BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(INOUT) :: SrvD                !< ServoDyn data
   TYPE(AeroDyn14_Data),     INTENT(INOUT) :: AD14                !< AeroDyn14 data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  !< AeroDyn data
   TYPE(InflowWind_Data),    INTENT(INOUT) :: IfW                 !< InflowWind data
   TYPE(OpenFOAM_Data),      INTENT(INOUT) :: OpFM                !< OpenFOAM data
   TYPE(SuperController_Data), INTENT(INOUT) :: SC                !< SuperController data
   TYPE(HydroDyn_Data),      INTENT(INOUT) :: HD                  !< HydroDyn data
   TYPE(SubDyn_Data),        INTENT(INOUT) :: SD                  !< SubDyn data
   TYPE(ExtPtfm_Data),       INTENT(INOUT) :: ExtPtfm             !< ExtPtfm_MCKF data
   TYPE(MAP_Data),           INTENT(INOUT) :: MAPp                !< MAP data
   TYPE(FEAMooring_Data),    INTENT(INOUT) :: FEAM                !< FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(INOUT) :: MD                  !< Data for the MoorDyn module
   TYPE(OrcaFlex_Data),      INTENT(INOUT) :: Orca                !< OrcaFlex interface data
   
   TYPE(IceFloe_Data),       INTENT(INOUT) :: IceF                !< IceFloe data
   TYPE(IceDyn_Data),        INTENT(INOUT) :: IceD                !< All the IceDyn data used in time-step loop

   TYPE(FAST_ModuleMapType), INTENT(INOUT) :: MeshMapData         !< Data for mapping between modules
      
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None
   CHARACTER(*), OPTIONAL,   INTENT(IN   ) :: InFile              !< A CHARACTER string containing the name of the primary FAST input file (if not present, we'll get it from the command line)
   
   TYPE(FAST_ExternInitType), OPTIONAL, INTENT(IN) :: ExternInitData !< Initialization input data from an external source (Simulink)
   
   ! local variables      
   CHARACTER(1024)                         :: InputFile           !< A CHARACTER string containing the name of the primary FAST input file
   TYPE(FAST_InitData)                     :: Init                !< Initialization data for all modules

       
   REAL(ReKi)                              :: AirDens             ! air density for initialization/normalization of OpenFOAM data
   REAL(DbKi)                              :: dt_IceD             ! tmp dt variable to ensure IceDyn doesn't specify different dt values for different legs (IceDyn instances)
   REAL(DbKi)                              :: dt_BD               ! tmp dt variable to ensure BeamDyn doesn't specify different dt values for different instances
   INTEGER(IntKi)                          :: ErrStat2
   INTEGER(IntKi)                          :: IceDim              ! dimension we're pre-allocating for number of IceDyn legs/instances
   INTEGER(IntKi)                          :: I                   ! generic loop counter
   INTEGER(IntKi)                          :: k                   ! blade loop counter
   logical                                 :: CallStart
   
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
                                           
   CHARACTER(*), PARAMETER                 :: RoutineName = 'FAST_InitializeAll'       
   
   
   !..........
   ErrStat = ErrID_None
   ErrMsg  = ""
      
   y_FAST%UnSum = -1                                                    ! set the summary file unit to -1 to indicate it's not open
   y_FAST%UnOu  = -1                                                    ! set the text output file unit to -1 to indicate it's not open
   y_FAST%UnGra = -1                                                    ! set the binary graphics output file unit to -1 to indicate it's not open
   
   p_FAST%WrVTK = VTK_Unknown                                           ! set this so that we can potentially output VTK information on initialization error
   p_FAST%VTK_tWidth = 1                                                ! initialize in case of error before reading the full file
   p_FAST%n_VTKTime  = 1                                                ! initialize in case of error before reading the full file
   y_FAST%VTK_LastWaveIndx = 1                                          ! Start looking for wave data at the first index
   y_FAST%VTK_count = 0                                                 ! first VTK file has 0 as output      
   y_FAST%n_Out = 0                                                     ! set the number of ouptut channels to 0 to indicate there's nothing to write to the binary file
   p_FAST%ModuleInitialized = .FALSE.                                   ! (array initialization) no modules are initialized 
   
      ! Get the current time
   CALL DATE_AND_TIME ( Values=m_FAST%StrtTime )                        ! Let's time the whole simulation
   CALL CPU_TIME ( m_FAST%UsrTime1 )                                    ! Initial time (this zeros the start time when used as a MATLAB function)
   m_FAST%UsrTime1 = MAX( 0.0_ReKi, m_FAST%UsrTime1 )                   ! CPU_TIME: If a meaningful time cannot be returned, a processor-dependent negative value is returned
   

   m_FAST%t_global        = t_initial - 20.                             ! initialize this to a number < t_initial for error message in ProgAbort
   m_FAST%calcJacobian    = .TRUE.                                      ! we need to calculate the Jacobian
   m_FAST%NextJacCalcTime = m_FAST%t_global                             ! We want to calculate the Jacobian on the first step
   p_FAST%TDesc           = ''
!   p_FAST%CheckHSSBrTrqC = .false.

   y_FAST%Lin%WindSpeed = 0.0_ReKi
   
   if (present(ExternInitData)) then
      CallStart = .not. ExternInitData%FarmIntegration ! .and. ExternInitData%TurbineID == 1
      if (ExternInitData%TurbineID > 0) p_FAST%TDesc = 'T'//trim(num2lstr(ExternInitData%TurbineID)) 
   else
      CallStart = .true.
   end if
           
   
      ! Init NWTC_Library, display copyright and version information:
   if (CallStart) then
      AbortErrLev = ErrID_Fatal                                 ! Until we read otherwise from the FAST input file, we abort only on FATAL errors
      CALL FAST_ProgStart( FAST_Ver )
      p_FAST%WrSttsTime = .TRUE.
   else
      ! if we don't call the start data (e.g., from FAST.Farm), we won't override AbortErrLev either 
      CALL DispNVD( FAST_Ver )
      p_FAST%WrSttsTime = .FALSE.
   end if
   
   IF (PRESENT(InFile)) THEN
      p_FAST%UseDWM = .FALSE.
      InputFile = InFile
   ELSE
      CALL GetInputFileName(InputFile,p_FAST%UseDWM,ErrStat2,ErrMsg2)            
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF (ErrStat >= AbortErrLev) THEN
            CALL Cleanup()
            RETURN
         END IF
   END IF
   
   ! ... Open and read input files ...
   ! also, set turbine reference position for graphics output
   if (PRESENT(ExternInitData)) then
      p_FAST%TurbinePos = ExternInitData%TurbinePos
      
      if (ExternInitData%FarmIntegration) then ! we're integrating with FAST.Farm
         CALL FAST_Init( p_FAST, m_FAST, y_FAST, t_initial, InputFile, ErrStat2, ErrMsg2, ExternInitData%TMax, OverrideAbortLev=.false., RootName=ExternInitData%RootName )         
      else
         CALL FAST_Init( p_FAST, m_FAST, y_FAST, t_initial, InputFile, ErrStat2, ErrMsg2, ExternInitData%TMax, ExternInitData%TurbineID )  ! We have the name of the input file and the simulation length from somewhere else (e.g. Simulink)         
      end if
      
   else
      p_FAST%TurbinePos = 0.0_ReKi
      CALL FAST_Init( p_FAST, m_FAST, y_FAST, t_initial, InputFile, ErrStat2, ErrMsg2 )                       ! We have the name of the input file from somewhere else (e.g. Simulink)
   end if
         
   CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF
      
      
   !...............................................................................................................................  
      
   p_FAST%dt_module = p_FAST%dt ! initialize time steps for each module   

   ! ........................
   ! initialize ElastoDyn (must be done first)
   ! ........................
   
   ALLOCATE( ED%Input( p_FAST%InterpOrder+1 ), ED%InputTimes( p_FAST%InterpOrder+1 ),STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating ED%Input and ED%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF
      
   Init%InData_ED%Linearize = p_FAST%Linearize
   Init%InData_ED%InputFile = p_FAST%EDFile
   IF ( p_FAST%CompAero == Module_AD14 ) THEN
      Init%InData_ED%ADInputFile = p_FAST%AeroFile
   ELSE
      Init%InData_ED%ADInputFile = ""
   END IF
   
   Init%InData_ED%RootName      = TRIM(p_FAST%OutFileRoot)//'.'//TRIM(y_FAST%Module_Abrev(Module_ED))
   Init%InData_ED%CompElast     = p_FAST%CompElast == Module_ED

   CALL ED_Init( Init%InData_ED, ED%Input(1), ED%p, ED%x(STATE_CURR), ED%xd(STATE_CURR), ED%z(STATE_CURR), ED%OtherSt(STATE_CURR), &
                  ED%y, ED%m, p_FAST%dt_module( MODULE_ED ), Init%OutData_ED, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
   p_FAST%ModuleInitialized(Module_ED) = .TRUE.
   CALL SetModuleSubstepTime(Module_ED, p_FAST, y_FAST, ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
      ! bjj: added this check per jmj; perhaps it would be better in ElastoDyn, but I'll leave it here for now:
   IF ( p_FAST%TurbineType == Type_Offshore_Floating ) THEN
      IF ( ED%p%TowerBsHt < 0.0_ReKi .AND. .NOT. EqualRealNos( ED%p%TowerBsHt, 0.0_ReKi ) ) THEN
         CALL SetErrStat(ErrID_Fatal,"ElastoDyn TowerBsHt must not be negative for floating offshore systems.",ErrStat,ErrMsg,RoutineName)
      END IF      
   END IF   

   allocate( y_FAST%Lin%Modules(MODULE_ED)%Instance(1), stat=ErrStat2)
   if (ErrStat2 /= 0 ) then
      call SetErrStat(ErrID_Fatal, "Error allocating Lin%Modules(ED).", ErrStat, ErrMsg, RoutineName )
   else
   
      if (allocated(Init%OutData_ED%LinNames_y)) call move_alloc(Init%OutData_ED%LinNames_y,y_FAST%Lin%Modules(MODULE_ED)%Instance(1)%Names_y)
      if (allocated(Init%OutData_ED%LinNames_x)) call move_alloc(Init%OutData_ED%LinNames_x,y_FAST%Lin%Modules(MODULE_ED)%Instance(1)%Names_x)
      if (allocated(Init%OutData_ED%LinNames_u)) call move_alloc(Init%OutData_ED%LinNames_u,y_FAST%Lin%Modules(MODULE_ED)%Instance(1)%Names_u)
      if (allocated(Init%OutData_ED%RotFrame_y)) call move_alloc(Init%OutData_ED%RotFrame_y,y_FAST%Lin%Modules(MODULE_ED)%Instance(1)%RotFrame_y)
      if (allocated(Init%OutData_ED%RotFrame_x)) call move_alloc(Init%OutData_ED%RotFrame_x,y_FAST%Lin%Modules(MODULE_ED)%Instance(1)%RotFrame_x)
      if (allocated(Init%OutData_ED%DerivOrder_x)) call move_alloc(Init%OutData_ED%DerivOrder_x,y_FAST%Lin%Modules(MODULE_ED)%Instance(1)%DerivOrder_x)
      if (allocated(Init%OutData_ED%RotFrame_u)) call move_alloc(Init%OutData_ED%RotFrame_u,y_FAST%Lin%Modules(MODULE_ED)%Instance(1)%RotFrame_u)
      if (allocated(Init%OutData_ED%IsLoad_u  )) call move_alloc(Init%OutData_ED%IsLoad_u  ,y_FAST%Lin%Modules(MODULE_ED)%Instance(1)%IsLoad_u  )
         
      if (allocated(Init%OutData_ED%WriteOutputHdr)) y_FAST%Lin%Modules(MODULE_ED)%Instance(1)%NumOutputs = size(Init%OutData_ED%WriteOutputHdr)
   end if

   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF
   
   if (p_FAST%CalcSteady) then
      if ( EqualRealNos(Init%OutData_ED%RotSpeed, 0.0_ReKi) ) then
         p_FAST%TrimCase = TrimCase_none
         p_FAST%NLinTimes = 1
         p_FAST%LinInterpOrder = 0 ! constant values
      elseif ( Init%OutData_ED%isFixed_GenDOF ) then
         p_FAST%TrimCase = TrimCase_none
      end if
   end if
   
   
   ! ........................
   ! initialize BeamDyn 
   ! ........................
   IF ( p_FAST%CompElast == Module_BD ) THEN      
      p_FAST%nBeams = Init%OutData_ED%NumBl          ! initialize number of BeamDyn instances = number of blades      
   ELSE
      p_FAST%nBeams = 0
   END IF

   ALLOCATE( BD%Input( p_FAST%InterpOrder+1, p_FAST%nBeams ), BD%InputTimes( p_FAST%InterpOrder+1, p_FAST%nBeams ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating BD%Input and BD%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF  
                        
   ALLOCATE( BD%x(           p_FAST%nBeams,2), &
             BD%xd(          p_FAST%nBeams,2), &
             BD%z(           p_FAST%nBeams,2), &
             BD%OtherSt(     p_FAST%nBeams,2), &
             BD%p(           p_FAST%nBeams  ), &
             BD%u(           p_FAST%nBeams  ), &
             BD%y(           p_FAST%nBeams  ), &
             BD%m(           p_FAST%nBeams  ), &
             Init%OutData_BD(p_FAST%nBeams  ), &
                                             STAT = ErrStat2 )                                                  
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating BeamDyn state, input, and output data.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF        
   
   IF (p_FAST%CompElast == Module_BD) THEN

      Init%InData_BD%DynamicSolve = .TRUE.       ! FAST can only couple to BeamDyn when dynamic solve is used.

      Init%InData_BD%Linearize = p_FAST%Linearize
      Init%InData_BD%gravity      = (/ 0.0_ReKi, 0.0_ReKi, -Init%OutData_ED%Gravity /)       ! "Gravitational acceleration" m/s^2
      
         ! now initialize BeamDyn for all beams
      dt_BD = p_FAST%dt_module( MODULE_BD )
                        
      Init%InData_BD%HubPos = ED%y%HubPtMotion%Position(:,1)
      Init%InData_BD%HubRot = ED%y%HubPtMotion%RefOrientation(:,:,1)
            
      p_FAST%BD_OutputSibling = .true.
      
      allocate( y_FAST%Lin%Modules(MODULE_BD)%Instance(p_FAST%nBeams), stat=ErrStat2)
      if (ErrStat2 /= 0 ) then
         call SetErrStat(ErrID_Fatal, "Error allocating Lin%Modules(BD).", ErrStat, ErrMsg, RoutineName )
         CALL Cleanup()
         RETURN
      end if

      DO k=1,p_FAST%nBeams
         Init%InData_BD%RootName     = TRIM(p_FAST%OutFileRoot)//'.'//TRIM(y_FAST%Module_Abrev(Module_BD))//TRIM( Num2LStr(k) )
         
         
         Init%InData_BD%InputFile    = p_FAST%BDBldFile(k)
         
         Init%InData_BD%GlbPos       = ED%y%BladeRootMotion(k)%Position(:,1)          ! {:}    - - "Initial Position Vector of the local blade coordinate system"
         Init%InData_BD%GlbRot       = ED%y%BladeRootMotion(k)%RefOrientation(:,:,1)  ! {:}{:} - - "Initial direction cosine matrix of the local blade coordinate system"
         
         Init%InData_BD%RootDisp     = ED%y%BladeRootMotion(k)%TranslationDisp(:,1)   ! {:}    - - "Initial root displacement"
         Init%InData_BD%RootOri      = ED%y%BladeRootMotion(k)%Orientation(:,:,1)     ! {:}{:} - - "Initial root orientation"
         Init%InData_BD%RootVel(1:3) = ED%y%BladeRootMotion(k)%TranslationVel(:,1)    ! {:}    - - "Initial root velocities and angular veolcities"                  
         Init%InData_BD%RootVel(4:6) = ED%y%BladeRootMotion(k)%RotationVel(:,1)       ! {:}    - - "Initial root velocities and angular veolcities"                  
                           
         CALL BD_Init( Init%InData_BD, BD%Input(1,k), BD%p(k),  BD%x(k,STATE_CURR), BD%xd(k,STATE_CURR), BD%z(k,STATE_CURR), &
                           BD%OtherSt(k,STATE_CURR), BD%y(k),  BD%m(k), dt_BD, Init%OutData_BD(k), ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)  
            
         !bjj: we're going to force this to have the same timestep because I don't want to have to deal with n BD modules with n timesteps.
         IF ( k == 1 ) THEN
            p_FAST%dt_module( MODULE_BD ) = dt_BD
            
            p_FAST%ModuleInitialized(Module_BD) = .TRUE. ! this really should be once per BD instance, but BD doesn't care so I won't go through the effort to track this
            CALL SetModuleSubstepTime(Module_BD, p_FAST, y_FAST, ErrStat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)            
         ELSEIF ( .NOT. EqualRealNos( p_FAST%dt_module( MODULE_BD ),dt_BD )) THEN
            CALL SetErrStat(ErrID_Fatal,"All instances of BeamDyn (one per blade) must have the same time step.",ErrStat,ErrMsg,RoutineName)
         END IF

            ! We're going to do fewer computations if the BD input and output meshes that couple to AD are siblings:
         if (BD%p(k)%BldMotionNodeLoc /= BD_MESH_QP) p_FAST%BD_OutputSibling = .false.
      
         if (ErrStat>=AbortErrLev) exit !exit this loop so we don't get p_FAST%nBeams of the same errors
         
         if (allocated(Init%OutData_BD(k)%LinNames_y)) call move_alloc(Init%OutData_BD(k)%LinNames_y, y_FAST%Lin%Modules(MODULE_BD)%Instance(k)%Names_y )
         if (allocated(Init%OutData_BD(k)%LinNames_x)) call move_alloc(Init%OutData_BD(k)%LinNames_x, y_FAST%Lin%Modules(MODULE_BD)%Instance(k)%Names_x )
         if (allocated(Init%OutData_BD(k)%LinNames_u)) call move_alloc(Init%OutData_BD(k)%LinNames_u, y_FAST%Lin%Modules(MODULE_BD)%Instance(k)%Names_u )
         if (allocated(Init%OutData_BD(k)%RotFrame_y)) call move_alloc(Init%OutData_BD(k)%RotFrame_y, y_FAST%Lin%Modules(MODULE_BD)%Instance(k)%RotFrame_y )
         if (allocated(Init%OutData_BD(k)%RotFrame_x)) call move_alloc(Init%OutData_BD(k)%RotFrame_x, y_FAST%Lin%Modules(MODULE_BD)%Instance(k)%RotFrame_x )
         if (allocated(Init%OutData_BD(k)%RotFrame_u)) call move_alloc(Init%OutData_BD(k)%RotFrame_u, y_FAST%Lin%Modules(MODULE_BD)%Instance(k)%RotFrame_u )
         if (allocated(Init%OutData_BD(k)%IsLoad_u  )) call move_alloc(Init%OutData_BD(k)%IsLoad_u  , y_FAST%Lin%Modules(MODULE_BD)%Instance(k)%IsLoad_u   )
         if (allocated(Init%OutData_BD(k)%DerivOrder_x  )) call move_alloc(Init%OutData_BD(k)%DerivOrder_x  , y_FAST%Lin%Modules(MODULE_BD)%Instance(k)%DerivOrder_x   )
         
         if (allocated(Init%OutData_BD(k)%WriteOutputHdr)) y_FAST%Lin%Modules(MODULE_BD)%Instance(k)%NumOutputs = size(Init%OutData_BD(k)%WriteOutputHdr)
         
      END DO
               
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF           
      
   END IF         
      

   ! ........................
   ! initialize AeroDyn 
   ! ........................
   ALLOCATE( AD14%Input( p_FAST%InterpOrder+1 ), AD14%InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating AD14%Input and AD14%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF
     
   ALLOCATE( AD%Input( p_FAST%InterpOrder+1 ), AD%InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating AD%Input and AD%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF
      
      
   IF ( p_FAST%CompAero == Module_AD14 ) THEN
               
      CALL AD_SetInitInput(Init%InData_AD14, Init%OutData_ED, ED%y, p_FAST, ErrStat2, ErrMsg2)            ! set the values in Init%InData_AD14
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
                                       
      CALL AD14_Init( Init%InData_AD14, AD14%Input(1), AD14%p, AD14%x(STATE_CURR), AD14%xd(STATE_CURR), AD14%z(STATE_CURR), &
                     AD14%OtherSt(STATE_CURR), AD14%y, AD14%m, p_FAST%dt_module( MODULE_AD14 ), Init%OutData_AD14, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      p_FAST%ModuleInitialized(Module_AD14) = .TRUE.            
      CALL SetModuleSubstepTime(Module_AD14, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
         ! bjj: this really shouldn't be in the FAST glue code, but I'm going to put this check here so people don't use an invalid model 
         !    and send me emails to debug numerical issues in their results.
      IF ( AD14%p%TwrProps%PJM_Version .AND. p_FAST%TurbineType == Type_Offshore_Floating ) THEN
         CALL SetErrStat(ErrID_Fatal,'AeroDyn v14 tower influence model "NEWTOWER" is invalid for models of floating offshore turbines.',ErrStat,ErrMsg,RoutineName)
      END IF         
            
      AirDens = Init%OutData_AD14%AirDens
      
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF       
      
   ELSEIF ( p_FAST%CompAero == Module_AD ) THEN
      
      
         ! set initialization data for AD
      CALL AllocAry( Init%InData_AD%BladeRootPosition,      3, Init%OutData_ED%NumBl, 'Init%InData_AD%BladeRootPosition', errStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL AllocAry( Init%InData_AD%BladeRootOrientation,3, 3, Init%OutData_ED%NumBl, 'Init%InData_AD%BladeRootOrientation', errStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF (ErrStat >= AbortErrLev) THEN
            CALL Cleanup()
            RETURN
         END IF
      Init%InData_AD%Gravity            = Init%OutData_ED%Gravity      
      Init%InData_AD%Linearize          = p_FAST%Linearize
      Init%InData_AD%InputFile          = p_FAST%AeroFile
      Init%InData_AD%NumBlades          = Init%OutData_ED%NumBl
      Init%InData_AD%RootName           = p_FAST%OutFileRoot
      Init%InData_AD%HubPosition        = ED%y%HubPtMotion%Position(:,1)
      Init%InData_AD%HubOrientation     = ED%y%HubPtMotion%RefOrientation(:,:,1)
      
      do k=1,Init%OutData_ED%NumBl
         Init%InData_AD%BladeRootPosition(:,k)      = ED%y%BladeRootMotion(k)%Position(:,1)
         Init%InData_AD%BladeRootOrientation(:,:,k) = ED%y%BladeRootMotion(k)%RefOrientation(:,:,1)
      end do
      
            
      CALL AD_Init( Init%InData_AD, AD%Input(1), AD%p, AD%x(STATE_CURR), AD%xd(STATE_CURR), AD%z(STATE_CURR), &
                    AD%OtherSt(STATE_CURR), AD%y, AD%m, p_FAST%dt_module( MODULE_AD ), Init%OutData_AD, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      p_FAST%ModuleInitialized(Module_AD) = .TRUE.            
      CALL SetModuleSubstepTime(Module_AD, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
                               
      allocate( y_FAST%Lin%Modules(MODULE_AD)%Instance(1), stat=ErrStat2)
      if (ErrStat2 /= 0 ) then
         call SetErrStat(ErrID_Fatal, "Error allocating Lin%Modules(AD).", ErrStat, ErrMsg, RoutineName )
      else
         if (allocated(Init%OutData_AD%LinNames_u  )) call move_alloc(Init%OutData_AD%LinNames_u  ,y_FAST%Lin%Modules(MODULE_AD)%Instance(1)%Names_u )
         if (allocated(Init%OutData_AD%LinNames_y  )) call move_alloc(Init%OutData_AD%LinNames_y  ,y_FAST%Lin%Modules(MODULE_AD)%Instance(1)%Names_y )
         if (allocated(Init%OutData_AD%LinNames_x  )) call move_alloc(Init%OutData_AD%LinNames_x  ,y_FAST%Lin%Modules(MODULE_AD)%Instance(1)%Names_x )
         if (allocated(Init%OutData_AD%RotFrame_u  )) call move_alloc(Init%OutData_AD%RotFrame_u  ,y_FAST%Lin%Modules(MODULE_AD)%Instance(1)%RotFrame_u )
         if (allocated(Init%OutData_AD%RotFrame_y  )) call move_alloc(Init%OutData_AD%RotFrame_y  ,y_FAST%Lin%Modules(MODULE_AD)%Instance(1)%RotFrame_y )
         if (allocated(Init%OutData_AD%RotFrame_x  )) call move_alloc(Init%OutData_AD%RotFrame_x  ,y_FAST%Lin%Modules(MODULE_AD)%Instance(1)%RotFrame_x )
         if (allocated(Init%OutData_AD%IsLoad_u    )) call move_alloc(Init%OutData_AD%IsLoad_u    ,y_FAST%Lin%Modules(MODULE_AD)%Instance(1)%IsLoad_u   )
         if (allocated(Init%OutData_AD%DerivOrder_x)) call move_alloc(Init%OutData_AD%DerivOrder_x,y_FAST%Lin%Modules(MODULE_AD)%Instance(1)%DerivOrder_x )

         if (allocated(Init%OutData_AD%WriteOutputHdr)) y_FAST%Lin%Modules(MODULE_AD)%Instance(1)%NumOutputs = size(Init%OutData_AD%WriteOutputHdr)
      end if
      
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF       
      
      AirDens = Init%OutData_AD%AirDens
      
   ELSE
      AirDens = 0.0_ReKi
   END IF ! CompAero
   
               
   ! ........................
   ! initialize InflowWind
   ! ........................   
   ALLOCATE( IfW%Input( p_FAST%InterpOrder+1 ), IfW%InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating IfW%Input and IfW%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF
              
   IF ( p_FAST%CompInflow == Module_IfW ) THEN
      
      Init%InData_IfW%Linearize        = p_FAST%Linearize
      Init%InData_IfW%InputFileName    = p_FAST%InflowFile
      Init%InData_IfW%RootName         = TRIM(p_FAST%OutFileRoot)//'.'//TRIM(y_FAST%Module_Abrev(Module_IfW))
      Init%InData_IfW%UseInputFile     = .TRUE.
   
      Init%InData_IfW%NumWindPoints = 0      
      IF ( p_FAST%CompServo == Module_SrvD ) Init%InData_IfW%NumWindPoints = Init%InData_IfW%NumWindPoints + 1
      IF ( p_FAST%CompAero  == Module_AD14 ) THEN
         Init%InData_IfW%NumWindPoints = Init%InData_IfW%NumWindPoints + Init%OutData_ED%NumBl * AD14%Input(1)%InputMarkers(1)%NNodes + AD14%Input(1)%Twr_InputMarkers%NNodes
      ELSEIF ( p_FAST%CompAero  == Module_AD ) THEN
         Init%InData_IfW%NumWindPoints = Init%InData_IfW%NumWindPoints + AD%Input(1)%TowerMotion%NNodes
         DO k=1,Init%OutData_ED%NumBl
            Init%InData_IfW%NumWindPoints = Init%InData_IfW%NumWindPoints + AD%Input(1)%BladeMotion(k)%NNodes
         END DO
         if (allocated(AD%OtherSt(STATE_CURR)%WakeLocationPoints)) then
            Init%InData_IfW%NumWindPoints = Init%InData_IfW%NumWindPoints + size(AD%OtherSt(STATE_CURR)%WakeLocationPoints,DIM=2)
         end if
      END IF
      
      ! lidar        
      Init%InData_IfW%lidar%Tmax                   = p_FAST%TMax
      Init%InData_IfW%lidar%HubPosition            = ED%y%HubPtMotion%Position(:,1) 
      
      IF ( PRESENT(ExternInitData) ) THEN
         Init%InData_IfW%Use4Dext = ExternInitData%FarmIntegration

         if (Init%InData_IfW%Use4Dext) then
            Init%InData_IfW%FDext%n      = ExternInitData%windGrid_n
            Init%InData_IfW%FDext%delta  = ExternInitData%windGrid_delta
            Init%InData_IfW%FDext%pZero  = ExternInitData%windGrid_pZero
         end if
         
         ! bjj: these lidar inputs should come from an InflowWind input file; I'm hard coding them here for now
         Init%InData_IfW%lidar%SensorType          = ExternInitData%SensorType   
         Init%InData_IfW%lidar%LidRadialVel        = ExternInitData%LidRadialVel   
         Init%InData_IfW%lidar%RotorApexOffsetPos  = 0.0         
         Init%InData_IfW%lidar%NumPulseGate        = 0
      ELSE
         Init%InData_IfW%lidar%SensorType          = SensorType_None
         Init%InData_IfW%Use4Dext                  = .false.
      END IF
                                     
      CALL InflowWind_Init( Init%InData_IfW, IfW%Input(1), IfW%p, IfW%x(STATE_CURR), IfW%xd(STATE_CURR), IfW%z(STATE_CURR),  &
                     IfW%OtherSt(STATE_CURR), IfW%y, IfW%m, p_FAST%dt_module( MODULE_IfW ), Init%OutData_IfW, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      p_FAST%ModuleInitialized(Module_IfW) = .TRUE.            
      CALL SetModuleSubstepTime(Module_IfW, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
      allocate( y_FAST%Lin%Modules(MODULE_IfW)%Instance(1), stat=ErrStat2)
      if (ErrStat2 /= 0 ) then
         call SetErrStat(ErrID_Fatal, "Error allocating Lin%Modules(IfW).", ErrStat, ErrMsg, RoutineName )
      else
         if (allocated(Init%OutData_IfW%LinNames_y)) call move_alloc(Init%OutData_IfW%LinNames_y,y_FAST%Lin%Modules(MODULE_IfW)%Instance(1)%Names_y )
         if (allocated(Init%OutData_IfW%LinNames_u)) call move_alloc(Init%OutData_IfW%LinNames_u,y_FAST%Lin%Modules(MODULE_IfW)%Instance(1)%Names_u )
         if (allocated(Init%OutData_IfW%RotFrame_y)) call move_alloc(Init%OutData_IfW%RotFrame_y,y_FAST%Lin%Modules(MODULE_IfW)%Instance(1)%RotFrame_y )
         if (allocated(Init%OutData_IfW%RotFrame_u)) call move_alloc(Init%OutData_IfW%RotFrame_u,y_FAST%Lin%Modules(MODULE_IfW)%Instance(1)%RotFrame_u )
         if (allocated(Init%OutData_IfW%IsLoad_u  )) call move_alloc(Init%OutData_IfW%IsLoad_u  ,y_FAST%Lin%Modules(MODULE_IfW)%Instance(1)%IsLoad_u   )

         if (allocated(Init%OutData_IfW%WriteOutputHdr)) y_FAST%Lin%Modules(MODULE_IfW)%Instance(1)%NumOutputs = size(Init%OutData_IfW%WriteOutputHdr)
         y_FAST%Lin%WindSpeed = Init%OutData_IfW%WindFileInfo%MWS
      end if
      
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF       
      
   ELSEIF ( p_FAST%CompInflow == Module_OpFM ) THEN
      
      IF ( PRESENT(ExternInitData) ) THEN
         Init%InData_OpFM%NumSC2Ctrl = ExternInitData%NumSC2Ctrl
         Init%InData_OpFM%NumCtrl2SC = ExternInitData%NumCtrl2SC  
         Init%InData_OpFM%NumActForcePtsBlade = ExternInitData%NumActForcePtsBlade
         Init%InData_OpFM%NumActForcePtsTower = ExternInitData%NumActForcePtsTower 
      ELSE
         CALL SetErrStat( ErrID_Fatal, 'OpenFOAM integration can be used only with external input data (not the stand-alone executable).', ErrStat, ErrMsg, RoutineName )
         CALL Cleanup()
         RETURN         
      END IF
      Init%InData_OpFM%BladeLength = Init%OutData_ED%BladeLength
      Init%InData_OpFM%TowerHeight = Init%OutData_ED%TowerHeight
      Init%InData_OpFM%TowerBaseHeight = Init%OutData_ED%TowerBaseHeight
      ALLOCATE(Init%InData_OpFM%StructBldRNodes( SIZE(Init%OutData_ED%BldRNodes)),  STAT=ErrStat2)
      Init%InData_OpFM%StructBldRNodes(:) = Init%OutData_ED%BldRNodes(:)
      ALLOCATE(Init%InData_OpFM%StructTwrHNodes( SIZE(Init%OutData_ED%TwrHNodes)),  STAT=ErrStat2)
      Init%InData_OpFM%StructTwrHNodes(:) = Init%OutData_ED%TwrHNodes(:)
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating OpFM%InitInput.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF
         ! set up the data structures for integration with OpenFOAM
      CALL Init_OpFM( Init%InData_OpFM, p_FAST, AirDens, AD14%Input(1), AD%Input(1), Init%OutData_AD, AD%y, ED%y, OpFM, Init%OutData_OpFM, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF       
                  
      !bjj: fix me!!! to do
      Init%OutData_IfW%WindFileInfo%MWS = 0.0_ReKi
      
   ELSE
      Init%OutData_IfW%WindFileInfo%MWS = 0.0_ReKi
   END IF   ! CompInflow
   
   ! ........................
   ! initialize SuperController
   ! ........................   
      IF ( PRESENT(ExternInitData) ) THEN
         Init%InData_SC%NumSC2Ctrl = ExternInitData%NumSC2Ctrl
         Init%InData_SC%NumCtrl2SC = ExternInitData%NumCtrl2SC  
      ELSE
         Init%InData_SC%NumSC2Ctrl = 0
         Init%InData_SC%NumCtrl2SC = 0
      END IF
      
         ! set up the data structures for integration with supercontroller
      CALL Init_SC( Init%InData_SC, SC, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF       

   ! ........................
   ! some checks for AeroDyn14's Dynamic Inflow with Mean Wind Speed from InflowWind:
   ! (DO NOT COPY THIS CODE!)
   ! bjj: AeroDyn14 should not need this rule of thumb; it should check the instantaneous values when the code runs
   ! ........................   
   
   IF ( p_FAST%CompAero == Module_AD14 ) THEN
      IF (AD14%p%DynInfl) THEN               
         IF ( Init%OutData_IfW%WindFileInfo%MWS  < 8.0 ) THEN
            CALL SetErrStat(ErrID_Fatal,'AeroDyn v14 "DYNINFL" InfModel is invalid for models with wind speeds less than 8 m/s.',ErrStat,ErrMsg,RoutineName)
            !CALL SetErrStat(ErrID_Info,'Estimated average inflow wind speed is less than 8 m/s. Dynamic Inflow will be turned off.',ErrStat,ErrMess,RoutineName )
         END IF
      END IF      
   END IF
   
   
   ! ........................
   ! initialize ServoDyn 
   ! ........................
   ALLOCATE( SrvD%Input( p_FAST%InterpOrder+1 ), SrvD%InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating SrvD%Input and SrvD%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF
      
   IF ( p_FAST%CompServo == Module_SrvD ) THEN
      Init%InData_SrvD%InputFile     = p_FAST%ServoFile
      Init%InData_SrvD%RootName      = TRIM(p_FAST%OutFileRoot)//'.'//TRIM(y_FAST%Module_Abrev(Module_SrvD))
      Init%InData_SrvD%NumBl         = Init%OutData_ED%NumBl
      Init%InData_SrvD%gravity       = Init%OutData_ED%gravity
      Init%InData_SrvD%r_N_O_G       = ED%Input(1)%NacelleLoads%Position(:,1)
      Init%InData_SrvD%r_TwrBase     = Init%OutData_ED%TwrBasePos
      Init%InData_SrvD%TMax          = p_FAST%TMax
      Init%InData_SrvD%AirDens       = AirDens
      Init%InData_SrvD%AvgWindSpeed  = Init%OutData_IfW%WindFileInfo%MWS
      Init%InData_SrvD%Linearize     = p_FAST%Linearize
      Init%InData_SrvD%TrimCase      = p_FAST%TrimCase
      Init%InData_SrvD%TrimGain      = p_FAST%TrimGain
      Init%InData_SrvD%RotSpeedRef   = Init%OutData_ED%RotSpeed
      
      IF ( PRESENT(ExternInitData) ) THEN
         Init%InData_SrvD%NumSC2Ctrl = ExternInitData%NumSC2Ctrl
         Init%InData_SrvD%NumCtrl2SC = ExternInitData%NumCtrl2SC
      ELSE
         Init%InData_SrvD%NumSC2Ctrl = 0
         Init%InData_SrvD%NumCtrl2SC = 0
      END IF      
            
      CALL AllocAry(Init%InData_SrvD%BlPitchInit, Init%OutData_ED%NumBl, 'BlPitchInit', ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      if (ErrStat >= abortErrLev) then ! make sure allocatable arrays are valid before setting them
         CALL Cleanup()
         RETURN
      end if

      Init%InData_SrvD%BlPitchInit   = Init%OutData_ED%BlPitch
      CALL SrvD_Init( Init%InData_SrvD, SrvD%Input(1), SrvD%p, SrvD%x(STATE_CURR), SrvD%xd(STATE_CURR), SrvD%z(STATE_CURR), &
                      SrvD%OtherSt(STATE_CURR), SrvD%y, SrvD%m, p_FAST%dt_module( MODULE_SrvD ), Init%OutData_SrvD, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      p_FAST%ModuleInitialized(Module_SrvD) = .TRUE.

      !IF ( Init%OutData_SrvD%CouplingScheme == ExplicitLoose ) THEN ...  bjj: abort if we're doing anything else!

      CALL SetModuleSubstepTime(Module_SrvD, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      !! initialize SrvD%y%ElecPwr and SrvD%y%GenTq because they are one timestep different (used as input for the next step)?
                  
      allocate( y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1), stat=ErrStat2)
      if (ErrStat2 /= 0 ) then
         call SetErrStat(ErrID_Fatal, "Error allocating Lin%Modules(SrvD).", ErrStat, ErrMsg, RoutineName )
      else
         if (allocated(Init%OutData_SrvD%LinNames_y)) call move_alloc(Init%OutData_SrvD%LinNames_y,y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%Names_y )
         if (allocated(Init%OutData_SrvD%LinNames_u)) call move_alloc(Init%OutData_SrvD%LinNames_u,y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%Names_u )
         if (allocated(Init%OutData_SrvD%RotFrame_y)) call move_alloc(Init%OutData_SrvD%RotFrame_y,y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%RotFrame_y )
         if (allocated(Init%OutData_SrvD%RotFrame_u)) call move_alloc(Init%OutData_SrvD%RotFrame_u,y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%RotFrame_u )
         if (allocated(Init%OutData_SrvD%IsLoad_u  )) call move_alloc(Init%OutData_SrvD%IsLoad_u  ,y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%IsLoad_u   )

         if (allocated(Init%OutData_SrvD%WriteOutputHdr)) y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%NumOutputs = size(Init%OutData_SrvD%WriteOutputHdr)
      end if
      
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF
      
   ! ........................
   ! some checks for AeroDyn and ElastoDyn inputs with the high-speed shaft brake hack in ElastoDyn:
   ! (DO NOT COPY THIS CODE!)
   ! ........................   
         ! bjj: this is a hack to get high-speed shaft braking in FAST v8
      
      IF ( Init%OutData_SrvD%UseHSSBrake ) THEN
         IF ( p_FAST%CompAero == Module_AD14 ) THEN
            IF ( AD14%p%DYNINFL ) THEN
               CALL SetErrStat(ErrID_Fatal,'AeroDyn v14 "DYNINFL" InfModel is invalid for models with high-speed shaft braking.',ErrStat,ErrMsg,RoutineName)
            END IF
         END IF
         

         IF ( ED%p%method == Method_RK4 ) THEN ! bjj: should be using ElastoDyn's Method_ABM4 Method_AB4 parameters
            CALL SetErrStat(ErrID_Fatal,'ElastoDyn must use the AB4 or ABM4 integration method to implement high-speed shaft braking.',ErrStat,ErrMsg,RoutineName)
         ENDIF
      END IF ! Init%OutData_SrvD%UseHSSBrake
      
      
   END IF

   ! ........................
   ! set some VTK parameters required before HydroDyn init (so we can get wave elevations for visualization)
   ! ........................
   
      ! get wave elevation data for visualization
   if ( p_FAST%WrVTK > VTK_None ) then   
      call SetVTKParameters_B4HD(p_FAST, Init%OutData_ED, Init%InData_HD, BD, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF (ErrStat >= AbortErrLev) THEN
            CALL Cleanup()
            RETURN
         END IF       
   end if
   
   
   ! ........................
   ! initialize HydroDyn 
   ! ........................
   ALLOCATE( HD%Input( p_FAST%InterpOrder+1 ), HD%InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating HD%Input and HD%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF
      
   IF ( p_FAST%CompHydro == Module_HD ) THEN

      Init%InData_HD%Gravity       = Init%OutData_ED%Gravity
      Init%InData_HD%UseInputFile  = .TRUE.
      Init%InData_HD%InputFile     = p_FAST%HydroFile
      Init%InData_HD%OutRootName   = p_FAST%OutFileRoot
      Init%InData_HD%TMax          = p_FAST%TMax
      Init%InData_HD%hasIce        = p_FAST%CompIce /= Module_None
      Init%InData_HD%Linearize     = p_FAST%Linearize
      
         ! if wave field needs an offset, modify these values (added at request of SOWFA developers):
      Init%InData_HD%PtfmLocationX = p_FAST%TurbinePos(1) 
      Init%InData_HD%PtfmLocationY = p_FAST%TurbinePos(2)
      
      CALL HydroDyn_Init( Init%InData_HD, HD%Input(1), HD%p,  HD%x(STATE_CURR), HD%xd(STATE_CURR), HD%z(STATE_CURR), &
                          HD%OtherSt(STATE_CURR), HD%y, HD%m, p_FAST%dt_module( MODULE_HD ), Init%OutData_HD, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
      p_FAST%ModuleInitialized(Module_HD) = .TRUE.
      CALL SetModuleSubstepTime(Module_HD, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
      allocate( y_FAST%Lin%Modules(MODULE_HD)%Instance(1), stat=ErrStat2)
      if (ErrStat2 /= 0 ) then
         call SetErrStat(ErrID_Fatal, "Error allocating Lin%Modules(HD).", ErrStat, ErrMsg, RoutineName )
      else
         if (allocated(Init%OutData_HD%LinNames_y)) call move_alloc(Init%OutData_HD%LinNames_y,y_FAST%Lin%Modules(MODULE_HD)%Instance(1)%Names_y )
         if (allocated(Init%OutData_HD%LinNames_u)) call move_alloc(Init%OutData_HD%LinNames_u,y_FAST%Lin%Modules(MODULE_HD)%Instance(1)%Names_u )
         if (allocated(Init%OutData_HD%LinNames_x)) call move_alloc(Init%OutData_HD%LinNames_x, y_FAST%Lin%Modules(MODULE_HD)%Instance(1)%Names_x )
! LIN-TODO: Determine if we need to create this data even though we don't have rotating frames in HD
         !if (allocated(Init%OutData_HD%RotFrame_y)) call move_alloc(Init%OutData_HD%RotFrame_y,y_FAST%Lin%Modules(MODULE_HD)%Instance(1)%RotFrame_y )
         !if (allocated(Init%OutData_HD%RotFrame_u)) call move_alloc(Init%OutData_HD%RotFrame_u,y_FAST%Lin%Modules(MODULE_HD)%Instance(1)%RotFrame_u )
         if (allocated(Init%OutData_HD%DerivOrder_x)) call move_alloc(Init%OutData_HD%DerivOrder_x,y_FAST%Lin%Modules(MODULE_HD)%Instance(1)%DerivOrder_x)
         if (allocated(Init%OutData_HD%IsLoad_u  )) call move_alloc(Init%OutData_HD%IsLoad_u  ,y_FAST%Lin%Modules(MODULE_HD)%Instance(1)%IsLoad_u   )

         if (allocated(Init%OutData_HD%WriteOutputHdr)) y_FAST%Lin%Modules(MODULE_HD)%Instance(1)%NumOutputs = size(Init%OutData_HD%WriteOutputHdr)
      end if
     
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF       
   END IF   ! CompHydro

   ! ........................
   ! initialize SubDyn or ExtPtfm_MCKF
   ! ........................
   ALLOCATE( SD%Input( p_FAST%InterpOrder+1 ), SD%InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating SD%Input and SD%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF

   ALLOCATE( ExtPtfm%Input( p_FAST%InterpOrder+1 ), ExtPtfm%InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating ExtPtfm%Input and ExtPtfm%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF
      
   IF ( p_FAST%CompSub == Module_SD ) THEN
          
      IF ( p_FAST%CompHydro == Module_HD ) THEN
         Init%InData_SD%WtrDpth = Init%OutData_HD%WtrDpth
      ELSE
         Init%InData_SD%WtrDpth = 0.0_ReKi
      END IF
            
      Init%InData_SD%g             = Init%OutData_ED%Gravity     
      !Init%InData_SD%UseInputFile = .TRUE. 
      Init%InData_SD%SDInputFile   = p_FAST%SubFile
      Init%InData_SD%RootName      = p_FAST%OutFileRoot
      Init%InData_SD%TP_RefPoint   = ED%y%PlatformPtMesh%Position(:,1)  ! bjj: not sure what this is supposed to be 
      Init%InData_SD%SubRotateZ    = 0.0                                ! bjj: not sure what this is supposed to be 
      
            
      CALL SD_Init( Init%InData_SD, SD%Input(1), SD%p,  SD%x(STATE_CURR), SD%xd(STATE_CURR), SD%z(STATE_CURR),  &
                    SD%OtherSt(STATE_CURR), SD%y, SD%m, p_FAST%dt_module( MODULE_SD ), Init%OutData_SD, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
      p_FAST%ModuleInitialized(Module_SD) = .TRUE.
      CALL SetModuleSubstepTime(Module_SD, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
               
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF   
   ELSE IF ( p_FAST%CompSub == Module_ExtPtfm ) THEN

      Init%InData_ExtPtfm%InputFile = p_FAST%SubFile
      Init%InData_ExtPtfm%RootName  = trim(p_FAST%OutFileRoot)//'.'//TRIM(y_FAST%Module_Abrev(Module_ExtPtfm))
      Init%InData_ExtPtfm%Linearize = p_FAST%Linearize
      Init%InData_ExtPtfm%PtfmRefzt = ED%p%PtfmRefzt ! Required
      
      CALL ExtPtfm_Init( Init%InData_ExtPtfm, ExtPtfm%Input(1), ExtPtfm%p,  &
                         ExtPtfm%x(STATE_CURR), ExtPtfm%xd(STATE_CURR), ExtPtfm%z(STATE_CURR),  ExtPtfm%OtherSt(STATE_CURR), &
                         ExtPtfm%y, ExtPtfm%m, p_FAST%dt_module( MODULE_ExtPtfm ), Init%OutData_ExtPtfm, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
      p_FAST%ModuleInitialized(MODULE_ExtPtfm) = .TRUE.
      CALL SetModuleSubstepTime(MODULE_ExtPtfm, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      allocate( y_FAST%Lin%Modules(MODULE_ExtPtfm)%Instance(1), stat=ErrStat2)
      if (ErrStat2 /= 0 ) then
         call SetErrStat(ErrID_Fatal, "Error allocating Lin%Modules(ExtPtfm).", ErrStat, ErrMsg, RoutineName )
      else
         if (allocated(Init%OutData_ExtPtfm%LinNames_y)) call move_alloc(Init%OutData_ExtPtfm%LinNames_y,y_FAST%Lin%Modules(MODULE_ExtPtfm)%Instance(1)%Names_y)
         if (allocated(Init%OutData_ExtPtfm%LinNames_x)) call move_alloc(Init%OutData_ExtPtfm%LinNames_x,y_FAST%Lin%Modules(MODULE_ExtPtfm)%Instance(1)%Names_x)
         if (allocated(Init%OutData_ExtPtfm%LinNames_u)) call move_alloc(Init%OutData_ExtPtfm%LinNames_u,y_FAST%Lin%Modules(MODULE_ExtPtfm)%Instance(1)%Names_u)
         if (allocated(Init%OutData_ExtPtfm%RotFrame_y)) call move_alloc(Init%OutData_ExtPtfm%RotFrame_y,y_FAST%Lin%Modules(MODULE_ExtPtfm)%Instance(1)%RotFrame_y)
         if (allocated(Init%OutData_ExtPtfm%RotFrame_x)) call move_alloc(Init%OutData_ExtPtfm%RotFrame_x,y_FAST%Lin%Modules(MODULE_ExtPtfm)%Instance(1)%RotFrame_x)
         if (allocated(Init%OutData_ExtPtfm%RotFrame_u)) call move_alloc(Init%OutData_ExtPtfm%RotFrame_u,y_FAST%Lin%Modules(MODULE_ExtPtfm)%Instance(1)%RotFrame_u)
         if (allocated(Init%OutData_ExtPtfm%IsLoad_u  )) call move_alloc(Init%OutData_ExtPtfm%IsLoad_u  ,y_FAST%Lin%Modules(MODULE_ExtPtfm)%Instance(1)%IsLoad_u  )
         if (allocated(Init%OutData_ExtPtfm%WriteOutputHdr)) y_FAST%Lin%Modules(MODULE_ExtPtfm)%Instance(1)%NumOutputs = size(Init%OutData_ExtPtfm%WriteOutputHdr)
      end if
               
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF   
      
   END IF

   ! ------------------------------
   ! initialize CompMooring modules 
   ! ------------------------------
   ALLOCATE( MAPp%Input( p_FAST%InterpOrder+1 ), MAPp%InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating MAPp%Input and MAPp%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF
   ALLOCATE( MD%Input( p_FAST%InterpOrder+1 ), MD%InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating MD%Input and MD%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF   
   ALLOCATE( FEAM%Input( p_FAST%InterpOrder+1 ), FEAM%InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating FEAM%Input and FEAM%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF   
   ALLOCATE( Orca%Input( p_FAST%InterpOrder+1 ), Orca%InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating Orca%Input and Orca%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF   
      
   ! ........................
   ! initialize MAP 
   ! ........................
   IF (p_FAST%CompMooring == Module_MAP) THEN
      !bjj: until we modify this, MAP requires HydroDyn to be used. (perhaps we could send air density from AeroDyn or something...)
      
      CALL WrScr(NewLine) !bjj: I'm printing two blank lines here because MAP seems to be writing over the last line on the screen.
      

!      Init%InData_MAP%rootname          =  p_FAST%OutFileRoot        ! Output file name 
      Init%InData_MAP%gravity           =  Init%OutData_ED%Gravity    ! This need to be according to g used in ElastoDyn
      Init%InData_MAP%sea_density       =  Init%OutData_HD%WtrDens    ! This needs to be set according to seawater density in HydroDyn
      Init%InData_MAP%depth             =  Init%OutData_HD%WtrDpth    ! This need to be set according to the water depth in HydroDyn
                  
   ! differences for MAP++
      Init%InData_MAP%file_name         =  p_FAST%MooringFile        ! This needs to be set according to what is in the FAST input file. 
      Init%InData_MAP%summary_file_name =  TRIM(p_FAST%OutFileRoot)//'.MAP.sum'        ! Output file name 
      Init%InData_MAP%depth             = -Init%OutData_HD%WtrDpth    ! This need to be set according to the water depth in HydroDyn
            
      Init%InData_MAP%LinInitInp%Linearize = p_FAST%Linearize   
      
      CALL MAP_Init( Init%InData_MAP, MAPp%Input(1), MAPp%p,  MAPp%x(STATE_CURR), MAPp%xd(STATE_CURR), MAPp%z(STATE_CURR), MAPp%OtherSt, &
                      MAPp%y, p_FAST%dt_module( MODULE_MAP ), Init%OutData_MAP, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
      p_FAST%ModuleInitialized(Module_MAP) = .TRUE.
      CALL SetModuleSubstepTime(Module_MAP, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
      allocate( y_FAST%Lin%Modules(Module_MAP)%Instance(1), stat=ErrStat2)
      if (ErrStat2 /= 0 ) then
         call SetErrStat(ErrID_Fatal, "Error allocating Lin%Modules(MAP).", ErrStat, ErrMsg, RoutineName )
      else
         if (allocated(Init%OutData_MAP%LinInitOut%LinNames_y)) call move_alloc(Init%OutData_MAP%LinInitOut%LinNames_y,y_FAST%Lin%Modules(Module_MAP)%Instance(1)%Names_y )
         if (allocated(Init%OutData_MAP%LinInitOut%LinNames_u)) call move_alloc(Init%OutData_MAP%LinInitOut%LinNames_u,y_FAST%Lin%Modules(Module_MAP)%Instance(1)%Names_u )
! LIN-TODO: Determine if we need to create this data even though we don't have rotating frames in MAP
         !if (allocated(Init%OutData_MAP%LinInitOut%RotFrame_y)) call move_alloc(Init%OutData_MAP%LinInitOut%RotFrame_y,y_FAST%Lin%Modules(Module_MAP)%Instance(1)%RotFrame_y )
         !if (allocated(Init%OutData_MAP%LinInitOut%RotFrame_u)) call move_alloc(Init%OutData_MAP%LinInitOut%RotFrame_u,y_FAST%Lin%Modules(Module_MAP)%Instance(1)%RotFrame_u )
         if (allocated(Init%OutData_MAP%LinInitOut%IsLoad_u  )) call move_alloc(Init%OutData_MAP%LinInitOut%IsLoad_u  ,y_FAST%Lin%Modules(Module_MAP)%Instance(1)%IsLoad_u   )

         if (allocated(Init%OutData_MAP%WriteOutputHdr)) y_FAST%Lin%Modules(Module_MAP)%Instance(1)%NumOutputs = size(Init%OutData_MAP%WriteOutputHdr)
      end if
      
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF              
   ! ........................
   ! initialize MoorDyn 
   ! ........................
   ELSEIF (p_FAST%CompMooring == Module_MD) THEN
                        
      Init%InData_MD%FileName  = p_FAST%MooringFile         ! This needs to be set according to what is in the FAST input file. 
      Init%InData_MD%RootName  = p_FAST%OutFileRoot
      
      Init%InData_MD%PtfmInit  = Init%OutData_ED%PlatformPos !ED%x(STATE_CURR)%QT(1:6)   ! initial position of the platform !bjj: this should come from Init%OutData_ED, not x_ED
      Init%InData_MD%g         = Init%OutData_ED%Gravity     ! This need to be according to g used in ElastoDyn 
      Init%InData_MD%rhoW      = Init%OutData_HD%WtrDens     ! This needs to be set according to seawater density in HydroDyn      
      Init%InData_MD%WtrDepth  = Init%OutData_HD%WtrDpth    ! This need to be set according to the water depth in HydroDyn
            
      CALL MD_Init( Init%InData_MD, MD%Input(1), MD%p, MD%x(STATE_CURR), MD%xd(STATE_CURR), MD%z(STATE_CURR), &
                    MD%OtherSt(STATE_CURR), MD%y, MD%m, p_FAST%dt_module( MODULE_MD ), Init%OutData_MD, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
      p_FAST%ModuleInitialized(Module_MD) = .TRUE.
      CALL SetModuleSubstepTime(Module_MD, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF
   ! ........................
   ! initialize FEAM 
   ! ........................
   ELSEIF (p_FAST%CompMooring == Module_FEAM) THEN
            
      Init%InData_FEAM%InputFile   = p_FAST%MooringFile         ! This needs to be set according to what is in the FAST input file. 
      Init%InData_FEAM%RootName    = TRIM(p_FAST%OutFileRoot)//'.'//TRIM(y_FAST%Module_Abrev(Module_FEAM))
      
      Init%InData_FEAM%PtfmInit    = Init%OutData_ED%PlatformPos !ED%x(STATE_CURR)%QT(1:6)   ! initial position of the platform !bjj: this should come from Init%OutData_ED, not x_ED
      Init%InData_FEAM%NStepWave   = 1                          ! an arbitrary number > 0 (to set the size of the wave data, which currently contains all zero values)     
      Init%InData_FEAM%gravity     = Init%OutData_ED%Gravity     ! This need to be according to g used in ElastoDyn 
      Init%InData_FEAM%WtrDens     = Init%OutData_HD%WtrDens     ! This needs to be set according to seawater density in HydroDyn      
!      Init%InData_FEAM%depth       =  Init%OutData_HD%WtrDpth    ! This need to be set according to the water depth in HydroDyn
            
      CALL FEAM_Init( Init%InData_FEAM, FEAM%Input(1), FEAM%p,  FEAM%x(STATE_CURR), FEAM%xd(STATE_CURR), FEAM%z(STATE_CURR), &
                      FEAM%OtherSt(STATE_CURR), FEAM%y, FEAM%m, p_FAST%dt_module( MODULE_FEAM ), Init%OutData_FEAM, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
      p_FAST%ModuleInitialized(Module_FEAM) = .TRUE.
      CALL SetModuleSubstepTime(Module_FEAM, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF              
   ! ........................
   ! initialize OrcaFlex Interface 
   ! ........................
   ELSEIF (p_FAST%CompMooring == Module_Orca) THEN
            
      Init%InData_Orca%InputFile = p_FAST%MooringFile
      Init%InData_Orca%RootName  = p_FAST%OutFileRoot
      Init%InData_Orca%TMax      = p_FAST%TMax 
                  
      CALL Orca_Init( Init%InData_Orca, Orca%Input(1), Orca%p,  Orca%x(STATE_CURR), Orca%xd(STATE_CURR), Orca%z(STATE_CURR), Orca%OtherSt(STATE_CURR), &
                      Orca%y, Orca%m, p_FAST%dt_module( MODULE_Orca ), Init%OutData_Orca, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
      p_FAST%ModuleInitialized(MODULE_Orca) = .TRUE.
      CALL SetModuleSubstepTime(MODULE_Orca, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF
   END IF

   ! ------------------------------
   ! initialize CompIce modules 
   ! ------------------------------
   ALLOCATE( IceF%Input( p_FAST%InterpOrder+1 ), IceF%InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating IceF%Input and IceF%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF  
      
      ! We need this to be allocated (else we have issues passing nonallocated arrays and using the first index of Input(),
      !   but we don't need the space of IceD_MaxLegs if we're not using it. 
   IF ( p_FAST%CompIce /= Module_IceD ) THEN   
      IceDim = 1
   ELSE
      IceDim = IceD_MaxLegs
   END IF
      
      ! because there may be multiple instances of IceDyn, we'll allocate arrays for that here
      ! we could allocate these after 
   ALLOCATE( IceD%Input( p_FAST%InterpOrder+1, IceDim ), IceD%InputTimes( p_FAST%InterpOrder+1, IceDim ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating IceD%Input and IceD%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF  
      
     ALLOCATE( IceD%x(           IceDim,2), &
               IceD%xd(          IceDim,2), &
               IceD%z(           IceDim,2), &
               IceD%OtherSt(     IceDim,2), &
               IceD%p(           IceDim  ), &
               IceD%u(           IceDim  ), &
               IceD%y(           IceDim  ), &
               IceD%m(           IceDim  ), &
                                             STAT = ErrStat2 )                                                  
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating IceD state, input, and output data.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF      
         
         
   ! ........................
   ! initialize IceFloe 
   ! ........................
   IF ( p_FAST%CompIce == Module_IceF ) THEN
                      
      Init%InData_IceF%InputFile     = p_FAST%IceFile
      Init%InData_IceF%RootName      = TRIM(p_FAST%OutFileRoot)//'.'//TRIM(y_FAST%Module_Abrev(Module_IceF))     
      Init%InData_IceF%simLength     = p_FAST%TMax  !bjj: IceFloe stores this as single-precision (ReKi) TMax is DbKi
      Init%InData_IceF%MSL2SWL       = Init%OutData_HD%MSL2SWL
      Init%InData_IceF%gravity       = Init%OutData_ED%Gravity
      
      CALL IceFloe_Init( Init%InData_IceF, IceF%Input(1), IceF%p,  IceF%x(STATE_CURR), IceF%xd(STATE_CURR), IceF%z(STATE_CURR), &
                         IceF%OtherSt(STATE_CURR), IceF%y, IceF%m, p_FAST%dt_module( MODULE_IceF ), Init%OutData_IceF, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
      p_FAST%ModuleInitialized(Module_IceF) = .TRUE.
      CALL SetModuleSubstepTime(Module_IceF, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
              
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF              
   ! ........................
   ! initialize IceDyn 
   ! ........................
   ELSEIF ( p_FAST%CompIce == Module_IceD ) THEN  
      
      Init%InData_IceD%InputFile     = p_FAST%IceFile
      Init%InData_IceD%RootName      = TRIM(p_FAST%OutFileRoot)//'.'//TRIM(y_FAST%Module_Abrev(Module_IceD))//'1'     
      Init%InData_IceD%MSL2SWL       = Init%OutData_HD%MSL2SWL      
      Init%InData_IceD%WtrDens       = Init%OutData_HD%WtrDens    
      Init%InData_IceD%gravity       = Init%OutData_ED%Gravity
      Init%InData_IceD%TMax          = p_FAST%TMax
      Init%InData_IceD%LegNum        = 1
      
      CALL IceD_Init( Init%InData_IceD, IceD%Input(1,1), IceD%p(1),  IceD%x(1,STATE_CURR), IceD%xd(1,STATE_CURR), IceD%z(1,STATE_CURR), &
                      IceD%OtherSt(1,STATE_CURR), IceD%y(1), IceD%m(1), p_FAST%dt_module( MODULE_IceD ), Init%OutData_IceD, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
      p_FAST%ModuleInitialized(Module_IceD) = .TRUE.
      CALL SetModuleSubstepTime(Module_IceD, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
         ! now initialize IceD for additional legs (if necessary)
      dt_IceD           = p_FAST%dt_module( MODULE_IceD )
      p_FAST%numIceLegs = Init%OutData_IceD%numLegs     
      
      IF (p_FAST%numIceLegs > IceD_MaxLegs) THEN
         CALL SetErrStat(ErrID_Fatal,'IceDyn-FAST coupling is supported for up to '//TRIM(Num2LStr(IceD_MaxLegs))//' legs, but ' &
                           //TRIM(Num2LStr(p_FAST%numIceLegs))//' legs were specified.',ErrStat,ErrMsg,RoutineName)
      END IF
                  

      DO i=2,p_FAST%numIceLegs  ! basically, we just need IceDyn to set up its meshes for inputs/outputs and possibly initial values for states
         Init%InData_IceD%LegNum = i
         Init%InData_IceD%RootName = TRIM(p_FAST%OutFileRoot)//'.'//TRIM(y_FAST%Module_Abrev(Module_IceD))//TRIM(Num2LStr(i))     
         
         CALL IceD_Init( Init%InData_IceD, IceD%Input(1,i), IceD%p(i),  IceD%x(i,STATE_CURR), IceD%xd(i,STATE_CURR), IceD%z(i,STATE_CURR), &
                            IceD%OtherSt(i,STATE_CURR), IceD%y(i), IceD%m(i), dt_IceD, Init%OutData_IceD, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            
         !bjj: we're going to force this to have the same timestep because I don't want to have to deal with n IceD modules with n timesteps.
         IF (.NOT. EqualRealNos( p_FAST%dt_module( MODULE_IceD ),dt_IceD )) THEN
            CALL SetErrStat(ErrID_Fatal,"All instances of IceDyn (one per support-structure leg) must be the same",ErrStat,ErrMsg,RoutineName)
         END IF
      END DO
            
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF           
      
   END IF   
   

   ! ........................
   ! Set up output for glue code (must be done after all modules are initialized so we have their WriteOutput information)
   ! ........................

   CALL FAST_InitOutput( p_FAST, y_FAST, Init, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)


   ! -------------------------------------------------------------------------
   ! Initialize mesh-mapping data
   ! -------------------------------------------------------------------------

   CALL InitModuleMappings(p_FAST, ED, BD, AD14, AD, HD, SD, ExtPtfm, SrvD, MAPp, FEAM, MD, Orca, IceF, IceD, MeshMapData, ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      ELSEIF (ErrStat /= ErrID_None) THEN
         ! a little work-around in case the mesh mapping info messages get too long
         CALL WrScr( NewLine//TRIM(ErrMsg)//NewLine )
         ErrStat = ErrID_None
         ErrMsg = ""
      END IF
      
   ! -------------------------------------------------------------------------
   ! Initialize for linearization:
   ! -------------------------------------------------------------------------
   if ( p_FAST%Linearize ) then      
      call Init_Lin(p_FAST, y_FAST, m_FAST, AD, ED, Init%OutData_ED%NumBl, ErrStat2, ErrMsg2)      
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

         if (ErrStat >= AbortErrLev) then
            call Cleanup()
            return
         end if      
   end if
   
      
   ! -------------------------------------------------------------------------
   ! Initialize data for VTK output
   ! -------------------------------------------------------------------------
   if ( p_FAST%WrVTK > VTK_None ) then
      call SetVTKParameters(p_FAST, Init%OutData_ED, Init%OutData_AD, Init%InData_HD, Init%OutData_HD, ED, BD, AD, HD, ErrStat2, ErrMsg2)      
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   end if
   
   ! -------------------------------------------------------------------------
   ! Write initialization data to FAST summary file:
   ! -------------------------------------------------------------------------
   if (p_FAST%SumPrint)  then
       CALL FAST_WrSum( p_FAST, y_FAST, MeshMapData, ErrStat2, ErrMsg2 )
          CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   endif
   
   
   ! -------------------------------------------------------------------------
   ! other misc variables initialized here:
   ! -------------------------------------------------------------------------
      
   m_FAST%t_global        = t_initial
         
   ! Initialize external inputs for first step  
   if ( p_FAST%CompServo == MODULE_SrvD ) then      
      m_FAST%ExternInput%GenTrq     = SrvD%Input(1)%ExternalGenTrq !0.0_ReKi
      m_FAST%ExternInput%ElecPwr    = SrvD%Input(1)%ExternalElecPwr
      m_FAST%ExternInput%YawPosCom  = SrvD%Input(1)%ExternalYawPosCom
      m_FAST%ExternInput%YawRateCom = SrvD%Input(1)%ExternalYawRateCom
      m_FAST%ExternInput%HSSBrFrac  = SrvD%Input(1)%ExternalHSSBrFrac
      
      do i=1,SIZE(SrvD%Input(1)%ExternalBlPitchCom)
         m_FAST%ExternInput%BlPitchCom(i) = SrvD%Input(1)%ExternalBlPitchCom(i)
      end do   
   end if
   
   m_FAST%ExternInput%LidarFocus = 1.0_ReKi  ! make this non-zero (until we add the initial position in the InflowWind input file)
         
   
   !...............................................................................................................................
   ! Destroy initializion data
   !...............................................................................................................................      
   CALL Cleanup()
   
CONTAINS
   SUBROUTINE Cleanup()
   !...............................................................................................................................
   ! Destroy initializion data
   !...............................................................................................................................
      CALL FAST_DestroyInitData( Init, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   
   END SUBROUTINE Cleanup

END SUBROUTINE FAST_InitializeAll
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine checks for command-line arguments, gets the root name of the input files
!! (including full path name), and creates the names of the output files.
SUBROUTINE FAST_Init( p, m_FAST, y_FAST, t_initial, InputFile, ErrStat, ErrMsg, TMax, TurbID, OverrideAbortLev, RootName )

   IMPLICIT                        NONE

   ! Passed variables

   TYPE(FAST_ParameterType), INTENT(INOUT)         :: p                 !< The parameter data for the FAST (glue-code) simulation
   TYPE(FAST_MiscVarType),   INTENT(INOUT)         :: m_FAST            !< Miscellaneous variables
   TYPE(FAST_OutputFileType),INTENT(INOUT)         :: y_FAST            !< The output data for the FAST (glue-code) simulation
   REAL(DbKi),               INTENT(IN)            :: t_initial         !< the beginning time of the simulation
   INTEGER(IntKi),           INTENT(OUT)           :: ErrStat           !< Error status
   CHARACTER(*),             INTENT(OUT)           :: ErrMsg            !< Error message
   CHARACTER(*),             INTENT(IN)            :: InputFile         !< A CHARACTER string containing the name of the primary FAST input file (if not present, we'll get it from the command line)
   REAL(DbKi),               INTENT(IN), OPTIONAL  :: TMax              !< the length of the simulation (from Simulink or FAST.Farm)
   INTEGER(IntKi),           INTENT(IN), OPTIONAL  :: TurbID            !< an ID for naming the tubine output file
   LOGICAL,                  INTENT(IN), OPTIONAL  :: OverrideAbortLev  !< whether or not we should override the abort error level (e.g., FAST.Farm)
   CHARACTER(*),             INTENT(IN), OPTIONAL  :: RootName          !< A CHARACTER string containing the root name of FAST output files, overriding normal naming convention
      ! Local variables

   INTEGER                      :: i                                    ! loop counter
   !CHARACTER(1024)              :: DirName                              ! A CHARACTER string containing the path of the current working directory


   LOGICAL                      :: OverrideAbortErrLev  
   CHARACTER(*), PARAMETER      :: RoutineName = "FAST_Init"

   INTEGER(IntKi)               :: ErrStat2
   CHARACTER(ErrMsgLen)         :: ErrMsg2

      ! Initialize some variables
   ErrStat = ErrID_None
   ErrMsg = ''

   IF (PRESENT(OverrideAbortLev)) THEN
      OverrideAbortErrLev = OverrideAbortLev
   ELSE
      OverrideAbortErrLev = .true.
   END IF



   !...............................................................................................................................
   ! Set the root name of the output files based on the input file name
   !...............................................................................................................................

   if (present(RootName)) then
      p%OutFileRoot = RootName
   else         
         ! Determine the root name of the primary file (will be used for output files)
      CALL GetRoot( InputFile, p%OutFileRoot )
      IF ( Cmpl4SFun )  p%OutFileRoot = TRIM( p%OutFileRoot )//'.SFunc'
      IF ( PRESENT(TurbID) ) THEN
         IF ( TurbID > 0 ) THEN
            p%OutFileRoot = TRIM( p%OutFileRoot )//'.T'//TRIM(Num2LStr(TurbID))
         END IF
      END IF

   end if
   p%VTK_OutFileRoot = p%OutFileRoot !initialize this here in case of error before it is set later


   !...............................................................................................................................
   ! Initialize the module name/date/version info:
   !...............................................................................................................................

   DO i=1,NumModules
      y_FAST%Module_Ver(i)%Date = 'unknown date'
      y_FAST%Module_Ver(i)%Ver  = 'unknown version'
   END DO       
   y_FAST%Module_Ver( Module_IfW    )%Name = 'InflowWind'
   y_FAST%Module_Ver( Module_OpFM   )%Name = 'OpenFOAM integration'
   y_FAST%Module_Ver( Module_ED     )%Name = 'ElastoDyn'
   y_FAST%Module_Ver( Module_BD     )%Name = 'BeamDyn'
   y_FAST%Module_Ver( Module_AD14   )%Name = 'AeroDyn14'
   y_FAST%Module_Ver( Module_AD     )%Name = 'AeroDyn'
   y_FAST%Module_Ver( Module_SrvD   )%Name = 'ServoDyn'
   y_FAST%Module_Ver( Module_HD     )%Name = 'HydroDyn'
   y_FAST%Module_Ver( Module_SD     )%Name = 'SubDyn'
   y_FAST%Module_Ver( Module_ExtPtfm)%Name = 'ExtPtfm_MCKF'
   y_FAST%Module_Ver( Module_MAP    )%Name = 'MAP'
   y_FAST%Module_Ver( Module_FEAM   )%Name = 'FEAMooring'
   y_FAST%Module_Ver( Module_MD     )%Name = 'MoorDyn'
   y_FAST%Module_Ver( Module_Orca   )%Name = 'OrcaFlexInterface'
   y_FAST%Module_Ver( Module_IceF   )%Name = 'IceFloe'
   y_FAST%Module_Ver( Module_IceD   )%Name = 'IceDyn'
         
   y_FAST%Module_Abrev( Module_IfW    ) = 'IfW'
   y_FAST%Module_Abrev( Module_OpFM   ) = 'OpFM'
   y_FAST%Module_Abrev( Module_ED     ) = 'ED'
   y_FAST%Module_Abrev( Module_BD     ) = 'BD'
   y_FAST%Module_Abrev( Module_AD14   ) = 'AD'
   y_FAST%Module_Abrev( Module_AD     ) = 'AD'
   y_FAST%Module_Abrev( Module_SrvD   ) = 'SrvD'
   y_FAST%Module_Abrev( Module_HD     ) = 'HD'
   y_FAST%Module_Abrev( Module_SD     ) = 'SD'
   y_FAST%Module_Abrev( Module_ExtPtfm) = 'ExtPtfm'
   y_FAST%Module_Abrev( Module_MAP    ) = 'MAP'
   y_FAST%Module_Abrev( Module_FEAM   ) = 'FEAM'
   y_FAST%Module_Abrev( Module_MD     ) = 'MD'
   y_FAST%Module_Abrev( Module_Orca   ) = 'Orca'
   y_FAST%Module_Abrev( Module_IceF   ) = 'IceF'
   y_FAST%Module_Abrev( Module_IceD   ) = 'IceD'   

   p%n_substeps = 1                                                ! number of substeps for between modules and global/FAST time
   p%BD_OutputSibling = .false.

   !...............................................................................................................................
   ! Read the primary file for the glue code:
   !...............................................................................................................................
   CALL FAST_ReadPrimaryFile( InputFile, p, m_FAST, OverrideAbortErrLev, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 

      ! make sure some linearization variables are consistant
   if (.not. p%Linearize)  p%CalcSteady = .false.
   if (.not. p%CalcSteady) p%TrimCase = TrimCase_none
   m_FAST%Lin%FoundSteady = .false.
   p%LinInterpOrder = p%InterpOrder ! 1 ! always use linear (or constant) interpolation on rotor?

      ! overwrite TMax if necessary)
   IF (PRESENT(TMax)) THEN
      p%TMax = TMax
      !p%TMax = MAX( TMax, p%TMax )
   END IF

   IF ( ErrStat >= AbortErrLev ) RETURN


   p%KMax = 1                 ! after more checking, we may put this in the input file...
   !IF (p%CompIce == Module_IceF) p%KMax = 2
   p%SizeJac_Opt1 = 0  ! initialize this vector to zero; after we figure out what size the ED/SD/HD/BD meshes are, we'll fill this

   p%numIceLegs = 0           ! initialize number of support-structure legs in contact with ice (IceDyn will set this later)

   p%nBeams = 0               ! initialize number of BeamDyn instances (will be set later)

      ! determine what kind of turbine we're modeling:
   IF ( p%CompHydro == Module_HD ) THEN
      IF ( p%CompSub == Module_SD ) THEN
         p%TurbineType = Type_Offshore_Fixed
      ELSE
         p%TurbineType = Type_Offshore_Floating
      END IF
   ELSEIF ( p%CompMooring == Module_Orca ) THEN
      p%TurbineType = Type_Offshore_Floating
   ELSEIF ( p%CompSub == Module_ExtPtfm ) THEN
      p%TurbineType = Type_Offshore_Fixed
   ELSE      
      p%TurbineType = Type_LandBased
   END IF   
         
   
   p%n_TMax_m1  = CEILING( ( (p%TMax - t_initial) / p%DT ) ) - 1 ! We're going to go from step 0 to n_TMax (thus the -1 here)

   if (p%TMax < 1.0_DbKi) then ! log10(0) gives floating point divide-by-zero error
      p%TChanLen = MinChanLen
   else
      p%TChanLen = max( MinChanLen, int(log10(p%TMax))+7 )
   end if
   p%OutFmt_t = 'F'//trim(num2lstr( p%TChanLen ))//'.4' ! 'F10.4'    
   
   !...............................................................................................................................
   ! Do some error checking on the inputs (validation):
   !...............................................................................................................................   
   call ValidateInputData(p, m_FAST, ErrStat2, ErrMsg2)    
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
   


   IF ( ErrStat >= AbortErrLev ) RETURN


   RETURN
   END SUBROUTINE FAST_Init
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine initializes the output for the glue code, including writing the header for the primary output file.
SUBROUTINE FAST_InitOutput( p_FAST, y_FAST, Init, ErrStat, ErrMsg )

   IMPLICIT NONE

      ! Passed variables
   TYPE(FAST_ParameterType),       INTENT(IN)           :: p_FAST                                !< Glue-code simulation parameters
   TYPE(FAST_OutputFileType),      INTENT(INOUT)        :: y_FAST                                !< Glue-code simulation outputs
   TYPE(FAST_InitData),            INTENT(IN)           :: Init                                  !< Initialization data for all modules

   INTEGER(IntKi),                 INTENT(OUT)          :: ErrStat                               !< Error status
   CHARACTER(*),                   INTENT(OUT)          :: ErrMsg                                !< Error message corresponding to ErrStat


      ! Local variables.

   INTEGER(IntKi)                   :: I, J                                            ! Generic index for DO loops.
   INTEGER(IntKi)                   :: indxNext                                        ! The index of the next value to be written to an array
   INTEGER(IntKi)                   :: NumOuts                                         ! number of channels to be written to the output file(s)



   !......................................................
   ! Set the description lines to be printed in the output file
   !......................................................
   y_FAST%FileDescLines(1)  = 'Predictions were generated on '//CurDate()//' at '//CurTime()//' using '//TRIM(GetVersion(FAST_Ver))
   y_FAST%FileDescLines(2)  = 'linked with ' //' '//TRIM(GetNVD(NWTC_Ver            ))  ! we'll get the rest of the linked modules in the section below
   y_FAST%FileDescLines(3)  = 'Description from the FAST input file: '//TRIM(p_FAST%FTitle)
   
   !......................................................
   ! We'll fill out the rest of FileDescLines(2), 
   ! and save the module version info for later use, too:
   !......................................................

   y_FAST%Module_Ver( Module_ED ) = Init%OutData_ED%Ver
   y_FAST%FileDescLines(2) = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_ED )  ))

   IF ( p_FAST%CompElast == Module_BD )  THEN
      y_FAST%Module_Ver( Module_BD ) = Init%OutData_BD(1)%Ver ! call copy routine for this type if it every uses dynamic memory     
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_BD ))) 
   END IF   
   
   
   IF ( p_FAST%CompInflow == Module_IfW )  THEN
      y_FAST%Module_Ver( Module_IfW ) = Init%OutData_IfW%Ver ! call copy routine for this type if it every uses dynamic memory     
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_IfW ))) 
   ELSEIF ( p_FAST%CompInflow == Module_OpFM )  THEN
      y_FAST%Module_Ver( Module_OpFM ) = Init%OutData_OpFM%Ver ! call copy routine for this type if it every uses dynamic memory     
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_OpFM ))) 
   END IF   
   
   IF ( p_FAST%CompAero == Module_AD14 )  THEN
      y_FAST%Module_Ver( Module_AD14  ) = Init%OutData_AD14%Ver     
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_AD14  ) ))                  
   ELSEIF ( p_FAST%CompAero == Module_AD )  THEN
      y_FAST%Module_Ver( Module_AD  ) = Init%OutData_AD%Ver     
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_AD  ) ))                  
   END IF

   IF ( p_FAST%CompServo == Module_SrvD ) THEN
      y_FAST%Module_Ver( Module_SrvD ) = Init%OutData_SrvD%Ver
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_SrvD )))
   END IF
         
   IF ( p_FAST%CompHydro == Module_HD ) THEN
      y_FAST%Module_Ver( Module_HD )   = Init%OutData_HD%Ver
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_HD )))
   END IF

   IF ( p_FAST%CompSub == Module_SD ) THEN
      y_FAST%Module_Ver( Module_SD )   = Init%OutData_SD%Ver
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_SD )))
   ELSE IF ( p_FAST%CompSub == Module_ExtPtfm ) THEN
      y_FAST%Module_Ver( Module_ExtPtfm )   = Init%OutData_ExtPtfm%Ver
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_ExtPtfm )))
   END IF

   IF ( p_FAST%CompMooring == Module_MAP ) THEN
      y_FAST%Module_Ver( Module_MAP )   = Init%OutData_MAP%Ver
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_MAP )))
   ELSEIF ( p_FAST%CompMooring == Module_MD ) THEN
      y_FAST%Module_Ver( Module_MD )   = Init%OutData_MD%Ver
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_MD )))
   ELSEIF ( p_FAST%CompMooring == Module_FEAM ) THEN
      y_FAST%Module_Ver( Module_FEAM )   = Init%OutData_FEAM%Ver
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_FEAM )))
   ELSEIF ( p_FAST%CompMooring == Module_Orca ) THEN
      y_FAST%Module_Ver( Module_Orca )   = Init%OutData_Orca%Ver
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_Orca)))
   END IF   
   
   IF ( p_FAST%CompIce == Module_IceF ) THEN
      y_FAST%Module_Ver( Module_IceF )   = Init%OutData_IceF%Ver
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_IceF )))
   ELSEIF ( p_FAST%CompIce == Module_IceD ) THEN
      y_FAST%Module_Ver( Module_IceD )   = Init%OutData_IceD%Ver
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_IceD )))   
   END IF      
   
   !......................................................
   ! Set the number of output columns from each module
   !......................................................
   y_FAST%numOuts = 0    ! Inintialize entire array
   
   
   
   !y_FAST%numOuts(Module_InfW)  = 3  !hack for now: always output 3 wind speeds at hub-height
   IF ( ALLOCATED( Init%OutData_IfW%WriteOutputHdr  ) ) y_FAST%numOuts(Module_IfW)  = SIZE(Init%OutData_IfW%WriteOutputHdr)
   IF ( ALLOCATED( Init%OutData_OpFM%WriteOutputHdr ) ) y_FAST%numOuts(Module_OpFM) = SIZE(Init%OutData_OpFM%WriteOutputHdr)
   IF ( ALLOCATED( Init%OutData_ED%WriteOutputHdr   ) ) y_FAST%numOuts(Module_ED)   = SIZE(Init%OutData_ED%WriteOutputHdr)
do i=1,p_FAST%nBeams
   IF ( ALLOCATED( Init%OutData_BD(i)%WriteOutputHdr) ) y_FAST%numOuts(Module_BD)   = y_FAST%numOuts(Module_BD) + SIZE(Init%OutData_BD(i)%WriteOutputHdr)
end do   
!ad14 doesn't have outputs:
                                                         y_FAST%numOuts(Module_AD14) = 0
                                                         
   IF ( ALLOCATED( Init%OutData_AD%WriteOutputHdr     ) ) y_FAST%numOuts(Module_AD)     = SIZE(Init%OutData_AD%WriteOutputHdr)
   IF ( ALLOCATED( Init%OutData_SrvD%WriteOutputHdr   ) ) y_FAST%numOuts(Module_SrvD)   = SIZE(Init%OutData_SrvD%WriteOutputHdr)
   IF ( ALLOCATED( Init%OutData_HD%WriteOutputHdr     ) ) y_FAST%numOuts(Module_HD)     = SIZE(Init%OutData_HD%WriteOutputHdr)
   IF ( ALLOCATED( Init%OutData_SD%WriteOutputHdr     ) ) y_FAST%numOuts(Module_SD)     = SIZE(Init%OutData_SD%WriteOutputHdr)
   IF ( ALLOCATED( Init%OutData_ExtPtfm%WriteOutputHdr) ) y_FAST%numOuts(Module_ExtPtfm)= SIZE(Init%OutData_ExtPtfm%WriteOutputHdr)
   IF ( ALLOCATED( Init%OutData_MAP%WriteOutputHdr    ) ) y_FAST%numOuts(Module_MAP)    = SIZE(Init%OutData_MAP%WriteOutputHdr)
   IF ( ALLOCATED( Init%OutData_FEAM%WriteOutputHdr   ) ) y_FAST%numOuts(Module_FEAM)   = SIZE(Init%OutData_FEAM%WriteOutputHdr)
   IF ( ALLOCATED( Init%OutData_MD%WriteOutputHdr     ) ) y_FAST%numOuts(Module_MD)     = SIZE(Init%OutData_MD%WriteOutputHdr)
   IF ( ALLOCATED( Init%OutData_Orca%WriteOutputHdr   ) ) y_FAST%numOuts(Module_Orca)   = SIZE(Init%OutData_Orca%WriteOutputHdr)
   IF ( ALLOCATED( Init%OutData_IceF%WriteOutputHdr   ) ) y_FAST%numOuts(Module_IceF)   = SIZE(Init%OutData_IceF%WriteOutputHdr)
   IF ( ALLOCATED( Init%OutData_IceD%WriteOutputHdr   ) ) y_FAST%numOuts(Module_IceD)   = SIZE(Init%OutData_IceD%WriteOutputHdr)*p_FAST%numIceLegs         
   
   !......................................................
   ! Initialize the output channel names and units
   !......................................................
   NumOuts   = 1 + SUM( y_FAST%numOuts )

   CALL AllocAry( y_FAST%ChannelNames,NumOuts, 'ChannelNames', ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( y_FAST%ChannelUnits,NumOuts, 'ChannelUnits', ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN

   y_FAST%ChannelNames(1) = 'Time'
   y_FAST%ChannelUnits(1) = '(s)'

   indxNext = 2
   DO i=1,y_FAST%numOuts(Module_IfW) !InflowWind
      y_FAST%ChannelNames(indxNext) = Init%OutData_IfW%WriteOutputHdr(i)
      y_FAST%ChannelUnits(indxNext) = Init%OutData_IfW%WriteOutputUnt(i)
      indxNext = indxNext + 1
   END DO

   DO i=1,y_FAST%numOuts(Module_OpFM) !OpenFOAM
      y_FAST%ChannelNames(indxNext) = Init%OutData_OpFM%WriteOutputHdr(i)
      y_FAST%ChannelUnits(indxNext) = Init%OutData_OpFM%WriteOutputUnt(i)
      indxNext = indxNext + 1
   END DO   

   DO i=1,y_FAST%numOuts(Module_ED) !ElastoDyn
      y_FAST%ChannelNames(indxNext) = Init%OutData_ED%WriteOutputHdr(i)
      y_FAST%ChannelUnits(indxNext) = Init%OutData_ED%WriteOutputUnt(i)
      indxNext = indxNext + 1
   END DO   

   IF ( y_FAST%numOuts(Module_BD) > 0_IntKi ) THEN !BeamDyn
      do i=1,p_FAST%nBeams
         if ( allocated(Init%OutData_BD(i)%WriteOutputHdr) ) then            
            do j=1,size(Init%OutData_BD(i)%WriteOutputHdr) 
               y_FAST%ChannelNames(indxNext) = 'B'//TRIM(Num2Lstr(i))//trim(Init%OutData_BD(i)%WriteOutputHdr(j))
               y_FAST%ChannelUnits(indxNext) = Init%OutData_BD(i)%WriteOutputUnt(j)
               indxNext = indxNext + 1
            end do ! j            
         end if         
      end do                 
   END IF
   
   
   ! none for AeroDyn14 
   
   DO i=1,y_FAST%numOuts(Module_AD) !AeroDyn
      y_FAST%ChannelNames(indxNext) = Init%OutData_AD%WriteOutputHdr(i)
      y_FAST%ChannelUnits(indxNext) = Init%OutData_AD%WriteOutputUnt(i)
      indxNext = indxNext + 1
   END DO   
   
   DO i=1,y_FAST%numOuts(Module_SrvD) !ServoDyn
      y_FAST%ChannelNames(indxNext) = Init%OutData_SrvD%WriteOutputHdr(i)
      y_FAST%ChannelUnits(indxNext) = Init%OutData_SrvD%WriteOutputUnt(i)
      indxNext = indxNext + 1
   END DO   

   DO i=1,y_FAST%numOuts(Module_HD) !HydroDyn
      y_FAST%ChannelNames(indxNext) = Init%OutData_HD%WriteOutputHdr(i)
      y_FAST%ChannelUnits(indxNext) = Init%OutData_HD%WriteOutputUnt(i)
      indxNext = indxNext + 1
   END DO      

   DO i=1,y_FAST%numOuts(Module_SD) !SubDyn
      y_FAST%ChannelNames(indxNext) = Init%OutData_SD%WriteOutputHdr(i)
      y_FAST%ChannelUnits(indxNext) = Init%OutData_SD%WriteOutputUnt(i)
      indxNext = indxNext + 1
   END DO      

   DO i=1,y_FAST%numOuts(Module_ExtPtfm) !ExtPtfm_MCKF
      y_FAST%ChannelNames(indxNext) = Init%OutData_ExtPtfm%WriteOutputHdr(i)
      y_FAST%ChannelUnits(indxNext) = Init%OutData_ExtPtfm%WriteOutputUnt(i)
      indxNext = indxNext + 1
   END DO      
   
   DO i=1,y_FAST%numOuts(Module_MAP) !MAP
      y_FAST%ChannelNames(indxNext) = Init%OutData_MAP%WriteOutputHdr(i)
      y_FAST%ChannelUnits(indxNext) = Init%OutData_MAP%WriteOutputUnt(i)
      indxNext = indxNext + 1
   END DO      
   
   DO i=1,y_FAST%numOuts(Module_MD) !MoorDyn
      y_FAST%ChannelNames(indxNext) = Init%OutData_MD%WriteOutputHdr(i)
      y_FAST%ChannelUnits(indxNext) = Init%OutData_MD%WriteOutputUnt(i)
      indxNext = indxNext + 1
   END DO      

   DO i=1,y_FAST%numOuts(Module_FEAM) !FEAMooring
      y_FAST%ChannelNames(indxNext) = Init%OutData_FEAM%WriteOutputHdr(i)
      y_FAST%ChannelUnits(indxNext) = Init%OutData_FEAM%WriteOutputUnt(i)
      indxNext = indxNext + 1
   END DO      

   DO i=1,y_FAST%numOuts(Module_Orca) !OrcaFlex
      y_FAST%ChannelNames(indxNext) = Init%OutData_Orca%WriteOutputHdr(i)
      y_FAST%ChannelUnits(indxNext) = Init%OutData_Orca%WriteOutputUnt(i)
      indxNext = indxNext + 1
   END DO         
   
   DO i=1,y_FAST%numOuts(Module_IceF) !IceFloe
      y_FAST%ChannelNames(indxNext) = Init%OutData_IceF%WriteOutputHdr(i)
      y_FAST%ChannelUnits(indxNext) = Init%OutData_IceF%WriteOutputUnt(i)
      indxNext = indxNext + 1
   END DO         
   
   IF ( y_FAST%numOuts(Module_IceD) > 0_IntKi ) THEN !IceDyn
      DO I=1,p_FAST%numIceLegs         
         DO J=1,SIZE(Init%OutData_IceD%WriteOutputHdr) 
            y_FAST%ChannelNames(indxNext) =TRIM(Init%OutData_IceD%WriteOutputHdr(J))//'L'//TRIM(Num2Lstr(I))  !bjj: do we want this "Lx" at the end?
            y_FAST%ChannelUnits(indxNext) = Init%OutData_IceD%WriteOutputUnt(J)
            indxNext = indxNext + 1
         END DO ! J
      END DO ! I
   END IF   
      
   
   !......................................................
   ! Open the text output file and print the headers
   !......................................................

   IF (p_FAST%WrTxtOutFile) THEN

      y_FAST%ActualChanLen = max( MinChanLen, p_FAST%FmtWidth )
      DO I=1,NumOuts
         y_FAST%ActualChanLen = max( y_FAST%ActualChanLen, LEN_TRIM(y_FAST%ChannelNames(I)) )
         y_FAST%ActualChanLen = max( y_FAST%ActualChanLen, LEN_TRIM(y_FAST%ChannelUnits(I)) )
      ENDDO ! I

      y_FAST%OutFmt_a = '"'//p_FAST%Delim//'"'//p_FAST%OutFmt      ! format for array elements from individual modules
      if (p_FAST%FmtWidth < y_FAST%ActualChanLen) then
         y_FAST%OutFmt_a = trim(y_FAST%OutFmt_a)//','//trim(num2lstr(y_FAST%ActualChanLen - p_FAST%FmtWidth))//'x'
      end if
      
      CALL GetNewUnit( y_FAST%UnOu, ErrStat, ErrMsg )
         IF ( ErrStat >= AbortErrLev ) RETURN

      CALL OpenFOutFile ( y_FAST%UnOu, TRIM(p_FAST%OutFileRoot)//'.out', ErrStat, ErrMsg )
         IF ( ErrStat >= AbortErrLev ) RETURN

         ! Add some file information:

      WRITE (y_FAST%UnOu,'(/,A)')  TRIM( y_FAST%FileDescLines(1) )
      WRITE (y_FAST%UnOu,'(1X,A)') TRIM( y_FAST%FileDescLines(2) )
      WRITE (y_FAST%UnOu,'()' )    !print a blank line
      WRITE (y_FAST%UnOu,'(A)'   ) TRIM( y_FAST%FileDescLines(3) )
      WRITE (y_FAST%UnOu,'()' )    !print a blank line


         !......................................................
         ! Write the names of the output parameters on one line:
         !......................................................
      if (p_FAST%Delim /= " ") then ! trim trailing spaces if not space delimited:

         CALL WrFileNR ( y_FAST%UnOu, trim(y_FAST%ChannelNames(1)) ) ! first one is time, with a special format

         DO I=2,NumOuts
            CALL WrFileNR ( y_FAST%UnOu, p_FAST%Delim//trim(y_FAST%ChannelNames(I)) )
         ENDDO ! I
      else
      
         CALL WrFileNR ( y_FAST%UnOu, y_FAST%ChannelNames(1)(1:p_FAST%TChanLen) ) ! first one is time, with a special format

         DO I=2,NumOuts
            CALL WrFileNR ( y_FAST%UnOu, p_FAST%Delim//y_FAST%ChannelNames(I)(1:y_FAST%ActualChanLen) )
         ENDDO ! I
      end if

      WRITE (y_FAST%UnOu,'()')

         !......................................................
         ! Write the units of the output parameters on one line:
         !......................................................

      if (p_FAST%Delim /= " ") then
      
         CALL WrFileNR ( y_FAST%UnOu, trim(y_FAST%ChannelUnits(1)) )

         DO I=2,NumOuts
            CALL WrFileNR ( y_FAST%UnOu, p_FAST%Delim//trim(y_FAST%ChannelUnits(I)) )
         ENDDO ! I
      else
      
         CALL WrFileNR ( y_FAST%UnOu, y_FAST%ChannelUnits(1)(1:p_FAST%TChanLen) )

         DO I=2,NumOuts
            CALL WrFileNR ( y_FAST%UnOu, p_FAST%Delim//y_FAST%ChannelUnits(I)(1:y_FAST%ActualChanLen) )
         ENDDO ! I
      end if

      WRITE (y_FAST%UnOu,'()')

   END IF

   !......................................................
   ! Allocate data for binary output file
   !......................................................
   IF (p_FAST%WrBinOutFile) THEN

         ! calculate the size of the array of outputs we need to store
      y_FAST%NOutSteps = CEILING ( (p_FAST%TMax - p_FAST%TStart) / p_FAST%DT_OUT ) + 1

      CALL AllocAry( y_FAST%AllOutData, NumOuts-1, y_FAST%NOutSteps, 'AllOutData', ErrStat, ErrMsg )
      y_FAST%AllOutData = 0.0_ReKi
      IF ( ErrStat >= AbortErrLev ) RETURN

      IF ( p_FAST%WrBinMod == FileFmtID_WithTime ) THEN   ! we store the entire time array
         CALL AllocAry( y_FAST%TimeData, y_FAST%NOutSteps, 'TimeData', ErrStat, ErrMsg )
         IF ( ErrStat >= AbortErrLev ) RETURN
      ELSE  
         CALL AllocAry( y_FAST%TimeData, 2_IntKi, 'TimeData', ErrStat, ErrMsg )
         IF ( ErrStat >= AbortErrLev ) RETURN

         y_FAST%TimeData(1) = 0.0_DbKi           ! This is the first output time, which we will set later
         y_FAST%TimeData(2) = p_FAST%DT_out      ! This is the (constant) time between subsequent writes to the output file
      END IF

      y_FAST%n_Out = 0  !number of steps actually written to the file

   END IF

   y_FAST%VTK_count = 0  ! first VTK file has 0 as output
   
RETURN
END SUBROUTINE FAST_InitOutput
!----------------------------------------------------------------------------------------------------------------------------------
!> This function returns a string describing the glue code and some of the compilation options we're using.
FUNCTION GetVersion(ThisProgVer)

   ! Passed Variables:

   TYPE(ProgDesc), INTENT( IN    ) :: ThisProgVer     !< program name/date/version description
   CHARACTER(1024)                 :: GetVersion      !< String containing a description of the compiled precision.
   
   CHARACTER(200)                  :: git_commit
   
   GetVersion = TRIM(GetNVD(ThisProgVer))//', compiled'

   IF ( Cmpl4SFun )  THEN     ! FAST has been compiled as an S-Function for Simulink
      GetVersion = TRIM(GetVersion)//' as a DLL S-Function for Simulink'
   ELSEIF ( Cmpl4LV )  THEN     ! FAST has been compiled as a DLL for Labview
      GetVersion = TRIM(GetVersion)//' as a DLL for LabVIEW'
   ENDIF   
   
   GetVersion = TRIM(GetVersion)//' as a '//TRIM(Num2LStr(BITS_IN_ADDR))//'-bit application using'
   
   ! determine precision

      IF ( ReKi == SiKi )  THEN     ! Single precision
         GetVersion = TRIM(GetVersion)//' single'
      ELSEIF ( ReKi == R8Ki )  THEN ! Double precision
         GetVersion = TRIM(GetVersion)// ' double'
      ELSE                          ! Unknown precision
         GetVersion = TRIM(GetVersion)//' unknown'
      ENDIF
      

!   GetVersion = TRIM(GetVersion)//' precision with '//OS_Desc
   GetVersion = TRIM(GetVersion)//' precision'

   ! add git info
   git_commit = QueryGitVersion()
   GetVersion = TRIM(GetVersion)//' at commit '//git_commit

   RETURN
END FUNCTION GetVersion
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine sets up the information needed to initialize AeroDyn, then initializes AeroDyn
SUBROUTINE AD_SetInitInput(InitInData_AD14, InitOutData_ED, y_ED, p_FAST, ErrStat, ErrMsg)

   ! Passed variables:
   TYPE(AD14_InitInputType),INTENT(INOUT) :: InitInData_AD14  !< The initialization input to AeroDyn14
   TYPE(ED_InitOutputType), INTENT(IN)    :: InitOutData_ED   !< The initialization output from structural dynamics module
   TYPE(ED_OutputType),     INTENT(IN)    :: y_ED             !< The outputs of the structural dynamics module (meshes with position/RefOrientation set)
   TYPE(FAST_ParameterType),INTENT(IN)    :: p_FAST           !< The parameters of the glue code
   INTEGER(IntKi)                         :: ErrStat          !< Error status of the operation
   CHARACTER(*)                           :: ErrMsg           !< Error message if ErrStat /= ErrID_None

      ! Local variables

   !TYPE(AD_InitOptions)       :: ADOptions                  ! Options for AeroDyn

   INTEGER                    :: K


   ErrStat = ErrID_None
   ErrMsg  = ""
   
   
      ! Set up the AeroDyn parameters
   InitInData_AD14%ADFileName   = p_FAST%AeroFile
   InitInData_AD14%OutRootName  = p_FAST%OutFileRoot
   InitInData_AD14%WrSumFile    = p_FAST%SumPrint      
   InitInData_AD14%NumBl        = InitOutData_ED%NumBl
   InitInData_AD14%UseDWM       = p_FAST%UseDWM
   
   InitInData_AD14%DWM%IfW%InputFileName   = p_FAST%InflowFile
   
      ! Hub position and orientation (relative here, but does not need to be)

   InitInData_AD14%TurbineComponents%Hub%Position(:)      = y_ED%HubPtMotion14%Position(:,1) - y_ED%HubPtMotion14%Position(:,1)  ! bjj: was 0; mesh was changed by adding p_ED%HubHt to 3rd component
   InitInData_AD14%TurbineComponents%Hub%Orientation(:,:) = y_ED%HubPtMotion14%RefOrientation(:,:,1)
   InitInData_AD14%TurbineComponents%Hub%TranslationVel   = 0.0_ReKi ! bjj: we don't need this field
   InitInData_AD14%TurbineComponents%Hub%RotationVel      = 0.0_ReKi ! bjj: we don't need this field

      ! Blade root position and orientation (relative here, but does not need to be)

   IF (.NOT. ALLOCATED( InitInData_AD14%TurbineComponents%Blade ) ) THEN
      ALLOCATE( InitInData_AD14%TurbineComponents%Blade( InitInData_AD14%NumBl ), STAT = ErrStat )
      IF ( ErrStat /= 0 ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg = ' Error allocating space for InitInData_AD%TurbineComponents%Blade.'
         RETURN
      ELSE
         ErrStat = ErrID_None !reset to ErrID_None, just in case ErrID_None /= 0
      END IF
   END IF

   DO K=1, InitInData_AD14%NumBl
      InitInData_AD14%TurbineComponents%Blade(K)%Position        = y_ED%BladeRootMotion14%Position(:,K)
      InitInData_AD14%TurbineComponents%Blade(K)%Orientation     = y_ED%BladeRootMotion14%RefOrientation(:,:,K)
      InitInData_AD14%TurbineComponents%Blade(K)%TranslationVel  = 0.0_ReKi ! bjj: we don't need this field
      InitInData_AD14%TurbineComponents%Blade(K)%RotationVel     = 0.0_ReKi ! bjj: we don't need this field      
   END DO
  

      ! Blade length
   IF (p_FAST%CompElast == Module_ED) THEN  ! note, we can't get here if we're using BeamDyn....
      InitInData_AD14%TurbineComponents%BladeLength = InitOutData_ED%BladeLength
   END IF
   
   
      ! Tower mesh ( here only because we currently need line2 meshes to contain the same nodes/elements )
      
   InitInData_AD14%NumTwrNodes = y_ED%TowerLn2Mesh%NNodes - 2
   IF (.NOT. ALLOCATED( InitInData_AD14%TwrNodeLocs ) ) THEN
      ALLOCATE( InitInData_AD14%TwrNodeLocs( 3, InitInData_AD14%NumTwrNodes ), STAT = ErrStat )
      IF ( ErrStat /= 0 ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg = ' Error allocating space for InitInData_AD%TwrNodeLocs.'
         RETURN
      ELSE
         ErrStat = ErrID_None
      END IF
   END IF   
   
   IF ( InitInData_AD14%NumTwrNodes > 0 ) THEN
      InitInData_AD14%TwrNodeLocs = y_ED%TowerLn2Mesh%Position(:,1:InitInData_AD14%NumTwrNodes)  ! ED has extra nodes at beginning and top and bottom of tower
   END IF
   
      ! hub height         
   InitInData_AD14%HubHt = InitOutData_ED%HubHt
             

   RETURN
END SUBROUTINE AD_SetInitInput
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine is called at the start (or restart) of a FAST program (or FAST.Farm). It initializes the NWTC subroutine library,
!! displays the copyright notice, and displays some version information (including addressing scheme and precision).
SUBROUTINE FAST_ProgStart(ThisProgVer)
   TYPE(ProgDesc), INTENT(IN) :: ThisProgVer     !< program name/date/version description
   
   ! ... Initialize NWTC Library (open console, set pi constants) ...
   ! sets the pi constants, open console for output, etc...
   CALL NWTC_Init( ProgNameIN=ThisProgVer%Name, EchoLibVer=.FALSE. )
   
   ! Display the copyright notice
   CALL DispCopyrightLicense( ThisProgVer%Name )
   
   CALL DispCompileRuntimeInfo

END SUBROUTINE FAST_ProgStart
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine sets the number of subcycles (substeps) for modules at initialization, checking to make sure that their requested 
!! time step is valid.
SUBROUTINE SetModuleSubstepTime(ModuleID, p_FAST, y_FAST, ErrStat, ErrMsg)
   INTEGER(IntKi),           INTENT(IN   ) :: ModuleID            !< ID of the module to check time step and set
   TYPE(FAST_ParameterType), INTENT(INOUT) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(IN   ) :: y_FAST              !< Output variables for the glue code
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

      
   ErrStat = ErrID_None
   ErrMsg  = "" 
   
   IF ( EqualRealNos( p_FAST%dt_module( ModuleID ), p_FAST%dt ) ) THEN
      p_FAST%n_substeps(ModuleID) = 1
   ELSE
      IF ( p_FAST%dt_module( ModuleID ) > p_FAST%dt ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg = "The "//TRIM(y_FAST%Module_Ver(ModuleID)%Name)//" module time step ("//&
                          TRIM(Num2LStr(p_FAST%dt_module( ModuleID )))// &
                    " s) cannot be larger than FAST time step ("//TRIM(Num2LStr(p_FAST%dt))//" s)."
      ELSE
            ! calculate the number of subcycles:
         p_FAST%n_substeps(ModuleID) = NINT( p_FAST%dt / p_FAST%dt_module( ModuleID ) )
            
            ! let's make sure THE module DT is an exact integer divisor of the global (FAST) time step:
         IF ( .NOT. EqualRealNos( p_FAST%dt, p_FAST%dt_module( ModuleID ) * p_FAST%n_substeps(ModuleID) )  ) THEN
            ErrStat = ErrID_Fatal
            ErrMsg  = "The "//TRIM(y_FAST%Module_Ver(ModuleID)%Name)//" module time step ("//&
                              TRIM(Num2LStr(p_FAST%dt_module( ModuleID )))// &
                              " s) must be an integer divisor of the FAST time step ("//TRIM(Num2LStr(p_FAST%dt))//" s)."
         END IF
            
      END IF
   END IF      
                 
   RETURN
      
END SUBROUTINE SetModuleSubstepTime   

END MODULE
