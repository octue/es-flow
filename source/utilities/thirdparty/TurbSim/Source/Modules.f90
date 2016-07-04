!=======================================================================
MODULE TSMods

USE                             NWTC_Library

IMPLICIT                        NONE
SAVE

REAL(ReKi)                   :: AnalysisTime                             ! Analysis Time. (amount of time for analysis, allows user to perform analysis using one time length, but output UsableTime
REAL(ReKi)                   :: ChebyCoef_WS(11)                         ! The Chebyshev coefficients for wind speed
REAL(ReKi)                   :: ChebyCoef_WD(11)                         ! The Chebyshev coefficients for wind direction
REAL(ReKi)                   :: COHEXP                                   ! Coherence exponent
REAL(ReKi)                   :: CTLy                                     ! Fractional location of tower centerline from right (looking downwind) to left side of the dataset.
REAL(ReKi)                   :: CTLz                                     ! Fractional location of hub height from the bottom of the dataset.
REAL(ReKi)                   :: CTStartTime                              ! Minimum time to add coherent structures
REAL(ReKi)                   :: DistScl                                  ! Disturbance scale for AeroDyn coherent turbulence events
REAL(ReKi)                   :: ETMc                                     ! The c parameter in IEC ETM, 61400-1, Ed 3. Section 6.3.2.3, Eq. 19.  Variable per last sentence in section 7.4.1
REAL(ReKi), ALLOCATABLE      :: EventLen   (:)                           ! The length of each event stored in EventStart() (non-dimensional time)
REAL(ReKi)                   :: EventTimeStep                            ! The average length of timesteps in output events
REAL(ReKi)                   :: Fc                                       ! Coriolis parameter in units (1/sec)
REAL(ReKi), ALLOCATABLE      :: Freq       (:)                           ! The array of frequencies (NumFreq).
REAL(ReKi), ALLOCATABLE      :: Freq_USR(:)                              ! frequencies for the user-defined spectra
REAL(ReKi)                   :: GridHeight                               ! Grid height
REAL(ReKi)                   :: GridRes_Y                                ! Distance between two consecutive horizontal points on the grid (Horizontal resolution)
REAL(ReKi)                   :: GridRes_Z                                ! Distance between two consecutive vertical points on the grid (Vertical resolution)
REAL(ReKi)                   :: GridWidth                                ! Grid width.
REAL(ReKi)                   :: h                                        ! Boundary layer depth
REAL(ReKi)                   :: HFlowAng                                 ! Horizontal flow angle.
REAL(ReKi)                   :: HubHt                                    ! Hub height.
REAL(ReKi)                   :: InCDec     (3)                           ! Contains the coherence decrements
REAL(ReKi)                   :: InCohB     (3)                           ! Contains the coherence b/L (offset) parameters
REAL(ReKi)                   :: L                                        ! M-O length
REAL(ReKi), ALLOCATABLE      :: L_USR      (:)                           ! User-specified von Karman length scale, varying with height
REAL(ReKi)                   :: Latitude                                 ! The site latitude in radians
REAL(ReKi), PARAMETER        :: Omega     = 7.292116E-05                 ! Angular speed of rotation of the earth (rad/s)
REAL(ReKi)                   :: PerTurbInt                               ! Percent Turbulence Intensity
REAL(ReKi)                   :: PC_UW                                    ! u'w' cross-correlation coefficient
REAL(ReKi)                   :: PC_UV                                    ! u'v' cross-correlation coefficient
REAL(ReKi)                   :: PC_VW                                    ! v'w' cross-correlation coefficient
REAL(ReKi), ALLOCATABLE      :: pkCTKE     (:)                           ! Array containing the peak CTKE of each coherent event
REAL(ReKi)                   :: PLExp                                    ! Rotor disk power law exponent
REAL(ReKi), PARAMETER        :: profileZmax = 140.                       ! Upper height limit for extrapolating GP_LLJ profiles of ustar and zl
REAL(ReKi), PARAMETER        :: profileZmin =  50.                       ! Lower height limit for extrapolating GP_LLJ profiles of ustar and zl
REAL(ReKi), ALLOCATABLE      :: RandNum    (:)                           ! The array that holds the random numbers.
REAL(ReKi)                   :: RICH_NO                                  ! Gradient Richardson number
REAL(ReKi)                   :: RotorDiameter                            ! The assumed diameter of the rotor
REAL(ReKi), ALLOCATABLE      :: S          (:,:,:)                       ! The turbulence PSD array (NumFreq,NTot,3).
REAL(ReKi), ALLOCATABLE      :: SDary      (:)                           ! The array of standard deviations (NumGrid_Z,NumGrid_Y).
REAL(ReKi)                   :: SigmaIEC                                 ! IEC standard deviation.
REAL(ReKi), ALLOCATABLE      :: Sigma_USR  (:)                           ! User-specified standard deviation of the wind speed components (isotropic), varying with height
REAL(ReKi)                   :: StdScale   (3)                           ! Scaling for the user-specified standard deviation
REAL(ReKi)                   :: TimeStep                                 ! Time step.
REAL(ReKi)                   :: Sigma_U2                                 ! Standard Deviation of U velocity, squared.
REAL(ReKi)                   :: Sigma_V2                                 ! Standard Deviation of V velocity, squared.
REAL(ReKi)                   :: Sigma_W2                                 ! Standard Deviation of W velocity, squared.
REAL(ReKi)                   :: TurbIntH20                               ! Turbulence intensity used for HYDRO module.
REAL(ReKi), PARAMETER        :: Tolerance = 0.0001                       ! The largest difference between two numbers that are assumed to be equal
REAL(ReKi), ALLOCATABLE      :: TRH        (:)                           ! The transfer function  matrix (NumSteps).
REAL(ReKi)                   :: TsclFact                                 ! Scale factor for time (h/U0) in coherent turbulence events
REAL(ReKi), ALLOCATABLE      :: U          (:)                           ! The steady u-component wind speeds for the grid (ZLim).
REAL(ReKi)                   :: H_ref                                    ! Height for reference wind speed.
REAL(ReKi), ALLOCATABLE      :: DUDZ       (:)                           ! The steady u-component wind shear for the grid (ZLim).
REAL(ReKi), ALLOCATABLE      :: U_USR      (:)                           ! User-specified total wind speed, varying with height
REAL(ReKi)                   :: UHub                                     ! Hub-height (total) wind speed (m/s)
REAL(ReKi)                   :: UJetMax                                  ! The (horizontal) wind speed at the height of the jet maximum (m/s)
REAL(ReKi)                   :: UsableTime                               ! Usable time.  Program adds GridWidth/MeanHHWS.
REAL(ReKi), ALLOCATABLE      :: Uspec_USR(:)                             ! user-defined u-component spectrum
REAL(ReKi)                   :: Ustar                                    ! Shear or friction velocity (m/s) -- rotor-disk average
REAL(ReKi), ALLOCATABLE      :: Ustar_profile(:)                         ! A profile of ustar (measure of friction velocity with height)
REAL(ReKi)                   :: UstarDiab                                ! The diabatic ustar value
REAL(ReKi)                   :: UstarOffset                              ! A scaling/offset value used with the Ustar_profile to ensure that the mean hub u'w' and ustar inputs agree with the profile values
REAL(ReKi)                   :: UstarSlope                               ! A scaling/slope value used with the Ustar_profile to ensure that the mean hub u'w' and ustar inputs agree with the profile values
REAL(ReKi), ALLOCATABLE      :: V          (:,:,:)                       ! An array containing the summations of the rows of H (NumSteps,NTot,3).
REAL(ReKi)                   :: VFlowAng                                 ! Vertical flow angle.
REAL(ReKi)                   :: Vave                                     ! The IEC Vave for ETM
REAL(ReKi)                   :: Vref                                     ! The IEC Vref for ETM
REAL(ReKi), ALLOCATABLE      :: Vspec_USR(:)                             ! user-defined v-component spectrum

REAL(ReKi), ALLOCATABLE      :: WindDir_profile(:)                       ! A profile of horizontal wind angle (measure of wind direction with height)
REAL(ReKi), ALLOCATABLE      :: WindDir_USR    (:)                       ! User-specified wind direction profile, varying with height
REAL(ReKi), ALLOCATABLE      :: Work       (:,:)                         ! A temporary work array (NumSteps+2,3).
REAL(ReKi), ALLOCATABLE      :: Wspec_USR(:)                             ! user-defined w-component spectrum
REAL(ReKi), ALLOCATABLE      :: Y          (:)                           ! The lateral locations of the points (YLim).
REAL(ReKi)                   :: Ym_max                                   ! The nondimensional lateral width of the coherent turbulence dataset
REAL(ReKi), ALLOCATABLE      :: Z          (:)                           ! The vertical locations of the points (ZLim).
REAL(ReKi), ALLOCATABLE      :: Z_USR      (:)                           ! Heights of user-specified variables
REAL(ReKi)                   :: Z0                                       ! Surface roughness length, meters
REAL(ReKi)                   :: ZI                                       ! Mixing layer depth
REAL(ReKi)                   :: ZJetMax                                  ! The height of the jet maximum (m)
REAL(ReKi)                   :: ZL                                       ! A measure of stability
REAL(ReKi), ALLOCATABLE      :: ZL_profile(:)                            ! A profile of z/l (measure of stability with height)
REAL(ReKi)                   :: ZLoffset                                 ! An offset to align the zl profile with the mean zl input parameter
REAL(ReKi)                   :: Zm_max                                   ! The nondimensional vertical height of the coherent turbulence dataset



INTEGER,    ALLOCATABLE      :: EventName  (:)                           ! The timestep where the event starts, which determines the name of the event file
INTEGER,    ALLOCATABLE      :: EventTS    (:)                           ! The length of each event stored in EventStart() (number of timesteps)
INTEGER                      :: IECedition                               ! The edition number of the IEC 61400-1 standard that is being used (determines the scaling)
INTEGER                      :: IECstandard                              ! The standard number (x) of the IEC 61400-x that is being used
INTEGER,    PARAMETER        :: IEC_ETM        = 1                       ! Number to indicate the IEC Normal Turbulence Model
INTEGER,    PARAMETER        :: IEC_EWM1       = 2                       ! Number to indicate the IEC Extreme Wind speed Model (50-year)
INTEGER,    PARAMETER        :: IEC_EWM50      = 3                       ! Number to indicate the IEC Extreme Wind speed Model ( 1-year)
INTEGER,    PARAMETER        :: IEC_NTM        = 4                       ! Number to indicate the IEC Extreme Turbulence Model
INTEGER                      :: IEC_WindType                             ! Number to indicate the IEC wind type
INTEGER,    ALLOCATABLE      :: IYmax      (:)                           ! A temporary variable holding the maximum number of horizontal positions at each z
INTEGER                      :: LuxLevel   = 3                           ! Luxury Level for RanLux RNG
INTEGER                      :: MaxDims                                  ! Maximum number of time steps plus 2.
INTEGER                      :: NTot                                     ! Number of points in grid, plus the hub center.
INTEGER                      :: NumCTEvents                              ! Number of events to be inserted into the .cts file
INTEGER                      :: NumCTt                                   ! Number of data points in the output coherent event timestep file
INTEGER                      :: NumEvents                                ! Number of events in the event data file
INTEGER                      :: NumFreq                                  ! Number of frequencies (=NumSteps/2).
INTEGER                      :: NumGrid_Y                                ! Grid dimension. (in horizontal direction)
INTEGER                      :: NumGrid_Z                                ! Grid dimension. (in vertical direction)
INTEGER                      :: NumOutSteps                              ! Number of output time steps.
INTEGER                      :: NumSteps                                 ! Number of time steps for the FFT.
INTEGER                      :: NumUSRf                                  ! Number of frequencies in the user-defined spectra
INTEGER                      :: NumUSRz                                  ! Number of heights defined in the user-defined profiles.
INTEGER                      :: RandSeed   (3)                           ! The array that holds the random seeds.
INTEGER,    ALLOCATABLE      :: RandSeedAry(:)                           ! The array that holds the random seeds.
INTEGER                      :: RandSeedTmp                              ! Holds the input random seed for the SNLWIND-3D RNG
INTEGER                      :: ScaleIEC                                 ! Flag to indicate if turbulence should be scaled to target value; 0 = NO scaling; 1 = scale based on hub; 2 = scale each point individually
INTEGER,    PARAMETER        :: UACT     = 14                            ! I/O unit for AeroDyn coherent turbulence
INTEGER,    PARAMETER        :: UACTTS   = 15                            ! I/O unit for coherent turbulence time step history file
INTEGER,    PARAMETER        :: UAFFW    = 9                             ! I/O unit for AeroDyn FF data (*.bts file).
INTEGER,    PARAMETER        :: UAHH     = 10                            ! I/O unit for AeroDyn HH data (*.hh  file).
INTEGER,    PARAMETER        :: UATWR    = 13                            ! I/O unit for AeroDyn tower data (*.twr file).
INTEGER,    PARAMETER        :: UBFFW    = 16                            ! I/O unit for BLADED FF data (*.wnd file).
INTEGER,    PARAMETER        :: UFFF     = 4                             ! I/O unit for formatted FF data.
INTEGER,    PARAMETER        :: UFTP     = 12                            ! I/O unit for formatted HH turbulence properties.
INTEGER,    PARAMETER        :: UGTP     = 11                            ! I/O unit for GenPro HH turbulence properties.
INTEGER,    PARAMETER        :: UI       = 1                             ! I/O unit for input file.
INTEGER,    PARAMETER        :: US       = 3                             ! I/O unit for summary file.
INTEGER,    PARAMETER        :: USpec    = 17                            ! I/O unit for user-defined spectra
INTEGER                      :: YLim                                     ! Number of horizontal positions in the grid
INTEGER                      :: ZLim                                     ! Number of vertical positions in the grid, plus extra hub point (if necessary), plus tower points

INTEGER,    PARAMETER        :: UC       = 22                            ! I/O unit for Coherence debugging file.
INTEGER,    PARAMETER        :: UD       = 20                            ! I/O unit for debugging data.
INTEGER,    PARAMETER        :: UP       = 21                            ! I/O unit for PSD debugging file.

LOGICAL                      :: Clockwise                                ! Flag to indicate clockwise rotation when looking downwind.
LOGICAL,    PARAMETER        :: COH_OUT  = .FALSE.                       ! This parameter has been added to replace the NON-STANDARD compiler directive previously used
LOGICAL,    PARAMETER        :: DEBUG_OUT= .FALSE.                       ! This parameter has been added to replace the NON-STANDARD compiler directive previously used
LOGICAL,    PARAMETER        :: PSD_OUT  = .FALSE. !                     ! This parameter has been added to replace the NON-STANDARD compiler directive previously used
LOGICAL,    PARAMETER        :: MVK      = .FALSE.                       ! This parameter has been added to replace the NON-STANDARD compiler directive previously used
LOGICAL                      :: ExtraHubPT                               ! Flag to indicate if the hub is on the regular grid or if an extra point must be added
LOGICAL                      :: ExtraTwrPT                               ! Flag to indicate if the tower is on the regular grid or if an extra point must be added

LOGICAL                      :: KHtest                                   ! Flag to indicate that turbulence should be extreme, to demonstrate effect of KH billows
LOGICAL                      :: NumTurbInp                               ! Flag to indicate if turbulence is user-specified (as opposed to IEC standard A, B, or C)
LOGICAL,    PARAMETER        :: PeriodicY = .FALSE. !.TRUE.
LOGICAL                      :: UVskip                                   ! Flag to determine if UV cross-feed term should be skipped or used
LOGICAL                      :: UWskip                                   ! Flag to determine if UW cross-feed term should be skipped or used
LOGICAL                      :: VWskip                                   ! Flag to determine if VW cross-feed term should be skipped or used
LOGICAL                      :: WrACT                                    ! Flag to output AeroDyn coherent turbulence
LOGICAL                      :: WrADFF                                   ! Flag to output AeroDyn FF data (binary).
LOGICAL                      :: WrADHH                                   ! Flag to output AeroDyn HH data (formatted).
LOGICAL                      :: WrADTWR                                  ! Flag to output AeroDyn tower data (binary).
LOGICAL                      :: WrBHHTP                                  ! Flag to output binary HH turbulence parameters.
LOGICAL                      :: WrBLFF                                   ! Flag to output BLADED FF data (binary)
LOGICAL                      :: WrFHHTP                                  ! Flag to output formatted HH turbulence parameters.
LOGICAL                      :: WrFmtFF                                  ! Flag to output formatted FF data (Traditional SNLWIND-3D format).

CHARACTER(200)               :: CTEventPath                              ! String used to store the name of the coherent event definition file
CHARACTER(200)               :: CTEventFile                              ! String used to store the name of the coherent event definition file
CHARACTER(  3)               :: CTExt                                    ! String used to determine the type of coherent structures ("dns" or "les")
CHARACTER(200)               :: DescStr                                  ! String used to describe the run (and the first line of the summary file)
CHARACTER(200)               :: FormStr                                  ! String used to store format specifiers.
CHARACTER(200)               :: FormStr1                                 ! String used to store format specifiers.
CHARACTER(200)               :: FormStr2                                 ! String used to store format specifiers.
CHARACTER( 23)               :: IECeditionStr (3) = &
                                (/'IEC 61400-1 Ed. 1: 1993', &
                                  'IEC 61400-1 Ed. 2: 1999', &
                                  'IEC 61400-1 Ed. 3: 2005'/)            ! The string description of the IEC 61400-1 standard being used
CHARACTER(  1)               :: IECTurbC                                 ! IEC turbulence characteristic.
CHARACTER(  1)               :: IECTurbE                                 ! IEC Extreme turbulence class.
CHARACTER( 35)               :: IEC_WindDesc                             ! The description of the IEC wind type
CHARACTER(200)               :: InFile = 'TurbSim.inp'                   ! Root name of the I/O files.
CHARACTER(  6)               :: RNG_type                                 ! Type of Random Number Generator to use
CHARACTER(197)               :: RootName                                 ! Root name of the I/O files.
CHARACTER( 50)               :: TMName                                   ! Turbulence model name.
CHARACTER(  6)               :: TurbModel                                ! Turbulence model.
CHARACTER(  3)               :: WindProfileType                          ! The wind profile type


TYPE                         :: Event                                    ! Coherent turbulent event to add to the background wind
   INTEGER                   :: EventNum                                 ! The event number (index into EventName() array)
   REAL(ReKi)                :: TStart                                   ! The time at which to add this event
   REAL(ReKi)                :: delt                                     ! The delta time before the event begins (for interpolation in AeroDyn)
   LOGICAL(1)                :: Connect2Prev = .FALSE.                   ! Whether this event is connected to the next, otherwise there is space between them
   TYPE(Event), POINTER      :: Next         => NULL()                   ! The next event to add
END TYPE

TYPE (Event), POINTER        :: PtrHead      => NULL()                   ! Pointer to the first event
TYPE (Event), POINTER        :: PtrTail      => NULL()                   ! Pointer to the last event

CONTAINS
!---------------------------------------------------------
SUBROUTINE get_coefs(JetHt,UH_coef,WD_coef)

      ! This subroutine just returns the coefficients that Neil calculated
      ! for getting the Chebyshev coefficients for jet wind profiles.
      
      ! The coefficients are
      ! Row 1 = Jet maximum wind speed coefficient
      ! Row 2 = Turbine layer Richardson number coefficient
      ! Row 3 = uStar over the rotor diameter coefficient
      ! Row 4 = constant coefficient
      ! Columns 1:11 = coefficients for 0-10th Chebyshev Basis Functions
      

   REAL(ReKi),INTENT(IN)     :: JetHt                                    ! The height of the jet
   REAL(ReKi),INTENT(OUT)    :: UH_coef(4,11)                            ! The coefficients for horizontal wind speed
   REAL(ReKi),INTENT(OUT)    :: WD_coef(4,11)                            ! The coefficients for (horizontal) wind direction

   INTEGER                   :: HtIndx

   HtIndx  = INT(JetHt - 50) / INT(20)
   HtIndx  = MIN( MAX( HtIndx, 1 ), 21 )

      ! The Horizontal Wind Speed coefficients
   SELECT CASE ( HtIndx )
      CASE ( 1 )     ! 70-90 m
         UH_coef(:, 1) = (/  0.856851,  7.51E-02,  1.39276,   0.894127 /)
         UH_coef(:, 2) = (/ -4.88E-02,  0.576344,  1.23582,   1.72687  /)
         UH_coef(:, 3) = (/  1.39E-02,  9.67E-02,  1.36737,  -0.723851 /)
         UH_coef(:, 4) = (/  0.100585,  0.234968, -1.06287,  -0.372353 /)
         UH_coef(:, 5) = (/ -7.69E-02, -0.154071, -0.301483,  0.150179 /)
         UH_coef(:, 6) = (/  8.53E-03,  0.104602, -0.382453,  0.520224 /)
         UH_coef(:, 7) = (/ -4.44E-03, -4.80E-02,  0.219135, -0.266775 /)
         UH_coef(:, 8) = (/  2.63E-02, -3.08E-02, -6.94E-02, -0.210521 /)
         UH_coef(:, 9) = (/ -2.01E-02, -5.61E-02,  0.220825,  0.179622 /)
         UH_coef(:,10) = (/  8.11E-03,  3.96E-02,  0.109793, -3.81E-02 /)
         UH_coef(:,11) = (/  4.99E-03,  5.00E-02, -0.124887, -0.11035  /)
      CASE ( 2 )     ! 90-110 m
         UH_coef(:, 1) = (/  0.741241, -0.122521,  0.875062,  1.43294  /)
         UH_coef(:, 2) = (/ -0.264131,  0.28827,   0.717571,  3.30541  /)
         UH_coef(:, 3) = (/ -5.92E-02,  3.86E-02,  1.09453,  -0.377399 /)
         UH_coef(:, 4) = (/  0.13792,   0.175628, -0.57163,  -0.539205 /)
         UH_coef(:, 5) = (/ -2.59E-02, -0.211126, -4.25E-02, -0.338308 /)
         UH_coef(:, 6) = (/ -1.02E-02,  0.153597, -0.197867,  0.570708 /)
         UH_coef(:, 7) = (/ -3.22E-02, -8.17E-02, -9.63E-02,  0.19095  /)
         UH_coef(:, 8) = (/  2.72E-02,  3.09E-02, -0.249399, -0.273684 /)
         UH_coef(:, 9) = (/ -1.60E-02,  8.88E-03,  0.132523,  9.58E-02 /)
         UH_coef(:,10) = (/ -5.29E-03,  2.98E-02,  0.205812,  9.27E-02 /)
         UH_coef(:,11) = (/  7.00E-03, -1.47E-02, -2.11E-02, -0.123083 /)
      CASE ( 3 )     ! 110-130 m
         UH_coef(:, 1) = (/  0.809492, -1.41752,  -0.817619,  1.64159  /)
         UH_coef(:, 2) = (/ -0.121866, -1.09012,  -2.60044,   3.63875  /)
         UH_coef(:, 3) = (/ -0.105142, -0.263657, -5.60E-02,  0.374811 /)
         UH_coef(:, 4) = (/  8.33E-02,  0.625103,  0.422112, -0.199598 /)
         UH_coef(:, 5) = (/ -1.69E-02, -7.09E-02,  1.76933,  -0.847721 /)
         UH_coef(:, 6) = (/  1.88E-02,  7.70E-02, -0.121062,  0.10533  /)
         UH_coef(:, 7) = (/ -3.15E-02,  2.50E-02, -7.39E-02,  0.299197 /)
         UH_coef(:, 8) = (/  3.48E-03,  4.25E-02, -6.52E-02, -4.29E-03 /)
         UH_coef(:, 9) = (/ -1.18E-02, -0.100754,  0.170602,  3.42E-02 /)
         UH_coef(:,10) = (/  2.09E-02,  3.36E-02, -0.104123, -8.49E-02 /)
         UH_coef(:,11) = (/ -2.91E-03, -3.52E-02, -0.258115,  4.81E-02 /)
      CASE ( 4 )     ! 130-150 m
         UH_coef(:, 1) = (/  0.694325, -0.463252,  2.11406,   1.28643  /)
         UH_coef(:, 2) = (/ -0.269118, -1.31381,   2.13374,   3.46187  /)
         UH_coef(:, 3) = (/ -8.40E-02, -5.97E-02,  2.09803,  -0.592335 /)
         UH_coef(:, 4) = (/  0.135657, -0.117732, -0.11134,  -0.28161  /)
         UH_coef(:, 5) = (/ -1.29E-02, -0.239685,  0.151264, -0.412806 /)
         UH_coef(:, 6) = (/  3.54E-02,  0.513824,  0.673662, -0.519536 /)
         UH_coef(:, 7) = (/ -1.55E-02,  7.49E-03,  0.393002,  2.07E-02 /)
         UH_coef(:, 8) = (/  2.37E-02,  0.225841,  3.84E-02, -0.202507 /)
         UH_coef(:, 9) = (/ -3.26E-02, -0.239615, -0.133893,  0.29135  /)
         UH_coef(:,10) = (/  1.52E-02,  7.15E-02,  0.25228,  -0.113016 /)
         UH_coef(:,11) = (/  7.19E-03,  9.79E-02,  0.252125, -0.173201 /)
      CASE ( 5 )     ! 150-170 m
         UH_coef(:, 1) = (/  0.909534,  0.581254, -2.90539,  -0.581377 /)
         UH_coef(:, 2) = (/  0.155834, -0.836954, -6.77075,   0.627044 /)
         UH_coef(:, 3) = (/ -8.99E-02, -5.28E-02, -2.0719,    2.44E-02 /)
         UH_coef(:, 4) = (/  7.01E-02, -0.152904, -0.348237,  0.460754 /)
         UH_coef(:, 5) = (/ -1.78E-02, -0.263166,  0.375798, -0.215738 /)
         UH_coef(:, 6) = (/  9.70E-03,  0.254932,  0.449286, -0.234    /)
         UH_coef(:, 7) = (/  7.46E-03, -0.304057, -0.122661, -7.14E-03 /)
         UH_coef(:, 8) = (/ -6.26E-03, -0.142341, -1.95E-02,  0.299841 /)
         UH_coef(:, 9) = (/ -2.59E-02,  0.174282,  0.193868, -5.81E-03 /)
         UH_coef(:,10) = (/  2.54E-03, -8.22E-02,  1.84E-02,  6.77E-02 /)
         UH_coef(:,11) = (/  5.77E-04, -5.43E-02, -7.69E-02,  2.96E-02 /)
      CASE ( 6 )     ! 170-190 m
         UH_coef(:, 1) = (/  0.885753, -1.15015,   0.155218, -0.707043 /)
         UH_coef(:, 2) = (/ -2.53E-02, -2.65126,   0.850151,  1.85279  /)
         UH_coef(:, 3) = (/ -7.23E-02, -0.399161,  0.142486, -0.917176 /)
         UH_coef(:, 4) = (/  3.78E-02,  0.178924,  0.227745,  0.528861 /)
         UH_coef(:, 5) = (/ -6.43E-03,  5.42E-02,  0.359052, -0.26111  /)
         UH_coef(:, 6) = (/  5.33E-02,  0.1546,   -0.335116, -0.602604 /)
         UH_coef(:, 7) = (/ -6.50E-03, -0.205907, -8.59E-02,  8.16E-02 /)
         UH_coef(:, 8) = (/  3.16E-02,  0.151199, -0.126411, -0.148609 /)
         UH_coef(:, 9) = (/ -3.95E-02,  0.127418,  0.158511,  0.20932  /)
         UH_coef(:,10) = (/ -2.53E-02, -5.32E-02,  0.36536,   0.214466 /)
         UH_coef(:,11) = (/  4.03E-03,  1.02E-02, -7.01E-03, -4.32E-02 /)
      CASE ( 7 )     ! 190-210 m
         UH_coef(:, 1) = (/  0.735269, -1.48574,   0.983734,  0.887351 /)
         UH_coef(:, 2) = (/  0.233065, -0.850536, -1.17754,  -0.880493 /)
         UH_coef(:, 3) = (/ -0.172346, -0.862128,  1.20075,   3.48E-02 /)
         UH_coef(:, 4) = (/  8.04E-02,  5.24E-02, -0.916548,  0.247144 /)
         UH_coef(:, 5) = (/  2.88E-02,  0.112064,  1.51E-04, -0.466186 /)
         UH_coef(:, 6) = (/ -2.75E-02, -9.01E-02, -0.321617,  0.379162 /)
         UH_coef(:, 7) = (/ -1.08E-02, -0.161368, -2.51E-04, -1.33E-02 /)
         UH_coef(:, 8) = (/  5.09E-02,  0.228507,  0.195942, -0.45807  /)
         UH_coef(:, 9) = (/ -1.98E-02, -7.23E-02,  6.66E-02,  0.133182 /)
         UH_coef(:,10) = (/ -5.57E-03, -5.31E-02,  2.44E-02,  5.60E-02 /)
         UH_coef(:,11) = (/  3.71E-03, -1.63E-02, -5.44E-02, -1.40E-02 /)
      CASE ( 8 )     ! 210-230 m
         UH_coef(:, 1) = (/  0.723721, -0.691359, -0.147971,  1.16041  /)
         UH_coef(:, 2) = (/  0.18799,   0.370199,  0.354538, -0.494962 /)
         UH_coef(:, 3) = (/ -0.204727, -0.166723,  0.682431,  0.367566 /)
         UH_coef(:, 4) = (/  1.40E-02,  0.334677,  0.169944,  0.494211 /)
         UH_coef(:, 5) = (/  3.84E-02,  0.258361,  0.389453, -0.625709 /)
         UH_coef(:, 6) = (/ -6.62E-03, -2.19E-02, -0.606278,  0.205521 /)
         UH_coef(:, 7) = (/ -2.54E-02, -0.17744,   7.49E-02,  7.61E-02 /)
         UH_coef(:, 8) = (/  5.03E-02,  7.97E-02, -9.98E-02, -0.312218 /)
         UH_coef(:, 9) = (/ -2.25E-02,  2.20E-02,  0.263227,  0.123311 /)
         UH_coef(:,10) = (/ -1.43E-02, -2.01E-02, -5.14E-02,  0.159391 /)
         UH_coef(:,11) = (/  2.64E-03,  3.46E-02, -0.12318,  -2.22E-02 /)
      CASE ( 9 )     ! 230-250 m
         UH_coef(:, 1) = (/  0.717665, -0.294178, -0.521541,  0.876418 /)
         UH_coef(:, 2) = (/  0.183182, -0.52658,  -1.34668,   0.414396 /)
         UH_coef(:, 3) = (/ -0.196162,  9.84E-02, -3.83E-02,  0.156018 /)
         UH_coef(:, 4) = (/  2.92E-02, -0.362193, -0.658593,  0.521854 /)
         UH_coef(:, 5) = (/  3.37E-02,  0.108203,  0.318667, -0.375309 /)
         UH_coef(:, 6) = (/ -8.24E-03,  0.128457, -0.149225,  0.1621   /)
         UH_coef(:, 7) = (/ -3.06E-02, -0.210106,  4.55E-02,  8.42E-02 /)
         UH_coef(:, 8) = (/  3.02E-02,  0.184626,  9.46E-02, -0.215191 /)
         UH_coef(:, 9) = (/  7.03E-03,  2.49E-02,  3.13E-02, -9.70E-02 /)
         UH_coef(:,10) = (/ -3.06E-03, -4.82E-02, -9.70E-02,  5.82E-02 /)
         UH_coef(:,11) = (/ -9.57E-03, -3.93E-02, -0.125623,  0.112639 /)
      CASE ( 10 )    ! 250-270 m
         UH_coef(:, 1) = (/  0.786229, -0.164848,  0.244948, -0.126263 /)
         UH_coef(:, 2) = (/  0.15218,  -0.153233, -0.558524,  0.84425  /)
         UH_coef(:, 3) = (/ -0.130716, -0.217411,  0.13439,  -0.536893 /)
         UH_coef(:, 4) = (/  1.70E-03,  5.49E-02,  0.551012,  0.335778 /)
         UH_coef(:, 5) = (/  2.47E-02,  2.82E-02,  0.290918, -0.223416 /)
         UH_coef(:, 6) = (/  1.48E-02,  5.94E-02, -0.277959,  3.91E-02 /)
         UH_coef(:, 7) = (/ -4.43E-02,  6.99E-03,  0.302386,  0.123719 /)
         UH_coef(:, 8) = (/  2.07E-02,  4.05E-02, -0.256155, -5.84E-02 /)
         UH_coef(:, 9) = (/  4.51E-03, -4.37E-02, -0.111911, -9.20E-03 /)
         UH_coef(:,10) = (/  4.05E-03, -6.90E-03,  0.14697,  -7.03E-02 /)
         UH_coef(:,11) = (/ -6.68E-03,  1.53E-02, -2.55E-02,  4.97E-02 /)
      CASE ( 11 )    ! 270-290 m
         UH_coef(:, 1) = (/  0.715734, -0.772062, -0.556396,  1.02929  /)
         UH_coef(:, 2) = (/  0.322509, -0.465616, -0.671711, -1.2413   /)
         UH_coef(:, 3) = (/ -0.166728, -0.281268,  0.924893, -0.282907 /)
         UH_coef(:, 4) = (/  1.27E-02, -0.342767, -1.10823,   0.516431 /)
         UH_coef(:, 5) = (/  3.80E-02,  5.35E-03,  0.833719, -0.510102 /)
         UH_coef(:, 6) = (/  1.97E-02,  0.279705, -0.179026, -4.36E-02 /)
         UH_coef(:, 7) = (/ -4.74E-02, -0.227673,  9.00E-02,  0.341958 /)
         UH_coef(:, 8) = (/  8.99E-03, -1.92E-02, -0.433969,  5.90E-02 /)
         UH_coef(:, 9) = (/  4.34E-03,  8.12E-02,  0.25764,  -0.148492 /)
         UH_coef(:,10) = (/  1.03E-02,  3.24E-02,  0.141971, -0.105207 /)
         UH_coef(:,11) = (/ -4.84E-03, -1.99E-02,  7.33E-02,  2.84E-02 /)
      CASE ( 12 )    ! 290-310 m
         UH_coef(:, 1) = (/  0.723348, -0.289581, -1.10618,   0.970713 /)
         UH_coef(:, 2) = (/  0.283383,  1.12986,  -0.152861, -0.653269 /)
         UH_coef(:, 3) = (/ -0.16513,   0.295047,  0.245326, -7.06E-02 /)
         UH_coef(:, 4) = (/  8.55E-03,  9.38E-02, -0.826824,  0.283436 /)
         UH_coef(:, 5) = (/  3.45E-02,  0.364581,  0.566317, -0.521081 /)
         UH_coef(:, 6) = (/  2.83E-02,  0.107252, -0.124867, -4.80E-02 /)
         UH_coef(:, 7) = (/ -3.57E-02, -0.230151, -6.88E-02,  0.231208 /)
         UH_coef(:, 8) = (/  5.62E-04,  1.40E-02, -0.334942,  0.121313 /)
         UH_coef(:, 9) = (/ -6.35E-03, -6.19E-02,  0.139396,  2.77E-02 /)
         UH_coef(:,10) = (/  1.14E-02, -2.67E-02,  0.24201,  -0.127337 /)
         UH_coef(:,11) = (/  1.71E-04, -6.37E-04,  4.39E-02, -5.61E-03 /)
      CASE ( 13 )    ! 310-330 m
         UH_coef(:, 1) = (/  0.736987, -0.103727,  9.95E-02,  0.343208 /)
         UH_coef(:, 2) = (/  0.28285,   0.370583,  1.17749,  -0.490259 /)
         UH_coef(:, 3) = (/ -0.130451, -0.557928, -0.272771, -0.230816 /)
         UH_coef(:, 4) = (/ -1.83E-02,  1.00E-01, -0.367321,  0.486971 /)
         UH_coef(:, 5) = (/  2.66E-02, -0.149206,  0.365342, -0.318809 /)
         UH_coef(:, 6) = (/  4.16E-02,  3.60E-02, -0.801161,  6.00E-06 /)
         UH_coef(:, 7) = (/ -2.36E-02,  1.96E-04,  0.340449,  2.72E-02 /)
         UH_coef(:, 8) = (/  1.30E-03,  0.214384,  0.125371, -8.47E-02 /)
         UH_coef(:, 9) = (/ -1.23E-02,  4.75E-02,  0.182118,  1.78E-02 /)
         UH_coef(:,10) = (/  4.63E-03, -0.1309,   -0.130584,  2.35E-02 /)
         UH_coef(:,11) = (/  9.03E-04, -6.18E-02, -7.85E-03,  1.17E-02 /)
      CASE ( 14 )    ! 330-350 m
         UH_coef(:, 1) = (/  0.706488, -1.21766,   1.08617,   0.674247 /)
         UH_coef(:, 2) = (/  0.341777,  2.27476,   3.81434,  -2.32363  /)
         UH_coef(:, 3) = (/ -0.112822,  7.53E-02,  0.221349, -0.700428 /)
         UH_coef(:, 4) = (/ -1.99E-02, -1.95E-02,  0.947788,  4.68E-02 /)
         UH_coef(:, 5) = (/  3.08E-02,  0.334947,  0.10847,  -0.534662 /)
         UH_coef(:, 6) = (/  5.21E-02,  0.349056, -1.14517,  -0.147474 /)
         UH_coef(:, 7) = (/ -1.67E-02, -0.143994, -0.409398,  0.228081 /)
         UH_coef(:, 8) = (/ -1.75E-03, -0.115198,  3.23E-03,  0.100094 /)
         UH_coef(:, 9) = (/ -2.30E-02, -5.63E-02,  0.168561,  0.159537 /)
         UH_coef(:,10) = (/ -6.41E-03, -8.48E-02,  0.135087,  8.81E-02 /)
         UH_coef(:,11) = (/  1.13E-03,  2.07E-02,  9.18E-02, -3.77E-02 /)
      CASE ( 15 )    ! 350-370 m
         UH_coef(:, 1) = (/  0.721629, -0.941544,  0.923908,  0.543678 /)
         UH_coef(:, 2) = (/  0.346956, -0.281582, -2.32358,  -0.244435 /)
         UH_coef(:, 3) = (/ -0.109484,  0.275053,  0.86928,  -0.771081 /)
         UH_coef(:, 4) = (/ -3.96E-02, -0.790621, -8.84E-02,  0.723378 /)
         UH_coef(:, 5) = (/  1.59E-02, -0.394222, -0.479505, -8.67E-02 /)
         UH_coef(:, 6) = (/  2.68E-02,  0.466895,  0.522378, -0.263669 /)
         UH_coef(:, 7) = (/ -9.57E-03, -8.52E-02,  1.11E-02,  3.20E-02 /)
         UH_coef(:, 8) = (/  3.46E-04, -5.34E-02,  0.15998,   0.108225 /)
         UH_coef(:, 9) = (/ -1.10E-02, -0.116864, -6.06E-02,  6.09E-02 /)
         UH_coef(:,10) = (/ -2.93E-03,  2.72E-02,  5.08E-02,  7.50E-03 /)
         UH_coef(:,11) = (/ -2.04E-03, -2.07E-02, -3.07E-02,  3.58E-02 /)
      CASE ( 16 )    ! 370-390 m
         UH_coef(:, 1) = (/  0.732127, -2.66819,  -7.94E-02,  0.676096 /)
         UH_coef(:, 2) = (/  0.285167,  3.89442,  -0.917426,  0.104248 /)
         UH_coef(:, 3) = (/ -8.38E-02,  0.235268, -2.19E-03, -0.914663 /)
         UH_coef(:, 4) = (/ -3.98E-02, -0.858603, -0.538194,  0.843739 /)
         UH_coef(:, 5) = (/ -1.64E-02,  0.287007, -5.39E-02,  0.108834 /)
         UH_coef(:, 6) = (/  3.31E-02,  0.218726,  0.175636, -0.329844 /)
         UH_coef(:, 7) = (/  3.10E-05, -6.89E-02,  3.76E-02, -4.73E-02 /)
         UH_coef(:, 8) = (/  1.06E-02, -5.03E-02,  1.99E-02,  3.74E-02 /)
         UH_coef(:, 9) = (/ -1.05E-02,  9.92E-02,  0.11293,   2.26E-02 /)
         UH_coef(:,10) = (/ -2.99E-03, -0.106831,  0.122628,  1.83E-02 /)
         UH_coef(:,11) = (/ -7.32E-03,  3.52E-02, -3.36E-02,  8.59E-02 /)
      CASE ( 17 )    ! 390-410 m
         UH_coef(:, 1) = (/  0.707698,  0.119876,  0.427545,  0.2468   /)
         UH_coef(:, 2) = (/  0.307273,  0.428003, -3.09224,   1.01117  /)
         UH_coef(:, 3) = (/ -7.33E-02,  0.51572,  -0.229086, -0.792402 /)
         UH_coef(:, 4) = (/ -4.73E-02,  8.49E-02, -0.52415,   0.571084 /)
         UH_coef(:, 5) = (/ -2.83E-02,  0.165455, -0.691726,  0.349932 /)
         UH_coef(:, 6) = (/  2.17E-02,  0.258434,  0.170597, -0.236707 /)
         UH_coef(:, 7) = (/ -4.59E-03, -0.130722,  0.182955, -3.40E-02 /)
         UH_coef(:, 8) = (/  1.82E-02,  9.79E-02,  0.189511, -0.158597 /)
         UH_coef(:, 9) = (/ -7.84E-04, -2.50E-02,  0.137171, -5.77E-02 /)
         UH_coef(:,10) = (/ -2.91E-03, -4.84E-02,  0.168698,  8.22E-03 /)
         UH_coef(:,11) = (/ -4.67E-03,  1.75E-03,  1.80E-02,  4.41E-02 /)
      CASE ( 18 )    ! 410-430 m
         UH_coef(:, 1) = (/  0.688761, -0.7286,   -1.55711,   1.27145  /)
         UH_coef(:, 2) = (/  0.300421,  0.633115,  0.881706, -8.38E-03 /)
         UH_coef(:, 3) = (/ -6.81E-02,  0.210301,  0.610772, -0.714435 /)
         UH_coef(:, 4) = (/ -5.93E-02, -0.373997, -0.593894,  1.01556  /)
         UH_coef(:, 5) = (/ -4.26E-02, -2.45E-02, -0.400705,  0.399717 /)
         UH_coef(:, 6) = (/  1.39E-02,  6.09E-02, -0.161239, -3.06E-02 /)
         UH_coef(:, 7) = (/ -4.41E-03, -1.98E-02,  0.293288, -0.110401 /)
         UH_coef(:, 8) = (/  1.42E-02,  8.22E-02, -1.50E-02, -1.54E-02 /)
         UH_coef(:, 9) = (/  6.30E-03, -1.50E-02, -7.57E-02, -7.10E-02 /)
         UH_coef(:,10) = (/  2.19E-03, -2.59E-02,  8.53E-02, -2.29E-02 /)
         UH_coef(:,11) = (/ -2.76E-03,  1.68E-02, -8.77E-02,  3.27E-02 /)
      CASE ( 19 )    ! 430-450 m
         UH_coef(:, 1) = (/  0.659495, -0.22327,  -1.75403,   1.65777  /)
         UH_coef(:, 2) = (/  0.384097,  1.06351,   2.53779,  -1.63428  /)
         UH_coef(:, 3) = (/ -2.42E-02,  0.113735, -1.42805,  -0.690773 /)
         UH_coef(:, 4) = (/ -3.30E-02,  8.60E-02, -1.00836,   0.764307 /)
         UH_coef(:, 5) = (/ -2.76E-02,  0.297567,  0.697445, -0.187071 /)
         UH_coef(:, 6) = (/  1.21E-02,  0.212621, -0.570822,  1.23E-02 /)
         UH_coef(:, 7) = (/ -2.22E-02,  0.166286,  0.50751,   1.87E-02 /)
         UH_coef(:, 8) = (/  1.52E-02,  5.81E-02, -0.256912, -5.10E-02 /)
         UH_coef(:, 9) = (/  2.11E-03, -1.45E-02, -8.94E-02, -2.00E-02 /)
         UH_coef(:,10) = (/  3.06E-03,  1.60E-02,  7.45E-02, -3.77E-02 /)
         UH_coef(:,11) = (/ -1.84E-04, -1.56E-02, -6.25E-02,  1.57E-02 /)
      CASE ( 20 )    ! 450-470 m
         UH_coef(:, 1) = (/  0.64099,  -2.02496,   0.427597,  1.52166  /)
         UH_coef(:, 2) = (/  0.391609,  2.03441,  -0.122486, -1.03579  /)
         UH_coef(:, 3) = (/  8.28E-03,  0.5942,   -0.42469,  -1.35655  /)
         UH_coef(:, 4) = (/ -2.54E-02, -0.826812, -0.812187,  0.911776 /)
         UH_coef(:, 5) = (/ -2.77E-02, -9.73E-03,  0.315974,  2.34E-02 /)
         UH_coef(:, 6) = (/  1.37E-02,  0.365984,  0.141952, -0.299349 /)
         UH_coef(:, 7) = (/ -1.95E-02, -0.406182,  2.32E-02,  0.184752 /)
         UH_coef(:, 8) = (/  7.34E-03,  8.54E-02, -0.255458,  7.08E-02 /)
         UH_coef(:, 9) = (/  1.54E-03,  5.82E-02, -5.72E-02, -6.37E-02 /)
         UH_coef(:,10) = (/  5.11E-03, -6.11E-02, -7.04E-03, -3.64E-02 /)
         UH_coef(:,11) = (/  1.97E-03, -1.09E-02, -8.18E-02, -6.03E-03 /)
      CASE ( 21 )    ! 470-490 m
         UH_coef(:, 1) = (/  0.547127, -0.327778,  2.00666,   2.67869  /)
         UH_coef(:, 2) = (/  0.427112,  8.56E-02, -1.61197,  -1.17989  /)
         UH_coef(:, 3) = (/  6.23E-02,  0.760714, -0.659927, -2.30882  /)
         UH_coef(:, 4) = (/ -4.04E-02, -0.873328, -0.118326,  1.19626  /)
         UH_coef(:, 5) = (/ -4.85E-03,  0.130813, -0.169613, -0.181674 /)
         UH_coef(:, 6) = (/  4.82E-03,  0.289038,  7.34E-02,  6.45E-03 /)
         UH_coef(:, 7) = (/ -2.49E-02, -0.375342,  0.15139,   0.208253 /)
         UH_coef(:, 8) = (/  9.48E-04,  5.23E-02, -0.213227,  0.137941 /)
         UH_coef(:, 9) = (/ -9.18E-03,  3.91E-02,  7.26E-02,  4.73E-02 /)
         UH_coef(:,10) = (/ -6.00E-05,  1.03E-02,  7.46E-03,  1.86E-02 /)
         UH_coef(:,11) = (/ -2.21E-03, -9.70E-05, -7.13E-02,  4.29E-02 /)
      CASE DEFAULT
         CALL ProgAbort ('Error getting UH coefficients' )
   END SELECT

   SELECT CASE ( HtIndx )
      CASE ( 1 )     ! 70-90 m
         WD_coef(:, 1) = (/  5.07735,   96.4785,   18.8465,   110.986    /)
         WD_coef(:, 2) = (/  0.75209,  -16.5103,  -25.9592,     9.05636  /)
         WD_coef(:, 3) = (/ -1.50806,    1.69319,  -7.7859,    13.3041   /)
         WD_coef(:, 4) = (/  1.11287,    3.711,    13.1084,   -11.9491   /)
         WD_coef(:, 5) = (/ -0.987363,  -2.93059,  -4.75454,    9.04282  /)
         WD_coef(:, 6) = (/  0.65727,    0.560223, -0.541911,  -5.33397  /)
         WD_coef(:, 7) = (/ -0.493572,  -0.455574,  2.03972,    3.53745  /)
         WD_coef(:, 8) = (/  0.244207,   0.390402,  1.5338,    -1.9793   /)
         WD_coef(:, 9) = (/ -1.26E-02,   0.19732,  -2.70454,    0.179412 /)
         WD_coef(:,10) = (/  9.13E-04,   9.65E-02,  0.304467,   4.79E-02 /)
         WD_coef(:,11) = (/ -7.71E-02,  -0.11096,   0.51028,    0.585717 /)
      CASE ( 2 )     ! 90-110 m
         WD_coef(:, 1) = (/  2.98622,   87.1045,   41.7453,   124.301    /)
         WD_coef(:, 2) = (/  0.241282, -10.9238,  -31.5696,    11.0764   /)
         WD_coef(:, 3) = (/ -0.380786,  -1.71395,  -8.35561,    3.68007  /)
         WD_coef(:, 4) = (/  0.287014,   6.76407,  17.1736,    -7.4345   /)
         WD_coef(:, 5) = (/ -0.682991,  -5.48805, -12.7947,    10.9313   /)
         WD_coef(:, 6) = (/  0.415999,   2.36938,   4.47285,   -5.47595  /)
         WD_coef(:, 7) = (/ -0.184533,  -7.04E-02,  0.81309,    1.06891  /)
         WD_coef(:, 8) = (/  0.152381,  -0.344921,  3.40496,   -1.81465  /)
         WD_coef(:, 9) = (/ -0.113556,  -1.02575,  -5.54619,    2.51668  /)
         WD_coef(:,10) = (/  3.87E-02,   1.0794,    0.98668,   -0.942351 /)
         WD_coef(:,11) = (/  7.37E-02,  -0.284347,  1.12315,   -1.04163  /)
      CASE ( 3 )     ! 110-130 m
         WD_coef(:, 1) = (/ -10.8064,   63.1523,   18.7751,   255.252    /)
         WD_coef(:, 2) = (/  1.89875,  -15.7662,  -27.2545,    -5.90699  /)
         WD_coef(:, 3) = (/ -1.81141,   -7.58E-03,  4.49E-02,  19.4007   /)
         WD_coef(:, 4) = (/ -0.420216,   4.54261,  16.6642,    -1.5632   /)
         WD_coef(:, 5) = (/  3.09E-02,   0.162346, -5.68196,    1.70168  /)
         WD_coef(:, 6) = (/  0.372585,  -0.888944, -0.400871,  -3.98736  /)
         WD_coef(:, 7) = (/  0.137532,  -1.86E-02, -1.97659,   -1.07897  /)
         WD_coef(:, 8) = (/  7.11E-02,   0.275322,  2.06716,   -0.99703  /)
         WD_coef(:, 9) = (/ -0.142081,   0.690143,  1.74256,    0.963168 /)
         WD_coef(:,10) = (/ -0.225792,  -0.215169,  0.660299,   1.89319  /)
         WD_coef(:,11) = (/  1.91E-02,  -0.23,     -1.69222,    0.190668 /)
      CASE ( 4 )     ! 130-150 m
         WD_coef(:, 1) = (/  0.270461, 107.786,   140.705,    143.549    /)
         WD_coef(:, 2) = (/  2.46519,   25.9261,   54.6629,   -43.2182   /)
         WD_coef(:, 3) = (/ -1.11746,   -4.09287,  -5.71316,   16.4144   /)
         WD_coef(:, 4) = (/ -0.104557,   2.88836,  14.657,     -5.58632  /)
         WD_coef(:, 5) = (/  1.4104,    -0.862421,  1.88282,  -13.3856   /)
         WD_coef(:, 6) = (/ -0.994103,   6.07897,   6.16378,    6.53327  /)
         WD_coef(:, 7) = (/  0.440338,  -7.14173, -12.2957,     0.653282 /)
         WD_coef(:, 8) = (/ -0.705677,   2.13336,   2.39331,    5.62277  /)
         WD_coef(:, 9) = (/  0.398742,  -3.5049,   -3.97854,   -1.68531  /)
         WD_coef(:,10) = (/ -7.72E-02,   2.14124,   3.42657,   -0.982025 /)
         WD_coef(:,11) = (/  0.120525,  -1.80518,  -3.44124,    0.391772 /)
      CASE ( 5 )     ! 150-170 m
         WD_coef(:, 1) = (/  10.3894,  203.711,    87.9736,     0.818669 /)
         WD_coef(:, 2) = (/  4.15105,   37.734,    56.1061,   -72.0928   /)
         WD_coef(:, 3) = (/ -1.60031,   -6.42686,   2.99983,   21.7355   /)
         WD_coef(:, 4) = (/  0.162421, -22.7335,    4.23498,    0.433394 /)
         WD_coef(:, 5) = (/ -1.00817,   -1.82237, -17.2291,    18.8346   /)
         WD_coef(:, 6) = (/  0.591051,   5.30019,  22.1782,   -15.2786   /)
         WD_coef(:, 7) = (/ -0.350898,  -1.35238, -14.9057,     9.09022  /)
         WD_coef(:, 8) = (/  0.512704,   5.33682,  12.0501,   -11.3284   /)
         WD_coef(:, 9) = (/ -0.294613,  -6.61282, -13.756,      9.48747  /)
         WD_coef(:,10) = (/  0.180824,   6.67558,   8.1748,    -6.39538  /)
         WD_coef(:,11) = (/ -0.168678,  -3.5973,   -2.92266,    3.62255  /)
      CASE ( 6 )     ! 170-190 m
         WD_coef(:, 1) = (/ -3.05838,   92.242,    -6.17694,  218.678    /)
         WD_coef(:, 2) = (/ -1.19176,   10.9436,    5.33317,   23.6574   /)
         WD_coef(:, 3) = (/  0.396791,   5.36609,  14.86,     -12.1807   /)
         WD_coef(:, 4) = (/ -0.260044,  -3.3155,   -1.83325,    3.07872  /)
         WD_coef(:, 5) = (/  0.147588,   3.54423,   2.61624,   -2.87076  /)
         WD_coef(:, 6) = (/ -3.09E-02,  -0.298005, -3.99378,    2.512    /)
         WD_coef(:, 7) = (/  3.52E-02,   0.476622,  0.917889,  -1.19482  /)
         WD_coef(:, 8) = (/ -0.10397,   -3.13393,  -1.34654,    2.38467  /)
         WD_coef(:, 9) = (/  0.111959,   0.768005,  1.09164,   -1.84864  /)
         WD_coef(:,10) = (/ -5.32E-02,  -0.753046,  0.517477,   0.77376  /)
         WD_coef(:,11) = (/  2.36E-02,  -0.255733, -0.765475,  -0.183366 /)
      CASE ( 7 )     ! 190-210 m
         WD_coef(:, 1) = (/  2.63747,   48.8574, -148.839,    198.635    /)
         WD_coef(:, 2) = (/  0.276349,   8.15568,  11.5466,     4.89475  /)
         WD_coef(:, 3) = (/ -0.161153,  -3.92434,  15.2465,    -2.75263  /)
         WD_coef(:, 4) = (/ -0.215546,  -6.05707,  -0.221136,   2.96778  /)
         WD_coef(:, 5) = (/ -0.174687,   0.722833,  2.58751,    1.43519  /)
         WD_coef(:, 6) = (/ -3.24E-03,   0.841219,  2.36677,   -0.541046 /)
         WD_coef(:, 7) = (/ -0.14379,   -0.422125,  6.03272,   -3.55E-02 /)
         WD_coef(:, 8) = (/  4.94E-02,  -0.165447, -1.64947,   -0.118004 /)
         WD_coef(:, 9) = (/  6.88E-03,   0.618011,  0.600728,  -0.312735 /)
         WD_coef(:,10) = (/ -2.96E-02,  -0.102388, -0.423526,   0.526055 /)
         WD_coef(:,11) = (/  3.77E-03,  -0.79762,  -1.48591,    0.487559 /)
      CASE ( 8 )     ! 210-230 m
         WD_coef(:, 1) = (/  1.25931,   81.7121,  -72.2497,   192.288    /)
         WD_coef(:, 2) = (/ -0.421425,   0.812039, 26.4136,    12.7087   /)
         WD_coef(:, 3) = (/ -0.477334,  -0.804493, 10.2938,     2.63738  /)
         WD_coef(:, 4) = (/  0.27025,   -1.48414,   6.44E-02,  -3.62925  /)
         WD_coef(:, 5) = (/ -0.206555,   2.60212,   4.78E-03,   1.41829  /)
         WD_coef(:, 6) = (/  0.199714,  -0.145286, -1.43609,   -1.0421   /)
         WD_coef(:, 7) = (/ -8.81E-02,  -1.11826,   0.562309,   0.568182 /)
         WD_coef(:, 8) = (/  4.38E-02,  -0.94946,  -1.20199,    0.184361 /)
         WD_coef(:, 9) = (/ -5.13E-02,  -0.157795, -0.596316,   0.747777 /)
         WD_coef(:,10) = (/  5.03E-02,   6.23E-02, -0.821348,  -0.411198 /)
         WD_coef(:,11) = (/ -2.45E-02,   3.66E-03,  0.61934,    0.147334 /)
      CASE ( 9 )     ! 230-250 m
         WD_coef(:, 1) = (/  4.99773,   45.439,   -22.9981,   142.166    /)
         WD_coef(:, 2) = (/  1.34923,   -0.690733,  1.11037,   -7.00256  /)
         WD_coef(:, 3) = (/ -4.58E-02,  -1.48399,   3.15438,   -1.20619  /)
         WD_coef(:, 4) = (/ -5.86E-02,  -0.324401, -0.520264,   0.827308 /)
         WD_coef(:, 5) = (/  6.67E-02,   1.95293,  -1.46579,   -1.66186  /)
         WD_coef(:, 6) = (/  2.23E-02,   1.10257,   1.61038,   -0.14154  /)
         WD_coef(:, 7) = (/  4.83E-02,  -0.46633,   0.318096,  -1.22718  /)
         WD_coef(:, 8) = (/ -3.56E-02,  -0.905797, -0.659337,   1.10221  /)
         WD_coef(:, 9) = (/ -6.54E-04,   0.514329,  0.38488,   -0.221416 /)
         WD_coef(:,10) = (/  2.40E-03,  -0.307029, -0.455799,   0.167602 /)
         WD_coef(:,11) = (/  5.79E-03,  -0.3575,   -6.82E-02,  -1.79E-02 /)
      CASE ( 10 )    ! 250-270 m
         WD_coef(:, 1) = (/  2.87491,   81.7603,  -14.221,    143.973    /)
         WD_coef(:, 2) = (/  0.176626,   0.711168, 14.3778,     3.41781  /)
         WD_coef(:, 3) = (/ -0.112353,  -4.44334,   5.01439,   -0.539061 /)
         WD_coef(:, 4) = (/  0.135496,   0.868787, -2.54952,   -1.4882   /)
         WD_coef(:, 5) = (/ -5.87E-02,   7.34E-02,  0.618705,   0.341871 /)
         WD_coef(:, 6) = (/  4.36E-02,   1.16076,  -2.2411,     0.371484 /)
         WD_coef(:, 7) = (/ -4.21E-03,  -0.219162,  3.07613,   -1.48294  /)
         WD_coef(:, 8) = (/  2.91E-02,  -7.90E-02, -2.06058,    0.637811 /)
         WD_coef(:, 9) = (/  6.84E-04,   0.398542, -0.227958,  -0.195655 /)
         WD_coef(:,10) = (/ -1.33E-02,  -0.148014,  0.112677,   0.28039  /)
         WD_coef(:,11) = (/  4.56E-02,  -0.4372,   -1.05259,   -0.39506  /)
      CASE ( 11 )    ! 270-290 m
         WD_coef(:, 1) = (/ -3.74E-02,   5.72313, -25.8459,   204.708    /)
         WD_coef(:, 2) = (/  0.387587,   5.70337,  37.0722,    -5.10619  /)
         WD_coef(:, 3) = (/  0.130067,   8.86213,   7.6219,    -6.77984  /)
         WD_coef(:, 4) = (/ -1.83E-02,  -4.80402,   1.26728,    1.1988   /)
         WD_coef(:, 5) = (/ -0.125984,   5.69111,  -2.4798,     0.370193 /)
         WD_coef(:, 6) = (/  7.02E-02,  -4.02809,   0.545202,   0.396538 /)
         WD_coef(:, 7) = (/ -4.89E-02,   1.99119,  -7.47E-02,  -0.617665 /)
         WD_coef(:, 8) = (/  7.28E-02,  -1.94844,  -0.9012,     0.174322 /)
         WD_coef(:, 9) = (/ -2.75E-02,   0.875895,  8.29E-02,   1.47E-02 /)
         WD_coef(:,10) = (/ -4.90E-03,  -0.26505,   0.684299,  -0.101304 /)
         WD_coef(:,11) = (/ -2.46E-03,  -9.03E-02, -0.25124,    0.130552 /)
      CASE ( 12 )    ! 290-310 m
         WD_coef(:, 1) = (/  4.48806,  101.681,   -24.2152,   108.849    /)
         WD_coef(:, 2) = (/  1.12228,  -11.8153,   -5.83094,   -3.59506  /)
         WD_coef(:, 3) = (/  0.152934,   0.610899, 10.1148,    -6.59595  /)
         WD_coef(:, 4) = (/  6.76E-02,   1.44362,  -8.36227,    1.70741  /)
         WD_coef(:, 5) = (/ -8.86E-02,   1.22016,   4.89384,   -1.422    /)
         WD_coef(:, 6) = (/  1.14E-02,  -0.801065, -4.6529,     2.29577  /)
         WD_coef(:, 7) = (/ -5.68E-03,  -0.156515,  3.48364,   -1.85745  /)
         WD_coef(:, 8) = (/  3.21E-02,   0.643855, -1.80571,    0.499593 /)
         WD_coef(:, 9) = (/ -5.96E-03,  -0.645,     1.0105,    -0.256849 /)
         WD_coef(:,10) = (/ -1.79E-02,   0.137457, -7.45E-03,   0.232805 /)
         WD_coef(:,11) = (/ -5.07E-04,  -1.20E-03, -0.280138,   9.13E-02 /)
      CASE ( 13 )    ! 310-330 m
         WD_coef(:, 1) = (/  0.253568,  43.3822,   42.3741,   166.917    /)
         WD_coef(:, 2) = (/ -0.210713,  14.3161,   12.187,      9.66539  /)
         WD_coef(:, 3) = (/  0.176871,  -3.28688,  -2.78059,   -1.64384  /)
         WD_coef(:, 4) = (/  0.30952,    2.34743,  -5.8261,    -3.72051  /)
         WD_coef(:, 5) = (/ -0.211586,  -1.38792,  -0.891686,   3.26282  /)
         WD_coef(:, 6) = (/  0.114874,  -1.0177,   -2.95833,   -0.285227 /)
         WD_coef(:, 7) = (/ -0.168163,   1.33608,   5.32715,    0.270668 /)
         WD_coef(:, 8) = (/  0.106821,   0.746965, -1.28128,   -1.11127  /)
         WD_coef(:, 9) = (/ -2.17E-02,   0.198171,  0.911532,   2.31E-02 /)
         WD_coef(:,10) = (/ -5.64E-03,   0.278658,  0.250055,  -9.16E-02 /)
         WD_coef(:,11) = (/  7.21E-03,   2.24E-02,  6.76E-02,  -0.1011   /)
      CASE ( 14 )    ! 330-350 m
         WD_coef(:, 1) = (/  1.4365,   104.113,    86.7884,   138.082    /)
         WD_coef(:, 2) = (/  1.01951,  -22.4231,    8.14651,   -3.0374   /)
         WD_coef(:, 3) = (/ -0.14238,    5.5217,   -8.37098,    1.9052   /)
         WD_coef(:, 4) = (/ -8.04E-02,   2.56411,   8.01756,    0.450076 /)
         WD_coef(:, 5) = (/  7.34E-03,  -3.31792,  -10.0037,    1.66433  /)
         WD_coef(:, 6) = (/ -3.82E-02,   3.00083,   6.14358,   -0.656165 /)
         WD_coef(:, 7) = (/  0.113861,  -4.41267,  -2.98194,   -1.24882  /)
         WD_coef(:, 8) = (/ -0.154066,   4.29174,   3.74587,    1.4816   /)
         WD_coef(:, 9) = (/  0.127996,  -2.88696,  -2.49795,   -1.24336  /)
         WD_coef(:,10) = (/ -6.71E-02,   1.70388,   0.935254,   0.748082 /)
         WD_coef(:,11) = (/  8.19E-03,  -4.50E-02, -0.263839,  -5.18E-02 /)
      CASE ( 15 )    ! 350-370 m
         WD_coef(:, 1) = (/ -0.675054, 121.016,     0.173435, 199.751    /)
         WD_coef(:, 2) = (/ -0.52795,   26.7663,   36.6465,     8.14164  /)
         WD_coef(:, 3) = (/  0.686068,  -2.58652,   1.37125,  -12.8021   /)
         WD_coef(:, 4) = (/ -0.115391,  -0.715049,  0.225913,   2.68255  /)
         WD_coef(:, 5) = (/  0.127924,   1.18619,  -3.81934,   -2.40047  /)
         WD_coef(:, 6) = (/ -0.201212,  -1.51136,   4.51548,    3.23679  /)
         WD_coef(:, 7) = (/  0.175571,  -0.664591, -5.74074,   -2.24143  /)
         WD_coef(:, 8) = (/ -0.107098,   0.889236,  3.25149,    1.18349  /)
         WD_coef(:, 9) = (/  3.15E-02,  -6.48E-02, -0.882842,  -0.404645 /)
         WD_coef(:,10) = (/ -9.69E-03,  -0.486174, -0.284323,   0.336898 /)
         WD_coef(:,11) = (/  1.04E-03,  -0.144399, -6.10E-02,   6.62E-02 /)
      CASE ( 16 )    ! 370-390 m
         WD_coef(:, 1) = (/  0.610558, -90.3161,  -86.1311,   221.346    /)
         WD_coef(:, 2) = (/ -0.878196,   0.234356, -1.96802,   30.3835   /)
         WD_coef(:, 3) = (/  0.536954,   2.31986,   0.611791, -11.624    /)
         WD_coef(:, 4) = (/ -0.203843,  -2.10521,  -1.77538,    5.20693  /)
         WD_coef(:, 5) = (/ -6.04E-02,  -1.53784,   0.391834,   1.09004  /)
         WD_coef(:, 6) = (/ -3.32E-02,   1.08307,   0.756223,   0.579045 /)
         WD_coef(:, 7) = (/  2.20E-03,   1.00851,   0.872176,  -1.24302  /)
         WD_coef(:, 8) = (/ -4.70E-02,   0.313443, -5.20E-02,   1.24129  /)
         WD_coef(:, 9) = (/  0.105906,   2.60251,  -0.805126,  -2.35033  /)
         WD_coef(:,10) = (/ -3.95E-02,  -0.866726,  0.244709,   0.996069 /)
         WD_coef(:,11) = (/  5.34E-02,   0.423689, -0.910358,  -0.888237 /)
      CASE ( 17 )    ! 390-410 m
         WD_coef(:, 1) = (/ -0.256694, -53.0924,  -28.899,    212.286    /)
         WD_coef(:, 2) = (/  0.368178,   0.200188,-15.1321,     9.40209  /)
         WD_coef(:, 3) = (/ -0.102825,  -4.83546,   9.24228,   -0.64019  /)
         WD_coef(:, 4) = (/  0.191961,   2.99238,  -4.8869,    -2.80575  /)
         WD_coef(:, 5) = (/ -9.33E-02,   0.237869,  3.72573,   -8.03E-02 /)
         WD_coef(:, 6) = (/  1.70E-02,   2.22246,  -0.874,      0.324301 /)
         WD_coef(:, 7) = (/ -4.39E-02,  -1.22545,   1.03253,   -7.41E-02 /)
         WD_coef(:, 8) = (/  9.07E-03,  -0.438369, -1.85468,    0.746178 /)
         WD_coef(:, 9) = (/ -2.97E-02,  -0.626331,  1.32958,    0.161941 /)
         WD_coef(:,10) = (/ -4.73E-03,  -0.639604, -0.50062,    0.398523 /)
         WD_coef(:,11) = (/  7.78E-04,   0.203885,  0.111938,  -9.66E-02 /)
      CASE ( 18 )    ! 410-430 m
         WD_coef(:, 1) = (/ -1.05454,   19.3432,   14.3866,   209.914    /)
         WD_coef(:, 2) = (/ -5.37E-02,  -6.69143,  -5.48868,   13.8188   /)
         WD_coef(:, 3) = (/  0.130461,   1.84379,  10.2975,    -6.85151  /)
         WD_coef(:, 4) = (/  0.120135,   3.25255,  -4.64527,   -0.957415 /)
         WD_coef(:, 5) = (/ -0.157071,  -1.87681,   4.37492,    1.52585  /)
         WD_coef(:, 6) = (/  0.220174,   1.14707,  -5.27774,   -2.10403  /)
         WD_coef(:, 7) = (/ -0.185849,  -8.73E-02,  4.5702,     1.45097  /)
         WD_coef(:, 8) = (/  5.77E-02,  -0.265271, -2.17262,    1.19E-02 /)
         WD_coef(:, 9) = (/ -3.19E-02,   0.159054,  1.11463,    9.91E-02 /)
         WD_coef(:,10) = (/ -9.31E-03,  -0.514427, -0.486658,   0.472324 /)
         WD_coef(:,11) = (/  5.84E-03,  -6.98E-02, -6.53E-02,  -7.68E-02 /)
      CASE ( 19 )    ! 430-450 m
         WD_coef(:, 1) = (/  0.624689,  63.9533, -115.139,    203.718    /)
         WD_coef(:, 2) = (/ -0.249911,   8.56489,  12.0426,    11.2274   /)
         WD_coef(:, 3) = (/  0.208499,  -2.38494,   8.76157,   -7.17681  /)
         WD_coef(:, 4) = (/ -0.205812,   3.60713,   5.60652,    2.51439  /)
         WD_coef(:, 5) = (/  0.320606,  -7.16713, -10.6408,    -3.32927  /)
         WD_coef(:, 6) = (/ -0.178674,   5.15743,   3.70481,    2.92097  /)
         WD_coef(:, 7) = (/  0.101549,  -5.22916,  -1.89887,   -1.64557  /)
         WD_coef(:, 8) = (/ -9.30E-02,   2.8729,    1.14221,    1.4604   /)
         WD_coef(:, 9) = (/  1.45E-02,  -1.29998,  -0.491218,  -6.91E-02 /)
         WD_coef(:,10) = (/ -6.95E-04,   0.830442,  1.25591,   -0.451134 /)
         WD_coef(:,11) = (/ -6.90E-04,   1.30E-02, -0.16423,    7.65E-02 /)
      CASE ( 20 )    ! 450-470 m
         WD_coef(:, 1) = (/  4.30205,   83.823,   -77.8869,   120.115    /)
         WD_coef(:, 2) = (/  0.11147,   -2.13123, -13.0305,    11.4506   /)
         WD_coef(:, 3) = (/  5.36E-02,  -9.82942,   3.21203,   -2.14437  /)
         WD_coef(:, 4) = (/  3.12E-02,  -0.694,    -2.56494,    0.846492 /)
         WD_coef(:, 5) = (/ -3.97E-02,   0.628515,  0.898384,  -0.403596 /)
         WD_coef(:, 6) = (/  0.187725,  -1.32489,  -3.10108,   -1.64756  /)
         WD_coef(:, 7) = (/ -8.75E-02,  -0.750003,  1.2358,     0.95118  /)
         WD_coef(:, 8) = (/  4.29E-02,   0.206995, -0.591777,  -0.495133 /)
         WD_coef(:, 9) = (/ -3.25E-02,   0.187007,  0.351131,   0.374602 /)
         WD_coef(:,10) = (/ -1.79E-02,  -0.651232, -0.437205,   0.653204 /)
         WD_coef(:,11) = (/  5.74E-03,   0.210108, -0.185616,  -8.91E-02 /)
      CASE ( 21 )    ! 470-490 m
         WD_coef(:, 1) = (/  0.685959,  76.5757,  -26.8137,   187.31     /)
         WD_coef(:, 2) = (/ -0.229648,   3.36903, -12.3466,    19.5787   /)
         WD_coef(:, 3) = (/  5.56E-02,  -6.33886,   2.64958,   -2.35925  /)
         WD_coef(:, 4) = (/ -3.42E-02,  -1.78314,   1.51304,    0.43034  /)
         WD_coef(:, 5) = (/  5.81E-02,   4.2818,   -1.08668,   -2.13185  /)
         WD_coef(:, 6) = (/ -1.94E-02,  -2.76039,  -0.573698,   1.97694  /)
         WD_coef(:, 7) = (/  1.26E-02,   0.932315,  0.974862,  -1.5273   /)
         WD_coef(:, 8) = (/  1.04E-02,  -0.143063, -0.728002,   0.464589 /)
         WD_coef(:, 9) = (/  1.21E-03,   0.262702, -0.133363,  -0.236706 /)
         WD_coef(:,10) = (/ -2.29E-04,  -0.162697, -0.138587,   0.17236  /)
         WD_coef(:,11) = (/  6.61E-03,  -5.47E-02, -0.104054,  -9.64E-02 /)
      CASE DEFAULT
         CALL ProgAbort ('Error getting WD coefficients' )
   END SELECT

   RETURN
END SUBROUTINE get_coefs

END MODULE TSMods
!=======================================================================
