MODULE CommonData
   !! Module for common parameters and variables.
   USE types
IMPLICIT NONE

    !! file for input and output
    character(200)          :: outfile, molfile, specref, moldat
    character(*), PARAMETER :: radat   = '../data/'
    character(*), PARAMETER :: version = 'version 2024'
    character(*), PARAMETER :: logfile = './radex.log'

    !Escape probability method (uncomment your choice)
    INTEGER  :: method=1 !1=uniform sphere,2=LVG, 3=slab

    !! how the initial value of the matrix is defined either ussing Tbg or Tkin or (combined Tkin and Tbg)
    INTEGER  :: imethod  !! imethod = 1-classic(tbg); imethod=2-modern(tkin), imethod=3-mixed
    REAL(dp) :: tol      !! Relative tolerance
    REAL(dp) :: alpha    !! under-relaxation factor to stabilize the convergence

    !!---------------------------------------------------------
    !!Physical and astronomical constants (CODATA 2002)
    REAL(dp), PARAMETER :: clight  = 2.99792458d10  !! speed of light     (cm/s)
    REAL(dp), PARAMETER :: hplanck = 6.6260963d-27  !! Planck constant    (erg/Hz)
    REAL(dp), PARAMETER :: kboltz  = 1.3806505d-16  !! Boltzmann constant (erg/K)
    REAL(dp), PARAMETER :: pi      = 3.14159265d0   !! pi
    REAL(dp), PARAMETER :: amu     = 1.67262171d-24 !! atomic mass unit   (g)
    REAL(dp), PARAMETER :: tcmb    = 2.725          !! CMB background temperature (K)
    !---------------------------------------------------------
      
    !Array sizes
    INTEGER, PARAMETER :: maxpart = 9     !! maximum no. of collision partners (seven defined)
    INTEGER, PARAMETER :: maxtemp = 99    !! maximum no. of collision temperatures
    INTEGER, PARAMETER :: maxlev  = 2999  !! maximum no. of energy levels
    INTEGER, PARAMETER :: maxline = 99999 !! maximum no. of radiative transitions
    INTEGER, PARAMETER :: maxcoll = 99999 !! maximum no. of collisional transitions
    INTEGER  :: nthick                !! counts optically thick lines
    !---------------------------------------------------------
    !!   Molecular data
    !---------------------------------------------------------
    INTEGER :: nlev             !! nlev : actual number of levels
    INTEGER :: nline            !! nline: actual number of lines
    INTEGER :: ncoll            !! ncoll: actual number of transitions
    INTEGER :: npart            !! npart: actual number of partners
    INTEGER :: ntemp            !! actual number of collision temperatures
    INTEGER :: iupp(maxline)    !! upper level of line i
    INTEGER :: ilow(maxline)    !! lower level of line i

    REAL(dp) :: amass           !! molecular mass              (amu)
    REAL(dp) :: eterm(maxlev)   !! energy levels               (1/cm)
    REAL(dp) :: gstat(maxlev)   !! statistical weights
    REAL(dp) :: aeinst(maxline) !! Einstein A coefficients     (1/s)
    REAL(dp) :: eup(maxline)    !! line upper level energy     (K)

    ! colld :  downward rate coefficients  (cm^3 /s)
    ! xpop  :  level populations

    !--------------------------------------------------------- 
    !! Physical conditions
    !---------------------------------------------------------
    REAL(dp) :: density(maxpart)  !! density:  number densities of collision partners  (cm^-3)
    REAL(dp) :: tkin              !! tkin   :  kinetic temperature                     (K)
    REAL(dp) :: tbg               !! tbg    :  temperature of background radiation     (K)
    REAL(dp) :: cdmol             !! cdmol  :  molecular column density                (cm^-2)
    REAL(dp) :: deltav            !! deltav :  FWHM line width                         (cm/s)
    REAL(dp) :: totdens           !! totdens:  total number density of all partners    (cm^-3)

    !---------------------------------------------------------      
    !!  Numerical parameters
    !---------------------------------------------------------
    INTEGER, PARAMETER  :: miniter = 10     !! minimum number of iterations
    INTEGER, PARAMETER  :: maxiter = 9999   !! maximum number of iterations
    REAL(dp)            :: fmin, fmax       !! minimum/maximum output frequency
    REAL(dp)            :: ccrit            !! relative tolerance on solution
    REAL(dp), PARAMETER :: eps    = 1.0d-30 !! round-off error
    REAL(dp), PARAMETER :: minpop = 1.0d-20 !! minimum level population

    !---------------------------------------------------------
    !! Radiative quantities
    !---------------------------------------------------------
    REAL(dp) :: taul(maxline)               !! line optical depth
    REAL(dp) :: backi(maxline)              !! background intensity [erg s-1 cm-2 Hz-1 sr-1]
    REAL(dp) :: xnu(maxline)                !! line frequency (cm^-1)
    REAL(dp) :: trj(maxline)                !! background brightness (RJ)
    REAL(dp) :: totalb(maxline)             !! background temperature (BB)
    REAL(dp) :: spfreq(maxline)             !! spectroscopic line frequency (GHz), not used in
                                            !! calculation but only to print output
    REAL(dp) :: antennaTemp(maxline)        !! line antenna temperature
    REAL(dp) :: upperPops(maxline)          !! upper level populations of line i
    REAL(dp) :: lowerPops(maxline)          !! lower level populations of line i
    REAL(dp) :: wavelength(maxline)         !! line wavelength (micron)
    REAL(dp) :: intensityKkms(maxline)      !! line integrated intensity (K km/s)
    REAL(dp) :: intensityErgs(maxline)      !! line flux (erg / s / cm^2)

    CHARACTER(6) :: qnum(maxlev)                          !! quantum numbers of levels
    CHARACTER(6) :: lowQNum(maxlev)                       !! quantum numbers of the lower level of the line
    CHARACTER(6) :: upperQNum(maxlev)                     !! quantum numbers of the upper level of the line
    REAL(dp), PARAMETER :: fk    = hplanck*clight/kboltz  !! help to calculate intensities
    REAL(dp), PARAMETER :: thc   = 2.d0*hplanck*clight    !! help to calculate intensities
    REAL(dp), PARAMETER :: fgaus = 1.0645*8.0*pi          !! accounts for Gaussian line shape

    ! tex  :  line excitation temperature (solution)
    ! btex :  line excitation temperature using the classic method
    ! itex :  line excitation temperature using the modern method
    ! lext :  line excitation temperature

    !---------------------------------------------------------
    !!   Collisional quantities
    !---------------------------------------------------------
    REAL(dp) :: ctot(maxlev)                !! ctot  : total collision rate 
    REAL(dp) :: crate(maxlev,maxlev)        !! crate : collision rate matrix (density * rate coefficient)
    REAL(dp) :: xpop(maxlev)                !! xpop  : level populations

    !For development / maintenance purposes:
    LOGICAL, PARAMETER :: debug =.false.

END MODULE CommonData
