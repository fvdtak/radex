MODULE CommonParams
   !! Module for common parameters and variables.
   USE types
IMPLICIT NONE

    !! file for input and output
    character(200)          :: outfile, molfile, specref
    character(*), PARAMETER :: radat   = 'data/'
    character(*), PARAMETER :: version = 'version 2024'
    character(*), PARAMETER :: logfile = './nradex.log'

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
!    REAL(dp) :: xpop(maxlev)                !! xpop  : level populations

    !For development / maintenance purposes:
    LOGICAL, PARAMETER :: debug =.false.

Contains
SUBROUTINE CalcOutputArrays(tex,  xpop, nlines)
    REAL(dp) :: tex(maxline)         !! Exicitation temperature
    INTEGER, INTENT(INOUT) :: nlines !! number of output lines
    REAL(dp), intent(out) :: xpop(maxlev)                !! xpop  : level populations

    !! local variables
    INTEGER :: iline    ! to loop over lines
    INTEGER :: m,n      ! upper & lower level of the line

    REAL(dp) ::  xt        ! frequency cubed
    REAL(dp) ::  hnu       ! photon energy
    REAL(dp) ::  bnutex    ! line source function
    REAL(dp) ::  ftau      ! exp(-tau)
    REAL(dp) ::  toti      ! background intensity
    REAL(dp) ::  tbl       ! black body temperature
    REAL(dp) ::  wh        ! Planck correction
    REAL(dp) ::  tback     ! background temperature
    REAL(dp) ::  ta        ! line antenna temperature
    REAL(dp) ::  tr        ! line radiation temperature
    REAL(dp) ::  beta ! escape probability
    REAL(dp) ::  bnu       ! Planck function
    REAL(dp) ::  kkms      ! line integrated intensity (K km/s)
    REAL(dp) ::  ergs      ! line flux (erg / s / cm^2)
    REAL(dp) ::  wavel     ! line wavelength (micron)
    nlines=0
    DO iline=1,nline
      m  = iupp(iline)
      n  = ilow(iline)
      xt = xnu(iline)**3.
      !Calculate source function
      hnu = fk*xnu(iline)/tex(iline)
      IF (hnu.ge.160.0d0) THEN
        bnutex = 0.0d0
      else
        bnutex = thc*xt/(dexp(fk*xnu(iline)/tex(iline))-1.d0)
      END IF
      !Calculate line brightness in excess of background
      ftau = 0.0d0
      IF (abs(taul(iline)).le.3.d2) ftau = dexp(-taul(iline))
      toti = backi(iline)*ftau+bnutex*(1.d0-ftau)
      IF (toti .eq. 0.0d0) THEN
        tbl = 0.0d0
      else
        wh = thc*xt/toti+1.d0
        if (wh .le. 0.d0) then
          tbl = toti/(thc*xnu(iline)*xnu(iline)/fk)
        else
          tbl = fk*xnu(iline)/dlog(wh)
        end if
      END IF
      if (backi(iline) .eq. 0.0d0) then
        tback = 0.0d0
      else
        tback = fk*xnu(iline)/dlog(thc*xt/backi(iline)+1.d0)
      end if
      !Calculate antenna temperature
      tbl = tbl-tback
      hnu = fk*xnu(iline)
      if (abs(tback/hnu) .le. 0.02) then
        ta = toti
      else
        ta = toti-backi(iline)
      endif
      ta = ta/(thc*xnu(iline)*xnu(iline)/fk)
      !Calculate radiation temperature
      beta = escprob(taul(iline))
      bnu  = totalb(iline)*beta+(1.d0-beta)*bnutex
      IF (bnu .eq. 0.0d0) THEN
        tr = totalb(iline)
      else
        wh = thc*xt/bnu + 1.0
        IF (wh .le. 0.0) THEN
          tr = bnu/(thc*xnu(iline)*xnu(iline)/fk)
        else
          tr = fk*xnu(iline)/dlog(wh)
        END IF
      END IF

      !! unit =  micron
      wavel = clight / spfreq(iline) / 1.0d5
      kkms  = 1.0645*deltav*ta
      ergs  = fgaus*kboltz*deltav*ta*(xnu(iline)**3.)
      antennaTemp(iline)   = ta
      wavelength(iline)    = wavel
      lowerPops(iline)     = (xpop(n))
      upperPops(iline)     = (xpop(m))
      intensityKkms(iline) = (kkms/1.0d5)
      intensityErgs(iline) = ergs
      lowQNum(iline)       = qnum(n)
      upperQNum(iline)     = qnum(m)
                nlines     = nlines+1
    END DO

  END SUBROUTINE CalcOutputArrays

   !! Calculate the escape probability
  FUNCTION EscProb(tau)

    REAL(dp) :: tau
    REAL(dp) :: EscProb, beta
    REAL(dp) :: taur  !optical radius

    !!new tau
    taur = tau/2.0

    SELECT CASE (method)
      CASE(1)
        !Uniform sphere formula from Osterbrock (Astrophysics of
        !Gaseous Nebulae and Active Galactic Nuclei) Appendix 2
        !with power law approximations for large and small tau
        IF (abs(taur) .lt. 0.1) THEN
          beta = 1.d0-0.75d0*taur+(taur**2.)/2.5d0&
            &-(taur**3.)/6.d0+(taur**4.)/17.5d0
        ELSE IF(abs(taur) .gt. 5.d1) THEN
            beta = 0.75d0/taur
        ELSE
          beta = 0.75d0/taur*(1.d0-1.d0/(2.d0*(taur**2.))+&
              &(1.d0/taur+1.d0/(2.d0*(taur**2.)))*dexp(-2.*taur))
        END IF


      CASE(2)
        !Expanding sphere = Large Velocity Gradient (LVG) or Sobolev case.
        !Formula from De Jong, Boland and Dalgarno (1980, A&A 91, 68)
        !corrected by factor 2 in order to match EscProb(TAU=0)=1
        IF (abs(taur) .lt. 0.01) THEN
          beta = 1.0
        else if(abs(taur) .lt. 7.0) THEN
          beta = 2.0*(1.0 - dexp(-2.34*taur))/(4.68*taur)
        else
          beta = 2.0/(taur*4.0*(sqrt(log(taur/sqrt(pi)))))
        END IF

      CASE (3)
        !Slab geometry (e.g., shocks): de Jong, Dalgarno & Chu 1975, 
        !ApJ 199, 69 (again with power law approximations)
        IF (abs(3.0*tau) .lt. 0.1) THEN
          beta = 1.0 - 1.5*(tau + tau**2.)
        ELSE IF (abs(3.0*tau) .gt. 50.0) THEN
          beta = 1.0d0/(3.0*tau)
        ELSE
          beta = (1.0d0 - dexp(-3.0*tau))/(3.0*tau)
        END IF
      CASE DEFAULT 
         WRITE(*,*) 'Error: Escape probability method undefined'
         STOP
    END SELECT
    EscProb = beta
    RETURN

  END FUNCTION EscProb

END MODULE CommonParams
