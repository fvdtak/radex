MODULE IO
   !! Module for input and output.
   !! {!src/io.md!}
   use CommonData
   use Solver
   IMPLICIT NONE

CONTAINS
  SUBROUTINE ReadData(status)
    INTEGER, INTENT(INOUT) :: status

    ! Reads molecular data files (2003 format)
    INTEGER  :: ilev,jlev                     !! to loop over energy levels
    INTEGER  :: iline                         !! to loop over lines
    INTEGER  :: ipart,jpart                   !! to loop over collision partners
    INTEGER  :: itemp                         !! to loop over collision temperatures
    INTEGER  :: icoll                         !! to loop over collisional transitions
    INTEGER  :: dummy                         !! to skip part of the file
    INTEGER  :: id(maxpart)                   !! to identify collision partners
    INTEGER  :: lcu(maxcoll),lcl(maxcoll)     !! upper/lower levels of collisional transition
    REAL(dp) :: coll(maxpart,maxcoll,maxtemp) !! collisional 3D arrays
    REAL(dp) :: colld(maxpart,maxlev,maxlev)  !! collisional 3D arrays after interpolation
    REAL(dp) :: ediff                         !! different between the energy levels
    REAL(dp) :: temp(maxtemp)                 !! collision temperatures
    character*120 :: collref                  !! text about source of collisional data

    !! to interpolate rate coeffs
    INTEGER  :: iup,ilo,nint
    REAL(dp) :: tupp,tlow,fint

    !! to verify matching of densities and rates
    logical found

    !! to calculate thermal o/p ratio for H2
    REAL(dp) :: opr

    !! Executable part begins here.
    WRITE(*,*) "Data file: " // molfile
    OPEN(unit=11,file=molfile,status='old',err=99)
    !! in the header, every second line is a comment
    101 format(a)
    READ(11,*) 
    READ(11,101) specref
    READ(11,*) 
    READ(11,*) amass
    READ(11,*) 
    READ(11,*) nlev

    IF (nlev.lt.1) THEN
      WRITE(*,*) 'error: too few energy levels defined'
      status=0
      RETURN
    END IF
    IF (nlev.gt.maxlev) THEN
      WRITE(*,*) 'error: too many energy levels defined' 
      status=0
      RETURN
    END IF
    IF (debug) write(*,*) 'readdata: basics'

    !! Term energies and statistical weights
    READ(11,*)
    DO ilev=1,nlev
     READ(11,*) dummy,eterm(ilev),gstat(ilev),qnum(ilev)
     IF ((dummy.lt.1).or.(dummy.gt.nlev)) stop 'error:illegal level number'
    END DO

    IF (debug) write(*,*) 'readdata: levels'!,(eterm(ilev),ilev=1,nlev)

    !! Radiative upper & lower levels and Einstein coefficients
    READ(11,*) 
    READ(11,*) nline
    READ(11,*)

    IF (nline .lt. 1) stop 'error: too few spectral lines defined' 
    IF (nline .gt. maxline) stop 'error: too many spectral lines defined'

    DO iline = 1, nline
     READ(11,*) dummy,iupp(iline),ilow(iline),aeinst(iline)&
              &,spfreq(iline),eup(iline)
     IF ((dummy .lt. 1) .or. (dummy .gt. nline)) STOP 'error:illegal line number'
     xnu(iline)=(eterm(iupp(iline))-eterm(ilow(iline)))
     IF ((xnu(iline) .lt. eps)) stop 'error:illegal line frequency'
    END DO

    IF (debug) write(*,*) 'readdata: lines',(xnu(iline),iline=1,nline)

    ! Number of collision partners
    READ(11,*)
    READ(11,*) npart
    IF (npart .lt. 1) stop 'error: too few collision partners defined'
    IF (npart .gt. maxpart) stop 'error: too many collision partners'

    102  format(i1,a)
    DO ipart=1,npart
     READ(11,*)
     READ(11,102) id(ipart),collref 
     READ(11,*)
     READ(11,*) ncoll
    IF (ncoll.lt.1) stop 'error: too few collision rates defined'
    IF (ncoll.gt.maxcoll) stop 'error: too many collision rates'
     READ(11,*)
     READ(11,*) ntemp
     IF (ntemp.lt.0) stop 'error: too few collision temperatures defined'
     IF (ntemp.gt.maxtemp) stop 'error: too many collision temperatures'
     READ(11,*)
     READ(11,*) (temp(itemp),itemp=1,ntemp)
     READ(11,*)

     IF (debug) write(*,*) 'ready to read ',ncoll,' rates for '&
            &,ntemp,' temperatures for partner ',ipart

     DO icoll=1,ncoll
        READ(11,*) dummy,lcu(icoll),lcl(icoll),&
               &(coll(ipart,icoll,itemp),itemp=1,ntemp)
        IF ((dummy .lt. 1) .or. (dummy .gt. ncoll)) stop 'error:illegal collision number'
     END DO

     ! interpolate array coll(ncol,ntemp) to desired temperature
     ! Must DO this now because generally, rates with different partners
     ! are calculated for a different set of temperatures
     IF (ntemp.le.1) THEN
        DO icoll=1,ncoll
           iup=lcu(icoll)
           ilo=lcl(icoll)
           colld(ipart,iup,ilo)=coll(ipart,icoll,1)
        END DO
     else
        IF (tkin.gt.temp(1)) THEN
           IF (tkin.lt.temp(ntemp)) THEN
              !===  interpolation :
              DO itemp=1,(ntemp-1)
                 IF (tkin.gt.temp(itemp).and.tkin.le.temp(itemp+1)) nint=itemp
              END DO
              tupp=temp(nint+1)
              tlow=temp(nint)
              fint=(tkin-tlow)/(tupp-tlow)

              DO icoll=1,ncoll
                 iup=lcu(icoll)
                 ilo=lcl(icoll)
                 colld(ipart,iup,ilo)=coll(ipart,icoll,nint)&
                        &+fint*(coll(ipart,icoll,nint+1)&
                        &-coll(ipart,icoll,nint))
                 IF (colld(ipart,iup,ilo).lt.0.0)&
                       & colld(ipart,iup,ilo)=coll(ipart,icoll,nint)
              END DO
           else
              !===  Tkin too high :
              IF (tkin.ne.temp(ntemp)) THEN
                  WRITE(*,*)' Warning : Tkin higher than temperatures '/&
                       &/'for which rates are present.'
              END IF
              DO icoll=1,ncoll
                 iup=lcu(icoll)
                 ilo=lcl(icoll)
                 colld(ipart,iup,ilo)=coll(ipart,icoll,ntemp)
              END DO
           END IF
        else
    ! Tkin too low :
           IF (tkin.ne.temp(1)) THEN
              WRITE(*,*)' Warning : Tkin lower than temperatures '//&
                     &'for which rates are present.'
           END IF
           DO icoll=1,ncoll
              iup=lcu(icoll)
              ilo=lcl(icoll)
              colld(ipart,iup,ilo)=coll(ipart,icoll,1)
           END DO
        END IF
     END IF

    ! Finished reading rate coefficients

    END DO

    IF (debug) write(*,*) 'readdata: rate coeffs'

    close(11)

    ! Combine rate coeffs of several partners, multiplying by partner density.
    DO iup=1,nlev
      DO ilo=1,nlev
        crate(iup,ilo)=0.0d0
      END DO
    END DO

    totdens = 0.0d0
    found   = .false.

    ! Special case (CO, atoms): user gives total H2 but data file has o/p-H2.
    ! Quite a big IF:
    IF ((npart .gt. 1) .and.&
        &(density(1).gt.eps).and.(density(2).lt.eps)&
         &.and.(density(3).lt.eps)) THEN
       opr        = min(3.d0,9.0*dexp(-170.6/tkin))
       density(2) = density(1)/(opr+1.d0)
       density(3) = density(1)/(1.d0+1.d0/opr)
       !WRITE(*,*)'*** Warning: Assuming thermal o/p ratio for H2 of ',opr
    END IF
    ! Note that for files like CN which have H2 and e- rates, the
    ! warning is given without reason. May fix this later if needed.

    DO ipart=1,maxpart
      totdens = totdens + density(ipart)
      DO jpart=1,maxpart
        IF ((id(jpart).eq.ipart).and.(density(ipart).gt.0.d0)) THEN
           found = .true.
           DO iup=1,nlev
              DO ilo=1,nlev
                 crate(iup,ilo) = crate(iup,ilo) +&
                        &density(ipart)*colld(jpart,iup,ilo)
              END DO
           END DO
        END IF
      END DO
    END DO

    IF (.not.found) THEN
     WRITE(*,*)'*** Warning: No rates found for any collision partner'
     WRITE(*,*)'It is likely your parameter dictionary does not have a non-zero density for any partner'
     WRITE(*,*)'in the collisional data file'
     status=0
     RETURN
    END IF

    ! Calculate upward rates from detailed balance
    DO iup = 1,nlev
      DO ilo = 1,nlev
        ediff = eterm(iup)-eterm(ilo)
        IF (ediff.gt.0.0d0) THEN
          IF ((fk*ediff/tkin).ge.160.0d0) THEN
              crate(ilo,iup) = 0.0d0
            else
              crate(ilo,iup) = gstat(iup)/gstat(ilo)&
            &              *dexp(-fk*ediff/tkin)*crate(iup,ilo)
          END IF
        END IF
      END DO
      ! initialize ctot array
      ctot(iup) = 0.0d0
    END DO

    ! Calculate total collision rates (inverse collisional lifetime)
    DO ilev=1,nlev
      DO jlev=1,nlev
        ctot(ilev)=crate(ilev,jlev)+ctot(ilev)
      END DO
   END DO

    RETURN
    !  IF (debug) WRITE(*,*)'ctot=',(ctot(ilev),ilev=1,nlev)
    99 write(*,*) 'error opening data file ',molfile
    STOP
  END SUBROUTINE ReadData

  SUBROUTINE GetInputs(status)
    INTEGER :: status

    !local variables
    INTEGER :: ipart        ! loop over collision partners
    CHARACTER(10) :: partner ! name of collision partner
    INTEGER :: id           ! ID code of collision partner
    CHARACTER(len=255) :: datdir

    status = 0
    !Set input parameters to default values
    CALL defaults

    !Write log file for later reference
    OPEN(unit=13,file=logfile,status='unknown',err=97)

    20   format(A)
    !must read file names in format(A): in free format, slashes are treated as separators
    21   format(A,$)
    22   format(A,I2,A,$)
    23   format(1pe10.3)
    24   format(1pe10.3,2x,1pe10.3)
    25   format(i2)
    
    CALL get_environment_variable("PWD", datdir)
    
    write(*,21) 'Molecular data file ? '
    READ(*,20) molfile
    moldat = trim(molfile)

    IF ((molfile(1:1).ne.'/').and.(molfile(1:1).ne.'.')) &
    &     molfile = datdir(1:length(datdir))//'/data/'//molfile(1:length(molfile))
    write(13,20) molfile(1:length(molfile))

    write(*,21) 'Name of output file ? '
    READ(*,20) outfile
    write(13,20) outfile(1:length(outfile))

    write(*,21) 'Minimum and maximum output frequency [GHz] ? '
    READ(*,*) fmin,fmax
    !Default values: DC / Lyman limit
    IF (fmin.eq.fmax) THEN
       fmin = 0.d0
       fmax = 3.d7    
    END IF
    IF (fmin.gt.fmax) THEN
       fmin = fmin + fmax
       fmax = fmin - fmax
       fmin = fmin - fmax
    END IF
    fmin = dmax1(0.d0,fmin)
    fmax = dmin1(3.d7,fmax)
    write(13,24) fmin,fmax

    !Read kinetic temp and loop until acceptable value entered
    WRITE(*,21) 'Kinetic temperature [K] ?  '
    READ(*,*) tkin
    DO WHILE ((tkin.lt.0.1).or.(tkin.gt.1.e4))
      WRITE(*,*) 'Please enter a value between 0.1 and 1e4'
      WRITE(*,21) 'Kinetic temperature [K] ?  '
      READ(*,*) tkin
    END DO
    write(13,23) tkin

    !Read N collisional partners and loop until acceptable value entered
    WRITE(*,21) 'Number of collision partners ?  '
    READ(*,*) npart
    DO WHILE ((npart.lt.1).or.(npart.gt.7))
      WRITE(*,*)'Please enter a value between 1 and 7'
      WRITE(*,21) 'Number of collision partners ?  '
      READ(*,*) npart
    END DO
    WRITE(13,25) npart

    !Loop over the collisional partners and read types
    DO ipart=1,npart
      WRITE(*,22) 'Type of partner',ipart,' ? '
      id = 0
      DO WHILE (id .eq. 0)
        READ(*,*) partner

        IF ((partner.eq.'h2').or.(partner.eq.'H2')) id=1
        IF ((partner(1:1).eq.'p').or.(partner(1:1).eq.'P')) id=2
        IF ((partner(1:1).eq.'o').or.(partner(1:1).eq.'O')) id=3
        IF ((partner(1:1).eq.'e').or.(partner(1:1).eq.'E')) id=4
        IF ((partner.eq.'h').or.(partner.eq.'H')) id=5
        IF ((partner.eq.'he').or.(partner.eq.'He')) id=6
        IF ((partner.eq.'h+').or.(partner.eq.'H+')) id=7
        IF (id.eq.0) WRITE(*,*) 'Unknown species. Choose from: H2 p-H2 o-H2 e H He H+'
      END DO
      
      WRITE(13,*) partner
      WRITE(*,22) 'Density of collision partner ',ipart,' [cm^-3] ? '
      READ(*,*) density(id)
      DO WHILE ((density(id).lt.1.0e-3).or.(density(id).gt.1.0e13))
        WRITE(*,*)'Please enter a value between 1e-3 and 1e13'
        WRITE(*,22) 'Density of collision partner ',ipart,' [cm^-3] ? '
        READ(*,*) density(id)
      END DO
      WRITE(13,23) density(id)
    END DO

    !Add ortho and para H2 densities if applicable
    IF ((density(2).gt.0.0).or.(density(3).gt.0.0)) &
      &density(1)=density(2)+density(3)

    WRITE(*,21) 'Background temperature [K] ?  '
    READ(*,*) tbg
    !Tbg > 0 means single blackbody such as CMB
    !Tbg = 0 means average ISRF
    !Tbg < 0 means use user-supplied routine
    DO WHILE ((tbg.lt.-1.e4).or.(tbg.gt.1.e4))
      WRITE(*,*)'Please enter a value between -1e4 and 1e4'
      WRITE(*,21) 'Background temperature [K] ?  '
      READ(*,*) tbg
    END DO
    WRITE(13,23) tbg

    WRITE(*,21) 'Molecular column density [cm^-2] ?  '
    READ(*,*) cdmol
    DO WHILE ((cdmol .lt. 1.e5).or.(cdmol.gt.1.e25))
      WRITE(*,*)'Please enter a value between 1e5 and 1e25'
      WRITE(*,21) 'Molecular column density [cm^-2] ?  '
      READ(*,*) cdmol
    END DO
    WRITE(13,23) cdmol

    WRITE(*,21) 'Line width [km/s] ?  '
    READ(*,*) deltav
    DO WHILE ((deltav.lt.1.e-3).or.(deltav.gt.1.e3))
      WRITE(*,*)'Please enter a value between 1e-3 and 1e3'
      WRITE(*,21) 'Line width [km/s] ?  '
      READ(*,*) deltav
    END DO
    WRITE(13,23) deltav
    !convert to cm/s
    deltav = deltav * 1.0e5
    
    WRITE(*,21) 'Choose a method to calculate Escape Probability!  '
    READ(*,*) method
    DO WHILE ((method .lt. 1) .or. (method .gt. 3))
      WRITE(*,*)'Please enter a number between 1 and 3 '
      WRITE(*,21) 'Which method 1=uniform sphere,2=LVG, 3=slab ?'
      READ(*,*) method
    END DO
    WRITE(13,25) method

    !Choose a method to get the first iteration
    WRITE(*,21) 'Choose a method to solve the system of equations!  '
    READ(*,*) imethod
    DO WHILE ((imethod < 1) .or. (imethod > 3))
      WRITE(*,*)'Please enter a number between 1 and 3'
      WRITE(*,21) 'Which method 1=Tbg, 2=Tkin, 3=mixed Tbg and Tkin?  '
      READ(*,*) imethod
    END DO
    WRITE(13,25) imethod
    write(*,*)' '

    WRITE(*,21) 'tolerance ?  '
    READ(*,*) tol
    DO WHILE ((tol.lt.1.e-6).or.(tol.gt.1.e-3))
      WRITE(*,*)'Please enter a value between 1e-6 and 1e-3'
      WRITE(*,21) 'tolerance error ?  '
      READ(*,*) tol
    END DO
    WRITE(13,23) tol

    WRITE(*,21) 'under-relaxation factor?  '
    READ(*,*) alpha
    DO WHILE ((alpha .lt. 0.) .or. (alpha .gt. 1.))
      WRITE(*,*)'Please enter a value between 0 and 1'
      WRITE(*,21) 'under-relaxation factor ?  '
      READ(*,*) alpha
    END DO
    WRITE(13,23) alpha

    Return
    97   WRITE(*,*) 'Error opening log file'
  END SUBROUTINE GetInputs


  SUBROUTINE Defaults
    !Set physical parameters to default values
    INTEGER :: ipart  ! to loop over collision partners

    tkin   = 30.0
    tbg    = 2.73
    cdmol  = 1.0e13
    deltav = 1.0
    alpha = 0.25

    density(1) = 1.0e5
    DO ipart=2,maxpart
      density(ipart) = 0.0
    END DO
  END SUBROUTINE Defaults

  FUNCTION length(str)
    !Returns the lengths of a string
    INTEGER :: length,i
    CHARACTER(*) :: str
    DO i=1,len(str)
       IF (str(i:i).eq.' ') THEN
          length=i-1
          RETURN
       END IF
    END DO
    length=len(str)
  END FUNCTION length

  SUBROUTINE output(tex, niter)
    !Writes results to file
    REAL(dp) :: tex(maxline) !! line excitation temperature
    INTEGER  :: niter        !! nr iterations

    integer ::iline    ! final number of iterations

    !Start with summary of input parameters
    OPEN(unit=8,file=outfile,status='unknown',err=98)
    30   format (a,f8.3)
    31   format (a,1pe10.3)
    32   format (a)

    write(8,32) '* Radex version        : '&
    &      //version(1:length(version))
    IF (method.eq.1)&
    &write(8,32) '* Geometry             : Uniform sphere'
    IF (method.eq.2)&
    &write(8,32) '* Geometry             : Expanding sphere'
    IF (method.eq.3) &
    &write(8,32) '* Geometry             : Plane parallel slab'
    write(8,32) '* Molecular data file  : '//moldat(1:80)
    write(8,30) '* T(kin)            [K]: ',tkin
    IF (density(1).gt.eps)&
      &write(8,31) '* Density of H2  [cm-3]: ',density(1)
    IF (density(2).gt.eps)&
      &write(8,31) '* Density of pH2 [cm-3]: ',density(2)
    IF (density(3).gt.eps)&
      &write(8,31) '* Density of oH2 [cm-3]: ',density(3)
    IF (density(4).gt.eps)&
      &write(8,31) '* Density of e-  [cm-3]: ',density(4)
    IF (density(5).gt.eps)&
      &write(8,31) '* Density of H   [cm-3]: ',density(5)
    IF (density(6).gt.eps)&
      &write(8,31) '* Density of He  [cm-3]: ',density(6)
    IF (density(7).gt.eps)&
      &write(8,31) '* Density of H+  [cm-3]: ',density(7)
    write(8,30) '* T(background)     [K]: ',tbg
    write(8,31) '* Column density [cm-2]: ',cdmol
    write(8,30) '* Line width     [km/s]: ',deltav/1.0d5

    write(8,33) '* Method to calculate the rate populations : '
    write(8,33) '  The escape Probability : ', method
    write(8,33) '  The initial matrix     : ', imethod
    write(8,31) '  The relative error-eps : ', tol
    write(8,30) '  The stab. factor       : ', alpha

    write(8,'(a,i5,a,1pe10.1)') '* Calculation finished in ',niter,' iterations with relative tolerance of ', tol
    33 format(a,i4,a)

    !Column header
    if (debug) then
    write(8,*)' '
    write(8,*) 'LINE,E_UP (K),FREQ (GHz),WAVEL (um),T_EX (K)'//&
      &',TAU,T_R (K), POP UP,POP LOW, FLUX (K*km/s), FLUX (erg/cm2/s)'


    DO iline=1,nline
      !Check if line within output freq range
      IF (spfreq(iline).lt.fmax.and.spfreq(iline).gt.fmin) THEN
        IF (dabs((tex(iline))).lt.1000.0) THEN
          write(8,113) upperQNum(iline),lowQNum(iline),eup(iline),spfreq(iline),wavelength(iline),&
            &tex(iline),taul(iline),antennaTemp(iline),upperPops(iline),lowerPops(iline),&
            &intensityKkms(iline),intensityErgs(iline)
        else
          write(8,114) upperQNum(iline),lowQNum(iline),eup(iline),spfreq(iline),wavelength(iline),&
            &tex(iline),taul(iline),antennaTemp(iline),upperPops(iline),lowerPops(iline),&
            &intensityKkms(iline),intensityErgs(iline)
        END IF
      END IF
    END DO
  else
 !c     Column header
      write(8,115)'     LINE         E_UP       FREQ        WAVEL     T_EX  TAU   &
      &     T_R       POP        POP      FLUX        FLUX'
!c     Units
      write(8,115)'                  (K)        (GHz)       (um)      (K)     (K)        UP        LOW     (K*km/s) (erg/cm2/s)'
115      format(11(1x, A))

    DO iline=1,nline
      !Check if line within output freq range
      IF (spfreq(iline).lt.fmax.and.spfreq(iline).gt.fmin) THEN
        IF (dabs((tex(iline))) .lt. 1000.0) THEN
          write(8,113) upperQNum(iline),lowQNum(iline),eup(iline),spfreq(iline),wavelength(iline),&
            &tex(iline),taul(iline),antennaTemp(iline),upperPops(iline),lowerPops(iline),&
            &intensityKkms(iline),intensityErgs(iline)
        else
          write(8,114) upperQNum(iline),lowQNum(iline),eup(iline),spfreq(iline),wavelength(iline),&
            &tex(iline),taul(iline),antennaTemp(iline),upperPops(iline),lowerPops(iline),&
            &intensityKkms(iline),intensityErgs(iline)
        END IF
      END IF
    END DO
  end if
113      format(a,' -- ', a, f8.1, 2(2x, f10.4), 1x, f8.3, 6(1x, 1pe10.3))
114      format(a,' -- ', a, f8.1, 2(2x, f10.4), 7(1x, 1pe10.3))


    RETURN
    98   write(*,*) 'error opening output file'
  END SUBROUTINE Output

  SUBROUTINE ParseInputDictionary(dictionary,status)
    INTEGER, INTENT(INOUT) :: status
    CHARACTER(LEN=*) :: dictionary
    INTEGER :: ipart          !! loop over collision partners
    CHARACTER(10) :: partner  !! name of collision partner
    INTEGER :: id             !! ID code of collision partner
    LOGICAL :: unreadParameters

    INTEGER :: posStart,posEnd !variables to help parse dictionary
    CHARACTER(LEN=200) :: inputParameter, inputValue !temporary storage for inputs
    unreadParameters=.True.
    !Set input parameters to default values
    CALL defaults
    posStart = scan(dictionary, '{')
    DO WHILE (unreadParameters)
      posEnd = scan(dictionary, ':')
      inputParameter = dictionary(posStart+2:posEnd-2)
      dictionary = dictionary(posEnd:)
      posStart = scan(dictionary, ' ')
      IF (scan(dictionary, ',') .EQ. 0) THEN
          posEnd = scan(dictionary, '}')
          unreadParameters=.False.
      ELSE
          posEnd = scan(dictionary, ',')
      END IF
      inputValue = dictionary(posStart+1:posEnd-1)
      dictionary = dictionary(posEnd:)
      SELECT CASE (inputParameter)
        CASE('molfile')
          READ(inputValue,*) molfile
          IF ((molfile(1:1).ne.'/').and.(molfile(1:1).ne.'.'))&
            &     molfile = radat(1:length(radat))//molfile(1:length(molfile))
        CASE('fmin')
          READ(inputValue,*) fmin
        CASE('fmax')
          READ(inputValue,*) fmax
        CASE('tkin')
          READ(inputValue,*) tkin
          IF ((tkin.lt.0.1).or.(tkin.gt.1.e4)) THEN
            WRITE(*,*) 'Please enter a TKin between 0.1 and 1e4'
            status=0
            RETURN
          END IF
        CASE("tbg")
          READ(inputValue,*) tbg
          IF ((tbg.lt.-1.e4).or.(tbg.gt.1.e4)) THEN
            WRITE(*,*)'Please enter a value between -1e4 and 1e4'
            status=0
            RETURN
          END IF
        CASE("linewidth")
          READ(inputValue,*) deltav
          IF ((deltav.lt.1.e-3).or.(deltav.gt.1.e3)) THEN
            WRITE(*,*)'Please enter a value between 1e-3 and 1e3'
            status=0
            RETURN
          END IF
          deltav = deltav * 1.0e5
        CASE("cdmol")
          READ(inputValue,*) cdmol
          IF ((cdmol.lt.1.e5).or.(cdmol.gt.1.e25)) THEN
            WRITE(*,*)'Please enter a column density between 1e5 and 1e25'
            status=0
            RETURN
          END IF
        CASE('h2')
          READ(inputValue,*) density(1)
        CASE('p-h2')
          READ(inputValue,*) density(2)
       CASE('o-h2')
          READ(inputValue,*) density(3)
        CASE('e-')
          READ(inputValue,*) density(4)
        CASE('h')
          READ(inputValue,*) density(5)
        CASE('he')
          READ(inputValue,*) density(6)
        CASE('h+')
          READ(inputValue,*) density(7)
        CASE('geometry')
          READ(inputValue,*) method
        CASE('Initial')
          READ(inputValue,*) imethod
        CASE('Tolerancy')
          READ(inputValue,*) tol
          if ((tol .lt. 1.e-6) .or. (tol .gt. 1.e-3)) then
            WRITE(*,*)'Please enter a tolerance between 1e-6 and 1e-3'
            status=0
            RETURN
          END IF
        CASE('Factor')
          READ(inputValue,*) alpha
          IF ((alpha .lt. 0.) .or. (alpha .gt. 1.)) then
            WRITE(*,*)'Please enter a value between 0 and 1'
            status=0
            RETURN
          END IF
        CASE DEFAULT
          write(*,*) "Unknown parameter in input dictionary, acceptable parameters are:"
          write(*,*) "molfile, tkin, tbg, cdmol, linewidth, h2, h, e-, p-h2, o-h2, h+"
          write(*,*) "Input dictionary was:"
          write(*,*) dictionary

      END SELECT
    END DO
    !Default values: DC / Lyman limit
    IF (fmin.eq.fmax) THEN
        write(*,*) "min and max frequencies equal, using defaults"
       fmin = 0.d0
       fmax = 3.d7    
    END IF
    IF (fmin.gt.fmax) THEN
      write(*,*) "fmin greater than fmax, switching them over and continuing"
       fmin = fmin + fmax
       fmax = fmin - fmax
       fmin = fmin - fmax
    END IF
    fmin = dmax1(0.d0,fmin)
    fmax = dmin1(3.d7,fmax)

  END SUBROUTINE ParseInputDictionary
END MODULE IO
