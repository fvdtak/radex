MODULE Solver
  USE CommonParams
  USE Slatec
  USE types

IMPLICIT NONE
Contains
  SUBROUTINE solve_matrix(niter, ival, xpb, tex, conv, status)
    !! Set up rate matrix
    INTEGER  :: niter                               !! Iteration counter
    REAL(dp)  :: rhs(maxlev)           !! RHS of rate equation
    REAL(dp), intent(in)   :: ival                  !! Initial value of the matrix niter=0
    REAL(dp), intent(out)  :: tex(maxline)          !! Line exitation temperature (solution)
    LOGICAL, intent(inout) :: conv                  !! are we converged?
    REAL(dp) :: xpop(maxlev)           !! xpop  : level populations

   ! REAL(dp), intent(out) :: xpop(maxlev)           !! xpop  : level populations
    REAL(dp), intent(out) :: xpb(maxlev)            !! xpop  : level populations
    REAL(dp) :: yrate(maxlev,maxlev)   !! rate matrix
    INTEGER, INTENT(INOUT) :: status                !! status error
    !! Local variables
    INTEGER  :: ilev,jlev,klev        !! to loop over energy levels
    INTEGER  :: nplus                 !! to solve statistical equilibrium
    INTEGER  :: iline                 !! to loop over lines
    INTEGER  :: m,n                   !! line upper/lower levels
    INTEGER  :: nfat                  !! counts highly optically thick lines
    INTEGER  :: nreduce               !! size of reduced rate matrix
    INTEGER  :: indx(maxlev),dsign    !! needed for NumRep equation solver
    INTEGER  :: terminate             !! terminate the procedure
    REAL(dp) :: etr,exr               !! to calculate radiative rates
    REAL(dp) :: xt                    !! frequency cubed
    REAL(dp) :: hnu                   !! photon energy
    REAL(dp) :: bnutex                !! line source function
    REAL(dp) :: cddv                  !! N(mol) / delta V
    REAL(dp) :: beta                  !! escape probability
    REAL(dp) :: bnu                   !! Planck function
    REAL(dp) :: uarray(maxlev,maxlev) !! reduced rate matrix
    REAL(dp) :: redcrit               !! reduction criterion
    REAL(dp) :: sumx                  !! summed radiative rate
    REAL(dp) :: total                 !! to normalize populations
    REAL(dp) :: tsum,thistex, dsum    !! to check convergence

    !! keep old population for underrelaxation procedure -- sb/fvdt 30nov2011
    REAL(dp) :: xpopold(maxlev), ltex
    !! an option to reduce the size of the matrix db
    LOGICAL  :: reduce
    integer :: n6, n7
    reduce = .false.
    if (imethod /=3) then
      n6 = 126
      n7 = 125
      open(n6,FILE='tbg_iter1.dat')
      open(n7,FILE='tkin_iter1.dat')
    end if
    !! default value of the relative error
    if (tol < 0) ccrit = 1.0e-6
    ccrit = tol

    !! Executable statements begin here
    IF (debug) write(*,*) 'niter = ', niter

    !! Clear array of level populations.
    DO ilev=1, nlev
      rhs(ilev) = 0.0
      DO jlev=1, nlev
        yrate(ilev,jlev) = 0.0
      END DO
    END DO

    !! Initialize rate matrix
    nplus = nlev + 1
    DO ilev=1, nlev
       DO jlev=1, nlev
         yrate(ilev,jlev) = -1.0d-30*totdens
       END DO
       !! Add conservation equation
       yrate(nplus,ilev) = 1.0d0
       rhs(ilev)         = 1.0e-30*totdens
       yrate(ilev,nplus) = 1.0d-30*totdens
    END DO

    !! rhs for jlev=nplus
    rhs(nplus) = 1.0e-30*totdens

    !! Contribution of radiative processes to the rate matrix.
    !! First iteration: use background intensity
    IF (niter .eq. 0) THEN
       DO iline = 1, nline
          if (ival .gt. 0.0d0) trj(iline) = ival
          m   = iupp(iline)
          n   = ilow(iline)
          etr = fk*xnu(iline)/trj(iline)
          if (etr .ge. 160.0d0) then
             exr = 0.0d0
          else
             exr = 1.0/(dexp(etr)-1.0d0)
          endif
          yrate(m,m) = yrate(m,m) + aeinst(iline)*(1.0 + exr)
          yrate(n,n) = yrate(n,n) + aeinst(iline)*(gstat(m)/gstat(n))*exr
          yrate(m,n) = yrate(m,n) - aeinst(iline)*(gstat(m)/gstat(n))*exr
          yrate(n,m) = yrate(n,m) - aeinst(iline)*(1.0 + exr)
          if (debug) write(135,*)n, m, yrate(m,m), yrate(m,n), yrate(n,m), yrate(n,n)
       END DO
    else
       !! Subsequent iterations: use escape probability.
       cddv = cdmol / deltav
       !! Count optically thick lines
       nthick = 0
       nfat   = 0
       DO iline = 1, nline
          xt  = xnu(iline)**3.0
          m   = iupp(iline)
          n   = ilow(iline)

          !! Calculate source function
          hnu = fk * xnu(iline) / tex(iline)
          if( debug .and. niter == 1 .and. imethod == 3) write(188,*) iline, tex(iline), hnu
          if(hnu .ge. 160.0) then
             bnutex = 0.0d0
          else
             bnutex = thc*xt/(dexp(fk*xnu(iline)/tex(iline))-1.0)
          endif

          !! Calculate line optical depth.
          taul(iline) = cddv*(xpop(n)*gstat(m)/gstat(n)-xpop(m)) &
               &      / (fgaus*xt/aeinst(iline))
          if(taul(iline) .gt. 1.d-2) nthick = nthick+1
          if(taul(iline) .gt. 1.d05) nfat   = nfat+1

          !! Use escape probability approximation for internal intensity.
          beta = EscProb(taul(iline))

          !! Split off local contribution to radiation field  sb/fvdt 30nov2011
          bnu  = totalb(iline)*beta

          exr  = bnu/(thc*xt)

          ! Radiative contribution to the rate matrix
          yrate(m,m) = yrate(m,m)+aeinst(iline)*(beta+exr)
          yrate(n,n) = yrate(n,n)+aeinst(iline)*(gstat(m)*exr/gstat(n))
          yrate(m,n) = yrate(m,n)-aeinst(iline)*(gstat(m)/gstat(n))*exr
          yrate(n,m) = yrate(n,m)-aeinst(iline)*(beta+exr)
       END DO
    END IF

    !!  Warn user if convergence problems expected
    IF ((niter.eq.1).and.(nfat.gt.0)) WRITE(*,*)&
         &"*** Warning: Some lines have very high optical depth"

    IF (debug) THEN
       WRITE(*,*) yrate(1,1),yrate(1,2),yrate(1,3),yrate(1,4)
       WRITE(*,*) yrate(2,1),yrate(2,2),yrate(2,3),yrate(2,4)
       WRITE(*,*) yrate(3,1),yrate(3,2),yrate(3,3),yrate(3,4)
       WRITE(*,*) yrate(4,1),yrate(4,2),yrate(4,3),yrate(4,4)
    END IF

    !! Contribution of collisional processes to the rate matrix.
    DO ilev=1,nlev
       yrate(ilev,ilev) = yrate(ilev,ilev) + ctot(ilev)
       DO jlev=1,nlev
          if(ilev .ne. jlev) yrate(ilev,jlev) = yrate(ilev,jlev) - crate(jlev,ilev)
       END DO
    END DO

    IF (debug) THEN
       WRITE(*,*) yrate(1,1),yrate(1,2),yrate(1,3),yrate(1,4)
       WRITE(*,*) yrate(2,1),yrate(2,2),yrate(2,3),yrate(2,4)
       WRITE(*,*) yrate(3,1),yrate(3,2),yrate(3,3),yrate(3,4)
       WRITE(*,*) yrate(4,1),yrate(4,2),yrate(4,3),yrate(4,4)
    END IF

    !! db
       IF (debug) WRITE(*,*) 'inverting non-reduced matrix...'
       if (niter .ge. 0) call lubksb(yrate, nplus, maxlev, rhs, status)
       IF (status .eq. 0) RETURN
    
    !! Level populations are the normalized RHS components
    total = 0.0d0
    DO ilev = 1, nlev
       total = rhs(ilev)+total
    END DO

    !! Debugging
    IF (debug) WRITE(*,*) 'total rhs=',total
    IF (debug) WRITE(*,*) 'rhs=',(rhs(ilev),ilev=1,nlev)

    !! Limit population to minpop
    DO ilev=1,nlev
       xpopold(ilev) = xpop(ilev) ! xpop_tb and xpop_tk (separate)
       xpop(ilev)    = dmax1(minpop,rhs(ilev)/total)
       if (niter == 0) xpopold(ilev) = xpop(ilev)
       if (niter == 0 .and. debug) write(125,*)ilev, rhs(ilev), xpop(ilev)
    END DO

    IF (debug) WRITE(*,*) 'computing T_ex...'

    !! Compute excitation temperatures of the lines
    tsum = 0.0
    dsum = 0.0
    DO iline = 1, nline
      m  = iupp(iline)
      n  = ilow(iline)
      xt = xnu(iline)**3.d0
      if (niter .eq. 0) then
         if ((xpop(n) .le. minpop) .or. (xpop(m) .le. minpop)) then
            tex(iline) = ival!totalb(iline)
         else
                  ltex = (dlog(xpop(n)*gstat(m)/(xpop(m)*gstat(n))))
            tex(iline) = ival!fk*xnu(iline)/ltex
         end if
         dsum = abs(xpopold(iline)-xpop(iline))
      else
         if ((xpop(n) .le. minpop) .or. (xpop(m) .le. minpop)) then
            thistex = tex(iline)
         else
            thistex = fk*xnu(iline)/&
                 &(dlog(xpop(n)*gstat(m)/(xpop(m)*gstat(n))))
         end if

         !! Update excitation temperature & optical depth
         if (imethod == 3) then
            tex(iline) = thistex
         else
            !! Only thick lines count for convergence
            if (taul(iline) .gt. 0.01) then
               tsum = tsum + abs((thistex-tex(iline))/thistex)
            endif
            tex(iline)  = 0.5*(thistex + tex(iline))
            dsum = abs(xpopold(iline)-xpop(iline))
         endif

         !Calculated optical depth  at the center of the spectral line
         taul(iline) = cddv*(xpop(n)*gstat(m)/gstat(n)-xpop(m)) &
                     / (fgaus*xt/aeinst(iline))
      end if
    END DO

    xpop = alpha*xpop + (1.-alpha)*xpopold
    xpb = xpop

    !!  Introduce a minimum number of iterations
    if(imethod /= 3) then
      if(niter .ge. 1) THEN
        if(nthick .eq. 0) conv = .true.
        if(tsum/nthick .lt. ccrit) conv = .true.
        IF (.not.debug .and. imethod == 1) write(n6,FMT='(I3,X,ES13.6)')niter, tsum/nthick
        IF (.not.debug .and. imethod == 2) write(n7,FMT='(I3,X,ES13.6)')niter, tsum/nthick
      else
        IF (.not.debug .and. imethod == 1) write(n6,FMT='(I3,X,ES13.6)')niter, tbg
        IF (.not.debug .and. imethod == 2) write(n7,FMT='(I3,X,ES13.6)')niter, tkin
      END IF
    endif

    !! sb301111 now DO the under-relaxation!
    !! to increase the stability.
    !! An under-relaxation factor a such that a*x^i + (1-a)*x^{i-1} with 0 < a <= 1
    !! here we choose alpha = 0.25 (default).
    
    if(debug) then
      do iline = 1, nline
        if (imethod == 3) write(227,*) iline, niter, alpha, xpopold(iline), xpop(iline), taul(iline)
      enddo
    endif
    return
  END SUBROUTINE solve_matrix


  SUBROUTINE Matrix_reduce(yrate, rhs, status)
    !!    Reduce the size of matrix
    REAL(dp) :: rhs(maxlev)           !! RHS of rate equation
    REAL(dp) :: yrate(maxlev,maxlev)  !! rate matrix
    INTEGER  :: status                !! status

    ! local variables
    INTEGER  :: ilev,jlev,klev        ! to loop over energy levels
    INTEGER  :: nreduce               ! size of reduced rate matrix
    REAL(dp) :: uarray(maxlev,maxlev) ! reduced rate matrix
    REAL(dp) :: redcrit               ! reduction criterion
    REAL(dp) :: sumx                  ! summed radiative rate

      IF (debug) WRITE(*,*) 'reducing matrix...' 

      DO jlev=1,nlev
         DO ilev=1,nlev
            uarray(ilev,jlev) = yrate(ilev,jlev)
         END DO
      END DO

      ! Now test whether the matrix should be reduced to exclude the radiatively coupled levels.

      redcrit = 10.0*tkin/fk
      nreduce = 0
      DO ilev = 1,nlev
         if(eterm(ilev) .le. redcrit) nreduce = nreduce+1
      END DO

      IF (debug) WRITE(*,*) 'nreduce=',nreduce

      ! We now separate the collisionally coupled levels from those that
      ! are coupled mainly by radiative processes, compute an effective
      ! cascade matrix for rates of transfer from one low-lying level
      ! to another and then solve this reduced system of equations
      ! explicitly for the low-lying levels only.
      DO jlev = 1,nreduce
         DO ilev = 1,nreduce
            DO klev = nreduce+1, nlev
               uarray(ilev,jlev) = abs(yrate(klev,jlev)*yrate(ilev,klev)&
               &                 / yrate(klev,klev)) + uarray(ilev,jlev) 
            END DO
         END DO
      END DO

      ! Invert this reduced matrix
      IF (debug) WRITE(*,*) 'inverting reduced matrix...'

      ! solving the system of equations      
      call lubksb(uarray, nreduce+1, maxlev, rhs, status)

      IF (status .eq.0) RETURN

      IF (debug) WRITE(*,*) 'computing cascade...'

      ! Compute the populations of the highly excited states
      if (nlev .gt. nreduce) then
         do klev = nreduce+1, nlev
            sumx = 0.0
            do jlev = 1,nreduce
               sumx = rhs(jlev)*yrate(klev,jlev) + sumx
            enddo
            rhs(klev) = abs(sumx/yrate(klev,klev))
         enddo
      endif

    return
  END SUBROUTINE Matrix_reduce
END MODULE Solver      
