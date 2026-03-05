!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! RADEX
!! Original code by Van der Tak et al. 2007
!! All publications using this code should reference their release paper: A&A 468, 627 (2007)
!!
!!   This amended version by Jon Holdship has been updated to Modern Fortran
!!   with a view to removing common blocks and making an F2PY module for python.
!!
!!  This code has been tested against the original but the authors do not
!!  guarantee there are no errors. We advise users to check results they intend
!!  to publish with another version of RADEX or another code.
!!
!!-Version code : 2024-The code has been modified by Christina Dwiriyanti based on the SpectralRadex code 
!! by Jonathan Holdship : https://github.com/uclchem/SpectralRadex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!Main program: controls program flow and drives subroutines
PROGRAM RADEX
   USE IO
   USE Solver, only :  solve_matrix, CalcOutputArrays
   USE Background

   IMPLICIT NONE
    INTEGER  :: niter, nlines   !! iteration counter
    INTEGER  :: imore = 1       !! are we running again?
    LOGICAL  :: conv            !! are we converged?
    integer  :: status          !! status error
    integer  :: iline           !! index iline
    REAL(dp) :: xb(maxlev)      !! RHS of rate equation-classical method
    REAL(dp) :: xt(maxlev)      !! RHS of rate equation-modern method
    REAL(dp) :: drhs(maxlev)    !! RHS of rate equation-mixed method
    REAL(dp) :: btex(maxline)   !! temporary line excitation temperature using the classic method
    REAL(dp) :: itex(maxline)   !! temporary line excitation temperature using the modern method
    REAL(dp) ::  tex(maxline)   !! line excitation temperature (solution)
    REAL(dp) :: dtex, tsum, tmp    !! diff
    REAL(dp) :: xpb(maxlev)     !! xpop  : level populations
    REAL(dp) :: xpk(maxlev)     !! xpop  : level populations
    REAL(dp) :: xpop(maxlev)                !! xpop  : level populations
    REAL(dp) :: yrate(maxlev,maxlev)  !! rate matrix
    REAL(dp) :: brate(maxlev,maxlev)  !! rate matrix


    character(len=10) :: imode
    !!temp
    integer  :: ii, n0, n1, n2, n3, n4, n5, nit, jj
    ! Begin executable statements
    write(*,*)
    write(*,*)'   Welcome to Radex, Moden Fortran Edition 2024'
    write(*,*)
    ii = 1

    DO WHILE (imore .eq. 1)
        !     Get input parameters
        IF (DEBUG) write(*,*) 'calling getinputs'

        CALL getinputs(status)

        ! Read data file
        IF (DEBUG) write(*,*) 'calling readdata'
        CALL ReadData(status)

        write(*,*) " Beginning Calculation", imethod

        ! Calculate background radiation field
        tmp = tbg

        ! Calculate background radiation field
        IF (DEBUG) write(*,*) 'calling backrad'
        CALL backrad(status)

        niter = 0
        conv  = .false.

        ! Set up rate matrix, splitting it in radiative and collisional parts
        !Invert rate matrix to get `thin' starting condition
        !!Choose the method to calculate the excitation temperature
        IF (DEBUG) write(*,*) 'calling matrix'
        ii = ii + 1
        n0 = 15+ii
        n1 = 200+ii
        n2 = 300+ii
        n3 = 350+ii
        n4 = 400+ii
        n5 = 124
        open(n1,FILE='tbg1.out')
        open(n2,FILE='tkin1.out')
        open(n3,FILE='mix1.out')
        open(n4,FILE='mixoutput.dat')
        open(n5,FILE='mix_iter1.dat')

        if (imethod == 1) then
            imode = 't_bg'
            call solve_matrix(niter, tbg, xb,  btex, conv,  status)
            do iline = 1, nline
               write(n1,FMT='(I3,X, 4(ES13.6,X), I3, X)')iline, xb(iline), btex(iline), xpop(iline), taul(iline),niter
            enddo
        elseif (imethod == 2) then
            imode = 't_kin'
            call solve_matrix(niter, tkin, xt, itex, conv, status)
            do iline = 1, nline
              write(n2,FMT='(I3,X, 4(ES13.6,X), I3,X )')iline, xt(iline), itex(iline), xpop(iline), taul(iline),niter
            enddo
        else
            imode = 't_bg&t_kin'
            call solve_matrix(niter, tbg,  xb, btex, conv, status)
            call solve_matrix(niter, tkin, xt, itex, conv, status)
            tsum = 0.0
            do iline = 1, nline !n3=350+ii; n4=400+ii
                 tex(iline)  = 0.5*(btex(iline) + itex(iline))
                       dtex  = abs(( 1.0 - btex(iline)/itex(iline)))
                       tsum  = tsum  + abs(dtex)
                write(n3,FMT='(I3,X, 5(ES13.6,X), I3, X)')iline, xb(iline), btex(iline), itex(iline), tex(iline), taul(iline), niter
                write(n4,FMT='(I3,X, 8(ES13.6,X))')iline, xb(iline), xt(iline), btex(iline), itex(iline), &
                 taul(iline), tex(iline)
            end do
            !if (tsum .lt. ccrit .or. dsum .lt. ccrit) conv = .true.
            if ( tsum/nthick .lt. ccrit) conv = .true.
            write(n5,FMT='(I3,X,ES13.6)')niter, tsum
        endif

        ! Start iterating
        DO niter = 1, maxiter
           !Invert rate solve_matrix using escape probability for line photons
           if (imethod == 1) then
              !! Classical method
              call solve_matrix(niter, tbg,  xb, tex, conv, status)

              do iline = 1, nline
                write(n1,FMT='(I3,X, 4(ES13.6,X), I3,X )')iline, xb(iline), tex(iline), xpop(iline), taul(iline),niter
              enddo
           else if (imethod == 2) then
              !! Modern method
              call solve_matrix(niter, tkin, xt, tex, conv, status)

              do iline = 1, nline
                write(n2,FMT='(I3,X, 4(ES13.6,X), I3,X )')iline, xt(iline), tex(iline), xpop(iline), taul(iline),niter
              enddo
           else
              !! Mixed method
              call solve_matrix(niter, tbg,  xb, btex, conv, status)
              call solve_matrix(niter, tkin, xt, itex, conv, status)
 
              tsum = 0.0
              do iline = 1, nline
                 if (taul(iline) .gt. 0.01) then
                  dtex  = abs(( 1.0 - itex(iline)/btex(iline)))
                  tsum  = tsum  + dtex
                 endif
                 tex(iline) = 0.5*(btex(iline) + itex(iline))
                 write(n3,FMT='(I3,X, 5(ES13.6,X), I3,X )')iline, xb(iline), btex(iline), &
                 itex(iline), tex(iline), taul(iline), niter
                 write(n4,FMT='(I3,X, 8(ES13.6,X))')iline, xb(iline), xt(iline), btex(iline), itex(iline), &
                                  taul(iline), tex(iline)
              end do
               ! if ( tsum/nthick .lt. ccrit .or. dsum/nthick .lt. ccrit) conv = .true.
               if ( tsum/nthick .lt. ccrit) conv = .true.
               write(n5,FMT='(I3,X,ES13.6)')niter, tsum/nthick
           end if
           if (conv) then
              write(*,'(A,X,I4,X,2A)') 'Finished in ',niter,' iterations using the initial value: ', imode
              exit
           endif
        END DO

        IF (.NOT. conv) write(*,*) ' Warning: Calculation did not converge in ', niter, maxiter, ' iterations.'

        !! Prepare output by calculating final quantities
        IF (DEBUG) write(*,*) 'calculating output summary variables'
        CALL CalcOutputArrays(tex, xt, nlines)

        !! Write output
        IF (DEBUG) write(*,*) 'calling output'
        call output(tex, niter)

        !! See if user wants more, else call it a day
        write(*,'(A)') '  Another calculation [0/1] ? '
        read(*,*) imore
        write(13, '(A,I5)') imore
    END DO

    write(*,*) '   Have a nice day.'

    ! Done! Now close log file.
    close(13)
    ! ...and output file.
    close(8)
    close(129)
END PROGRAM RADEX
