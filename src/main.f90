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
!!- Modern version code is based on the SpectralRadex code by Jonathan Holdship : https://github.com/uclchem/SpectralRadex
!!- The code has been modified by Christina Dwiriyanti (sron)
!!- New functionalities has been added 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!Main program: controls program flow and drives subroutines
PROGRAM RADEX
   USE IO
   USE Solver
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
    REAL(dp) :: tex(maxline)    !! line excitation temperature (solution)
    REAL(dp) :: dtex, tsum, dsum  !! diff

    ! Begin executable statements
    write(*,*)
    write(*,*)'   Welcome to Radex, Moden Fortran Edition 2024'
    write(*,*)

    DO WHILE (imore .eq. 1)
        !     Get input parameters
        IF (DEBUG) write(*,*) 'calling getinputs'

        CALL getinputs(status)

        !     Read data file
        IF (DEBUG) write(*,*) 'calling readdata'
        CALL ReadData(status)

        write(*,*) " Beginning Calculation", imethod, method

        !     Calculate background radiation field
        IF (DEBUG) write(*,*) 'calling backrad'
        CALL backrad(status)

        niter = 0
        conv  = .false.

        !Set up rate matrix, splitting it in radiative and collisional parts
        !Invert rate matrix to get `thin' starting condition
        !!Choose the method to calculate the excitation temperature
        IF (DEBUG) write(*,*) 'calling matrix'

        if (imethod == 1) then
            call matrix(niter, tbg, xb,  tex, conv, status)
        elseif (imethod == 2) then
            call matrix(niter, tkin, xt, tex, conv, status)
        else
            call matrix(niter, tbg,  xb, btex, conv, status)
            call matrix(niter, tkin, xt, itex, conv, status)
            tsum = 0.0
            do iline = 1, nline
                 drhs(iline) = abs(xb(iline) - xt(iline))
                       dtex  = abs(btex(iline) - itex(iline))
                 tex(iline)  = 0.5*(btex(iline) + itex(iline))
                       tsum  = tsum  + abs(drhs(iline))
                if(debug) call output_dummy(iline, xb(iline), xt(iline), btex(iline), itex(iline), spfreq(iline))
            end do
            if (tsum .lt. ccrit) conv = .true.
        endif

        ! Start iterating
        DO niter = 1, maxiter
           !Invert rate matrix using escape probability for line photons
           if (imethod == 1) then
              !! Classical method
              call matrix(niter, tbg,  xb, tex, conv, status)              
           else if (imethod == 2) then
              !! Modern method
              call matrix(niter, tkin, xt, tex, conv, status)
           else
              !! Mixed method
              call matrix(niter, tbg,  xb, btex, conv, status)
              call matrix(niter, tkin, xt, itex, conv, status)
              tsum = 0.0
              dsum = 0.0
              do iline = 1, nline
                   drhs(iline) = abs(xb(iline) - xt(iline))
                    dsum  = dsum  + abs(drhs(iline))
                    if (taul(iline) .gt. 0.01) then
                       dtex  = abs((btex(iline) - itex(iline))/btex(iline))
                    endif
                    tex(iline) = 0.5*(btex(iline) + itex(iline))
                         tsum  = tsum  + dtex
                   if(debug) call output_dummy(iline, xb(iline), xt(iline), btex(iline), itex(iline), spfreq(iline))
               end do
               if ( tsum/nthick .lt. ccrit .or. dsum .lt. ccrit) conv = .true.
           end if
           if (conv) then
              write(*,*) 'Finished in ',niter,' iterations.'
              exit
           endif
        END DO

        IF (.NOT. conv) write(*,*) '   Warning: Calculation did not converge in ', niter, maxiter, ' iterations.'

        !! Prepare output by calculating final quantities
        IF (DEBUG) write(*,*) 'calculating output summary variables'
        CALL CalcOutputArrays(tex, nlines)

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
