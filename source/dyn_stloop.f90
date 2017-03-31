subroutine stloop(istep)
    !   subroutine stloop (istep)
    !
    !   Purpose: Perform a series of time steps calling 
    !            post-processing/output routines at selected steps
    !   Input/output : istep = time step index
    !   Updated common block : lflag2
      
    use mod_lflags, only: lradsw, lrandf
    use mod_tsteps
    use mod_date, only: ihour, newdate
    use mod_dynvar

    implicit none

    integer, intent(inout) :: istep
    integer :: iitest = 0, j, jj

    ! Break up each day into four 6-hour windows
    do jj = 1, 4
        ! Each 6-hour window has nsteps/4 actual timesteps
        do j = 1, nsteps/4
            if (iitest == 1) print*, 'stloop: calling step ', istep
    
            ! Set logical flags
            lradsw = (mod(istep,nstrad) == 1)
            lrandf = ((istep <= nstrdf) .or. (nstrdf < 0))
    
            ! Perform one leapfrog time step
            call step(2, 2, delt2, alph, rob, wil)   
    
            ! Do diagnostic, post-processing and I/O tasks 
            call diagns(2, istep)
    
            if (ihout .eqv. .false.) then
                if (mod(istep, nstppr) == 0) call tminc
                if (nstout > 0 .and. mod(istep, nstout) == 0) call tmout(1)
            end if
    
            istep = istep + 1
        end do

        ! Increment hour timer (takes values of 0, 6, 12 or 18)
        ihour = mod(ihour + 6, 24)

        ! If it's a new day...
        if (ihour .eq. 0) then
            ! Compute new date
            call newdate(1)
        end if

        if (ihout) then
            if (ipout) call iogrid (2) !output for every 6 hours
            call iogrid (4) !gridded data output for every 6 hours
        end if

        if (sixhrrun) then
            call restart (2)
            print *,'normal end with 6-hr fcst (yeahhhhhhh!!!!)'
            stop 1111
        end if
    end do
end
