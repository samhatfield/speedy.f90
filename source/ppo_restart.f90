subroutine restart(jday)
    !  subroutine restart (jday)
    !
    !  Purpose : read or write a restart file
    !  Input :   JDAY  = 0 : read model variables from a restart file
    !                  > 0 : write model variables  to a restart file
    !                        at selected dates and at the end of run 
    !
    
    use mod_tsteps, only: nmonrs, iyear0, imont0
    use mod_atparam
    use mod_dynvar
    use mod_date, only: iyear, imonth, iday, ndaytot, ihour

    implicit none

    integer, intent(in) :: jday
    integer :: jrec
    real :: adummy

    if (jday.eq.0) then
        100 CONTINUE

        ! 1. Read the restart dataset corresponding to the specified initial date
        ! [Modified:] Read the restart dataset for any initial date
!       read (3,end=200) iyear, imonth
        read (3,end=200) iyear, imonth, iday, ihour

!        if (iyear.eq.iyear0.and.imonth.eq.imont0) then
!           print*, 'read restart dataset for year/month: ', iyear,imonth
            print '(A,I4.4,A,I2.2,A,I2.2,A,I2.2)',&
                & 'Read restart dataset for year/month/date/hour: ', &
                & iyear,'/',imonth,'/',iday,'/',ihour
            
            read (3) vor
            read (3) div
            read (3) t
            read (3) ps
            read (3) tr

            call rest_land(0)
            call rest_sea(0)
!        else
!            print *, 'Skip restart dataset for year/month: ', iyear,imonth
!            
!            do jrec=1,5
!              read (3) adummy
!            end do
!  
!            CALL REST_LAND(0)
!            CALL REST_SEA(0)
!  
!            go to 100
!        end if
    ! Check for write-up dates
!    else if ( (iday.eq.1) .and.&
!        & (mod(imonth-1,nmonrs).eq.0.or.jday.eq.ndaytot) ) then
    else
        ! 2. Write date and model variables to the restart file
!        print*, 'Write restart dataset for year/month: ', IYEAR,IMONTH
         print*, 'Write restart dataset for year/month/date/hour: ', &
             & iyear,'/',imonth,'/',iday,'/',ihour

!        write (10) iyear, imonth
         write (10) iyear, imonth, iday, ihour

         write (10) vor
         write (10) div
         write (10) t
         write (10) ps
         write (10) tr

         call rest_land(1)
         call rest_sea(1)
    end if

    return

    ! 4. Stop integration if restart file is not found
    200 continue

    print*, ' No restart dataset for the specified initial date'

    stop 'invalid restart'
end
