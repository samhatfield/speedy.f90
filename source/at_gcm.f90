program agcm_main
    implicit none

    ! program : agcm_main

    ! experiment identifier
    character(len=3) :: cexp = 'exp'
    integer :: jday, ndays

    ! 1. initialization
    ! ndays = no. of integration days, set by agcm_init
 
    call agcm_init (cexp, 0, 0, 0, ndays)

    print *, 'integration length in days: ', ndays

    ! 2. do loop over total no. of integration days
    do jday = 1, ndays

        ! 2.2 run atmospheric model for 1 day
        call agcm_1day(jday, cexp)

        ! 2.1 exchange data with coupler
        call agcm_to_coupler(jday)
        call coupler_to_agcm(jday)

        ! 2.3 write restart file at the end of selected months and at the end of
        ! the integration 
        call restart(jday)
    enddo
end

subroutine agcm_1day(jday, cexp)
    ! subroutine agcm_1day (jday)
    !
    ! perform atm. model integration for 1 day, 
    ! post-proc. and i/o at selected times 

    implicit none

    include "atparam.h"
    include "atparam1.h"

    include "com_tsteps.h"
    include "com_date.h"

    integer, intent(in) :: jday
    character(len=3), intent(in) :: cexp
    integer :: istep

    if (iday == 1) print *, ' start of year/month = ', iyear, imonth

    istep = 1 + (jday - 1) * nsteps

    ! 1. set forcing terms according to date
    call fordate(1)

    ! 2. set daily-average flux arrays to zero
    call dmflux(0)

    ! 3. integrate the atmospheric model for 1 day
    call stloop(istep)

    ! 4. write daily-mean output
    call dmout(idout)

    ! 5. compute new date
    call newdate(1)

    ! 6. write time-mean output files and restart file at the end of selected
    ! months
    if (iday == 1) then
        ! write monthly-mean output for previous month
        if (nstout < 0) call tmout(1)
        
        ! open new output files at the beginning of each year
        if (imonth == 1 .and. jday < ndaytot) call setgrd(1, cexp)
    endif
    return
end
