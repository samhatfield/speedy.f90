subroutine agcm_init(cexp, inidate, ntimes, irstart, ndays)
    !   subroutine agcm_init (cexp,inidate,ntimes,irstart,
    !  &                      ndays)
    !
    !   purpose: initialization of atmos. model and coupling interface 
    !
    !c      include "atparam.h"
    !c      include "atparam1.h"
    
    !c      parameter ( ngp = ix*il )

    use cpl_flags, only: icsea, isstan

    implicit none

    include "com_tsteps.h"
    include "com_date.h"

    include "com_lflags.h"

    ! input (reset by input/include files if inidate = 0):
    character(len=3), intent(inout) :: cexp        ! experiment identifier
    integer, intent(in) :: inidate     ! initial date yyyymm
    integer, intent(in) :: ntimes      ! integr. length in months (< 0) or days (> 0)
    integer, intent(in) :: irstart     ! restart flag: 0 = no, > 0 = yes

    ! output:
    integer, intent(inout) :: ndays       ! total no. of integration days

    print *, ' hallo from speedy_agcm'

    ! 1. set run initial time, duration, time-stepping and coupling options
    if (inidate <= 0) then
       read (2,*) istart
       read (2,'(a3)') cexp
    else
       istart = irstart
    endif

    if (istart /= 0) istart = 1

    include "cls_instep.h"

    if (inidate > 0) then
       iyear0 = inidate/100
       imont0 = mod(inidate,100)

       isst0 = (iyear0 - issty0) * 12 + imont0

       if (ntimes < 0) then
          nmonts = -ntimes
          ndaysl = 0
       else
          nmonts = 0
          ndaysl = ntimes
       endif
    endif

    call newdate(0)

    ndays = ndaytot

    ! check consistency of coupling and prescribed SST anomaly flags
    if (icsea >= 4) isstan = 1

    ! 2. initialization of atmospheric model constants and variables 
    call ini_atm(cexp)

    ! 3. initialization of coupled modules (land, sea, ice)
    call ini_coupler(istart)

    ! 4. set up the forcing fields for the first time step
    call fordate(0)

    ! 5. do the initial (2nd-order) time step, initialize the semi-impl. scheme
    call stepone
end subroutine

subroutine newdate(imode)
    !--   subroutine newdate (imode)
    !--   purpose:   initilialize and update date variables 
    !--   input :    imode = 0 for initialization, > 0 for update  

    implicit none

    integer, intent(in) :: imode
    integer, parameter :: ncal = 365
    integer :: jm, im
      
    include "com_date.h"
    include "com_tsteps.h"

    ! 365-day calendar
    integer :: ncal365(12) = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)

    if (imode <= 0) then
        ! calendar
        if (ncal == 365) then
            ndaycal(:,1) = ncal365(:)
        else
            ndaycal(:,1) = 30
        endif
 
        ndaycal(1,2) = 0
        do jm = 2, 12
            ndaycal(jm,2) = ndaycal(jm-1,1)+ndaycal(jm-1,2)
        enddo

        ! total no. of integration days
        ndaytot = ndaysl
        im = imont0
      
        do jm=1,nmonts
            ndaytot = ndaytot+ndaycal(im,1)
            im = im+1
            if (im.eq.13) im=1
        enddo

        ! initial date
        iyear  = iyear0
        imonth = imont0
        iday   = 1
    else
        ! set new date
        iday = iday+1

        if (iday > ndaycal(imonth,1)) then
            iday   = 1
            imonth = imonth+1
        endif

        if (imonth > 12) then
            imonth = 1
            iyear  = iyear+1
        endif
    endif

    ! additional variables to define forcing terms and boundary cond.
    if (iseasc >= 1) then
        imont1 = imonth
        tmonth = (iday-0.5)/float(ndaycal(imonth,1))
        tyear  = (ndaycal(imonth,2)+iday-0.5)/float(ncal)
    else
        imont1 = imont0
        tmonth = 0.5
        tyear  = (ndaycal(imont1,2)+0.5*ndaycal(imont1,2))/float(ncal)
    endif
end subroutine
