module mod_date
    implicit none

    private
    public iyear, imonth, iday, imont1, tmonth, tyear, ndaycal, ndaytot
    public newdate

    ! Date and time variables (updated in NEWDATE)
    integer :: iyear, imonth, iday, imont1
    real :: tmonth, tyear

    ! Calendar set-up (initialized in NEWDATE)
    integer :: ndaycal(12,2), ndaytot

    contains
        subroutine newdate(imode)
            !--   subroutine newdate (imode)
            !--   purpose:   initilialize and update date variables 
            !--   input :    imode = 0 for initialization, > 0 for update  
        
            use mod_tsteps
        
            implicit none
        
            integer, intent(in) :: imode
            integer, parameter :: ncal = 365
            integer :: jm, im
              
            ! 365-day calendar
            integer :: ncal365(12) = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31,&
                & 30, 31 /)
        
            if (imode <= 0) then
                ! calendar
                if (ncal == 365) then
                    ndaycal(:,1) = ncal365(:)
                else
                    ndaycal(:,1) = 30
                end if
         
                ndaycal(1,2) = 0
                do jm = 2, 12
                    ndaycal(jm,2) = ndaycal(jm-1,1)+ndaycal(jm-1,2)
                end do
        
                ! total no. of integration days
                ndaytot = ndaysl
                im = imont0
              
                do jm=1,nmonts
                    ndaytot = ndaytot+ndaycal(im,1)
                    im = im+1
                    if (im.eq.13) im=1
                end do
        
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
                end if
        
                if (imonth > 12) then
                    imonth = 1
                    iyear  = iyear+1
                end if
            end if
        
            ! additional variables to define forcing terms and boundary cond.
            if (iseasc >= 1) then
                imont1 = imonth
                tmonth = (iday-0.5)/float(ndaycal(imonth,1))
                tyear  = (ndaycal(imonth,2)+iday-0.5)/float(ncal)
            else
                imont1 = imont0
                tmonth = 0.5
                tyear  = (ndaycal(imont1,2)&
                    & +0.5*ndaycal(imont1,2))/float(ncal)
            end if
        end subroutine
end module
