!> author: Sam Hatfield, Fred Kucharski, Franco Molteni
!  date: 01/05/2019
!  For keeping track of the model's date and time.
module date
    implicit none

    private
    public model_datetime, start_datetime, end_datetime
    public imont1, tmonth, tyear, ndaycal
    public isst0
    public datetime_equal, newdate

    !> For storing dates and times.
    type datetime
        integer :: year
        integer :: month
        integer :: day
        integer :: hour
        integer :: minute
    end type

    ! Date and time variables
    type(datetime)     :: model_datetime !! The model's current datetime (continuously updated)
    type(datetime)     :: start_datetime !! The start datetime
    type(datetime)     :: end_datetime   !! The end datetime
    integer            :: imont1         !! The month used for computing seasonal forcing fields
    real               :: tmonth         !! The fraction of the current month elapsed
    real               :: tyear          !! The fraction of the current year elapsed
    integer            :: isst0          !! Initial month of SST anomalies
    integer            :: ndaycal(12,2)  !! The model calendar
    integer, parameter :: ncal = 365     !! The number of days in a year

    !> The number of days in each month
    integer :: ncal365(12) = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)

    contains
        !> Checks whether two datetimes are equal.
        logical function datetime_equal(datetime1, datetime2)
            type(datetime), intent(in) :: datetime1, datetime2

            if (datetime1%year == datetime2%year .and. &
                datetime1%month == datetime2%month .and. &
                datetime1%day == datetime2%day .and. &
                datetime1%hour == datetime2%hour .and. &
                datetime1%minute == datetime2%minute) then
                datetime_equal = .true.
            else
                datetime_equal = .false.
            end if
        end function

        !> Updates the current datetime and related date variables.
        subroutine newdate(imode)
            use params, only: iseasc, nsteps

            integer, intent(in) :: imode !! Mode -> 0 for initialization, > 0 for update

            integer :: jm

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
            else
                ! Increment minute counter
                model_datetime%minute = model_datetime%minute + int(24*60/nsteps)

                ! Increment hour counter if necessary
                if (model_datetime%minute >= 60) then
                    model_datetime%minute = mod(model_datetime%minute,60)
                    model_datetime%hour = model_datetime%hour + 1
                end if

                ! Increment day counter if necessary
                if (model_datetime%hour >= 24) then
                    model_datetime%hour = mod(model_datetime%hour,24)
                    model_datetime%day = model_datetime%day+1
                end if

                ! Increment month counter if necessary
                ! Leap year and February?
                if (mod(model_datetime%year,4) == 0 .and. model_datetime%month == 2) then
                    if (model_datetime%day > 29) then
                        model_datetime%day   = 1
                        model_datetime%month = model_datetime%month + 1
                    end if
                else
                    if (model_datetime%day > ndaycal(model_datetime%month,1)) then
                        model_datetime%day   = 1
                        model_datetime%month = model_datetime%month+1
                    end if
                end if

                ! Increment year counter if necessary
                if (model_datetime%month > 12) then
                    model_datetime%month = 1
                    model_datetime%year  = model_datetime%year+1
                end if
            end if

            ! additional variables to define forcing terms and boundary cond.
            if (iseasc >= 1) then
                imont1 = model_datetime%month
                tmonth = (model_datetime%day-0.5)/float(ndaycal(model_datetime%month,1))
                tyear  = (ndaycal(model_datetime%month,2)+model_datetime%day-0.5)/float(ncal)
            else
                imont1 = start_datetime%month
                tmonth = 0.5
                tyear  = (ndaycal(imont1,2)&
                    & +0.5*ndaycal(imont1,2))/float(ncal)
            end if
        end subroutine
end module
