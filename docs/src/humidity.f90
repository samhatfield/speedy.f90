!> author: Sam Hatfield, Fred Kucharski, Franco Molteni
!  date: 07/05/2019
!  For converting between specific and relative humidity, and computing the
!  saturation specific humidity.
module humidity
    use params, only: ix, il

    implicit none

    private
    public spec_hum_to_rel_hum, rel_hum_to_spec_hum, get_qsat

contains
    !> Converts specific humidity to relative humidity, and also returns
    !  saturation specific humidity.
    subroutine spec_hum_to_rel_hum(ta, ps, sig, qa, rh, qsat)
        real, intent(in) :: ta(ix,il)      !! Absolute temperature
        real, intent(in) :: ps(ix,il)      !! Normalized pressure (p/1000 hPa)
        real, intent(in) :: sig            !! Sigma level
        real, intent(in) :: qa(ix,il)      !! Specific humidity
        real, intent(inout) :: rh(ix,il)   !! Relative humidity
        real, intent(inout) :: qsat(ix,il) !! Saturation specific humidity

        qsat = get_qsat(ta, ps, sig)
        rh = qa/qsat
    end subroutine

    !> Converts relative humidity to specific humidity, and also returns
    !  saturation specific humidity.
    subroutine rel_hum_to_spec_hum(ta, ps, sig, rh, qa, qsat)
        real, intent(in) :: ta(ix,il)      !! Absolute temperature
        real, intent(in) :: ps(ix,il)      !! Normalized pressure (p/1000 hPa)
        real, intent(in) :: sig            !! Sigma level
        real, intent(in) :: rh(ix,il)      !! Relative humidity
        real, intent(inout) :: qa(ix,il)   !! Specific humidity
        real, intent(inout) :: qsat(ix,il) !! Saturation specific humidity

        qsat = get_qsat(ta, ps, sig)
        qa = rh*qsat
    end subroutine

    !> Computes saturation specific humidity.
    function get_qsat(ta, ps, sig) result(qsat)
        real, intent(in) :: ta(ix,il) !! Absolute temperature
        real, intent(in) :: ps(ix,il) !! Normalized pressure (p/1000 hPa)
        real, intent(in) :: sig       !! Sigma level
        real :: qsat(ix,il)           !! Saturation specific humidity in g/kg

        real :: e0, c1, c2, t0, t1, t2

        integer :: i, j

        ! 1. Compute Qsat (g/kg) from T (degK) and normalized pres. P (= p/1000_hPa)
        ! If sig > 0, P = Ps * sigma, otherwise P = Ps(1) = const.
        e0 = 6.108d-3
        c1 = 17.269
        c2 = 21.875
        t0 = 273.16
        t1 = 35.86
        t2 = 7.66

        do i = 1, ix
            do j = 1, il
                if (ta(i,j) >= t0) then
                    qsat(i,j) = e0*exp(c1*(ta(i,j) - t0)/(ta(i,j) - t1))
                else
                    qsat(i,j) = e0*exp(c2*(ta(i,j) - t0)/(ta(i,j) - t2))
                end if
            end do
        end do

        if (sig <= 0.0) then
            qsat = 622.0*qsat/(ps(1,1) - 0.378*qsat)
        else
            qsat = 622.0*qsat/(sig*ps - 0.378*qsat)
        end if
    end function
end module
