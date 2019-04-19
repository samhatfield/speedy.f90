module humidity
    use mod_atparam, only: ix, il

    implicit none

    private
    public spec_hum_to_rel_hum, rel_hum_to_spec_hum, get_qsat

contains
    subroutine spec_hum_to_rel_hum(ta, ps, sig, qa, rh, qsat)
        real, intent(in) :: ta(ix,il), ps(ix,il), sig, qa(ix,il)
        real, intent(inout) :: rh(ix,il), qsat(ix,il)

        qsat = get_qsat(ta, ps, sig)
        rh = qa/qsat
    end subroutine

    subroutine rel_hum_to_spec_hum(ta, ps, sig, rh, qa, qsat)
        real, intent(in) :: ta(ix,il), ps(ix,il), sig, rh(ix,il)
        real, intent(inout) :: qa(ix,il), qsat(ix,il)

        qsat = get_qsat(ta, ps, sig)
        qa = rh*qsat
    end subroutine

    ! Compute saturation specific humidity
    ! Input:   ta     : abs. temperature
    !          ps     : normalized pressure   (=  p/1000_hPa)
    !                 : normalized sfc. pres. (= ps/1000_hPa)
    !          sig    : sigma level
    !          qa     : specific humidity in g/kg
    !          rh     : relative humidity
    ! Output:  qsat   : saturation spec. hum. in g/kg
    function get_qsat(ta, ps, sig) result(qsat)
        real, intent(in) :: ta(ix,il), ps(ix,il), sig
        real :: qsat(ix,il), e0, c1, c2, t0, t1, t2

        integer :: i,j

        ! 1. Compute Qsat (g/kg) from T (degK) and normalized pres. P (= p/1000_hPa)
        ! If sig > 0, P = Ps * sigma, otherwise P = Ps(1) = const.
        e0 = 6.108e-3
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
