subroutine shtorh(imode,ngp,ta,ps,sig,qa,rh,qsat)
    ! subroutine shtorh (imode,ngp,ta,ps,sig,qa,rh,qsat)
    !
    ! Purpose: compute saturation specific humidity and 
    !          relative hum. from specific hum. (or viceversa)
    ! Input:   imode  : mode of operation
    !          ngp    : no. of grid-points
    !          ta     : abs. temperature
    !          ps     : normalized pressure   (=  p/1000_hPa) [if sig < 0]
    !                 : normalized sfc. pres. (= ps/1000_hPa) [if sig > 0]
    !          sig    : sigma level
    !          qa     : specific humidity in g/kg [if imode > 0]
    !          rh     : relative humidity         [if imode < 0]
    !          qsat   : saturation spec. hum. in g/kg 
    ! Output:  rh     : relative humidity         [if imode > 0] 
    !          qa     : specific humidity in g/kg [if imode < 0]
    !      

    implicit none

    integer, intent(in) :: imode, ngp
    real, intent(in) :: ta(ngp), ps(*), sig
    real :: qsat(ngp), e0, c1, c2, t0, t1, t2
    real, intent(inout) :: qa(ngp), rh(ngp)

    integer :: j

    ! 1. Compute Qsat (g/kg) from T (degK) and normalized pres. P (= p/1000_hPa)
    ! If sig > 0, P = Ps * sigma, otherwise P = Ps(1) = const. 
    e0 = 6.108e-3
    c1 = 17.269
    c2 = 21.875
    t0 = 273.16
    t1 = 35.86
    t2 = 7.66
      
    do j=1,ngp
        if (ta(j).ge.t0) then
          qsat(j)=e0*exp(c1*(ta(j)-t0)/(ta(j)-t1))
        else
          qsat(j)=e0*exp(c2*(ta(j)-t0)/(ta(j)-t2))
        end if
    end do

    if (sig.le.0.0) then
        do j=1,ngp
            qsat(j)=622.*qsat(j)/(ps(1)-0.378*qsat(j))
        end do
    else
        do j=1,ngp
            qsat(j)=622.*qsat(j)/(sig*ps(j)-0.378*qsat(j))
        end do
    end if

    ! 2. Compute rel.hum. RH=Q/Qsat (imode>0), or Q=RH*Qsat (imode<0)
    if (imode.gt.0) then
        do j=1,ngp
            rh(j)=qa(j)/qsat(j)
        end do
    else if (imode.lt.0) then
        do j=1,ngp
            qa(j)=rh(j)*qsat(j)
        end do
    end if
end

subroutine zmeddy(nlon,nlat,ff,zm,eddy)
    ! Decompose a field into zonal-mean and eddy component

    implicit none

    integer, intent(in) :: nlon, nlat
    real, intent(in) :: ff(nlon,nlat)
    real, intent(inout) :: zm(nlat), eddy(nlon,nlat)
    integer :: i, j
    real :: rnlon

    rnlon=1./nlon

    do j=1,nlat
        zm(j)=0.
        do i=1,nlon
            zm(j)=zm(j)+ff(i,j)
        end do
        zm(j)=zm(j)*rnlon

        do i=1,nlon
            eddy(i,j)=ff(i,j)-zm(j)
        end do
    end do
end
