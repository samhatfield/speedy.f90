subroutine implic(divdt,tdt,psdt)
    !
    !  subroutine implic (divdt,tdt,psdt)
    !
    !  Purpose : Correct tendencies for implicit gravity wave model
    !  Input/output : divdt = divergence tendency
    !                 tdt   = temperature tendency
    !                 psdt  = tendency of log(surf.pressure)
    !

    use mod_atparam
    use mod_dyncon1, only: dhs
    use mod_dyncon2, only: tref1, xc, xd, xj, dhsx, elz

    implicit none

    integer, parameter :: mxnxkx = mx*nx*kx

    complex, intent(inout) :: divdt(mx,nx,kx), tdt(mx,nx,kx), psdt(mx,nx)
    complex ::  ye(mx,nx,kx), yf(mx,nx,kx), zero
    integer :: k1, k, m, n, ll, mm

    zero = (0.,0.)

    ye(:,:,:) = zero

    do k1=1,kx
        do k=1,kx
            ye(:,:,k) = ye(:,:,k) + xd(k,k1) * tdt(:,:,k1)
        end do
    end do

    do k=1,kx
        ye(:,:,k) = ye(:,:,k) + tref1(k) * psdt
    end do

    do k=1,kx
        do m=1,mx
            do n=1,nx
                yf(m,n,k)=divdt(m,n,k)+elz(m,n)*ye(m,n,k)
            end do
        end do
    end do

    divdt(:,:,:) = zero

    do n=1,nx
        do m=1,mx
            mm=isc*(m-1)+1
            ll=mm+n-2
            if(ll.ne.0) then
                do k1=1,kx
                    divdt(m,n,:) = divdt(m,n,:) + xj(:,k1,ll) * yf(m,n,k1)
                end do
            endif
        end do
    end do

    do k=1,kx
        psdt = psdt - divdt(:,:,k) * dhsx(k)
    end do

    do k=1,kx
        do k1=1,kx
            tdt(:,:,k) = tdt(:,:,k) + xc(k,k1) * divdt(:,:,k1)
        end do
    end do
end
