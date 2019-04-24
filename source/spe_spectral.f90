!*********************************************************************
subroutine grad(psi,psdx,psdy)
    use mod_atparam
    use spectral, only: gradx, gradyp, gradym

    implicit none

    complex, dimension(mx,nx), intent(inout) :: psi
    complex, dimension(mx,nx), intent(inout) :: psdx, psdy

    integer :: k, n, m

    do n = 1, nx
        psdx(:,n) = gradx*psi(:,n)*(0.0, 1.0)
    end do

    do m=1,mx
        psdy(m,1)  =  gradyp(m,1)*psi(m,2)
        psdy(m,nx) = -gradym(m,nx)*psi(m,ntrun1)
    end do

    do n=2,ntrun1
        do m=1,mx
            psdy(m,n) = -gradym(m,n)*psi(m,n-1) + gradyp(m,n)*psi(m,n+1)
        end do
    end do
end
!******************************************************************
subroutine vds(ucosm,vcosm,vorm,divm)
    use mod_atparam
    use spectral, only: gradx, vddyp, vddym

    implicit none

    complex, dimension(mx,nx) :: ucosm, vcosm
    complex, dimension(mx,nx), intent(inout) :: vorm, divm
    complex, dimension(mx,nx) :: zc, zp

    integer :: n, m, k

    do n=1,nx
        zp(:,n) = gradx*ucosm(:,n)*(0.0, 1.0)
        zc(:,n) = gradx*vcosm(:,n)*(0.0, 1.0)
    end do

    do m=1,mx
        vorm(m,1)  = zc(m,1) - vddyp(m,1)*ucosm(m,2)
        vorm(m,nx) = vddym(m,nx)*ucosm(m,ntrun1)
        divm(m,1)  = zp(m,1) + vddyp(m,1)*vcosm(m,2)
        divm(m,nx) = -vddym(m,nx)*vcosm(m,ntrun1)
    end do

    do n=2,ntrun1
        do m=1,mx
            vorm(m,n) =  vddym(m,n)*ucosm(m,n-1) - vddyp(m,n)*ucosm(m,n+1) + zc(m,n)
            divm(m,n) = -vddym(m,n)*vcosm(m,n-1) + vddyp(m,n)*vcosm(m,n+1) + zp(m,n)
        end do
    end do
end
!******************************************************************
subroutine uvspec(vorm,divm,ucosm,vcosm)
    use mod_atparam
    use spectral, only: uvdx, uvdyp, uvdym

    complex, dimension(mx,nx), intent(in) :: vorm,divm
    complex, dimension(mx,nx), intent(inout) :: ucosm,vcosm
    complex, dimension(mx,nx) :: zc,zp

    integer :: k, n, m

    zp = uvdx*vorm*(0.0, 1.0)
    zc = uvdx*divm*(0.0, 1.0)

    do m=1,mx
        ucosm(m,1)  =  zc(m,1) - uvdyp(m,1)*vorm(m,2)
        ucosm(m,nx) =  uvdym(m,nx)*vorm(m,ntrun1)
        vcosm(m,1)  =  zp(m,1) + uvdyp(m,1)*divm(m,2)
        vcosm(m,nx) = -uvdym(m,nx)*divm(m,ntrun1)
    end do

    do n=2,ntrun1
        do m=1,mx
          vcosm(m,n) = -uvdym(m,n)*divm(m,n-1) + uvdyp(m,n)*divm(m,n+1) + zp(m,n)
          ucosm(m,n) =  uvdym(m,n)*vorm(m,n-1) - uvdyp(m,n)*vorm(m,n+1) + zc(m,n)
        end do
    end do
end
!*********************************************************************
subroutine vdspec(ug,vg,vorm,divm,kcos)
    use mod_atparam
    use geometry, only: cosgr, cosgr2
    use spectral, only: grid_to_spec

    implicit none

    real, intent(in) :: ug(ix,il), vg(ix,il)
    complex, intent(out) :: vorm(mx,nx), divm(mx,nx)
    integer, intent(in) :: kcos
    integer :: i, j
    real :: ug1(ix,il), vg1(ix,il)
    complex :: specu(mx,nx), specv(mx,nx)

    if (kcos.eq.2) then
        do j=1,il
            do i=1,ix
                ug1(i,j)=ug(i,j)*cosgr(j)
                vg1(i,j)=vg(i,j)*cosgr(j)
            end do
        end do
    else
        do j=1,il
            do i=1,ix
                ug1(i,j)=ug(i,j)*cosgr2(j)
                vg1(i,j)=vg(i,j)*cosgr2(j)
            end do
        end do
    end if

    specu = grid_to_spec(ug1)
    specv = grid_to_spec(vg1)
    call vds(specu,specv,vorm,divm)
end
!******************************************************************
subroutine trunct(vor)
    use mod_atparam
    use spectral, only: trfilt

    implicit none

    complex, intent(inout) :: vor(mx,nx)

    vor = vor * trfilt
end
