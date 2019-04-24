module geometry
    use mod_atparam

    implicit none

    private
    public initialize_geometry
    public hsg, dhs, fsg, dhsr, fsgr, radang, gsin, gcos, coriol, sia, coa, sia_half, coa_half, &
        cosg, cosgr, cosgr2

    ! Vertical level parameters
    real :: hsg(kxp), dhs(kx), fsg(kx), dhsr(kx), fsgr(kx)

    ! Functions of latitude and longitude
    real, dimension(il) :: radang, gsin, gcos, coriol, sia, coa
    real, dimension(iy) :: sia_half, coa_half
    real, dimension(il) :: cosg, cosgr, cosgr2

contains
    subroutine initialize_geometry
        use mod_dyncon1, only: akap, omega

        integer i, j, jj, k

        ! Definition of model levels
        ! Half (vertical velocity) levels
        if (kx == 5) then
            hsg(:6) = (/ 0.000, 0.150, 0.350, 0.650, 0.900, 1.000 /)
        else if (kx == 7) then
            hsg(:8) = (/ 0.020, 0.140, 0.260, 0.420, 0.600, 0.770, 0.900, 1.000 /)
        else if (kx == 8) then
            hsg(:9) = (/ 0.000, 0.050, 0.140, 0.260, 0.420, 0.600, 0.770, 0.900, 1.000 /)
        end if

        ! Layer thicknesses and full (u,v,T) levels
        do k = 1, kx
            dhs(k) = hsg(k+1)-hsg(k)
            fsg(k) = 0.5*(hsg(k+1)+hsg(k))
        end do

        ! Additional functions of sigma
        do k = 1, kx
            dhsr(k) = 0.5/dhs(k)
            fsgr(k) = akap/(2.*fsg(k))
        end do

        ! Horizontal functions

        ! Latitudes and functions of latitude
        ! NB: J=1 is Southernmost point!
        do j = 1, iy
            jj = il + 1 - j
            sia_half(j) = cos(3.141592654*(j - 0.25)/(il + 0.5))
            coa_half(j) = sqrt(1.0 - sia_half(j)**2.0)
            sia(j)  = -sia_half(j)
            sia(jj) =  sia_half(j)
            coa(j)  = coa_half(j)
            coa(jj) = coa_half(j)
            radang(j)  = -asin(sia_half(j))
            radang(jj) =  asin(sia_half(j))
            gsin(j)    = -sia_half(j)
            gsin(jj)   =  sia_half(j)
        end do

        ! Expand cosine and its reciprocal to cover both hemispheres
        do j=1,iy
            jj=il+1-j
            cosg(j)=coa_half(j)
            cosg(jj)=coa_half(j)
            cosgr(j)=1./coa_half(j)
            cosgr(jj)=1./coa_half(j)
            cosgr2(j)=1./(coa_half(j)*coa_half(j))
            cosgr2(jj)=1./(coa_half(j)*coa_half(j))
        end do

        do j = 1, il
            gcos(j) = cosg(j)
            coriol(j) = 2.*omega*gsin(j)
        end do
    end subroutine
end module