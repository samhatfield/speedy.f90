module geometry
    use mod_atparam

    implicit none

    private
    public initialize_geometry
    public hsg, dhs, fsg, dhsr, fsgr, radang, gsin, gcos, coriol, sia, coa, cosg, cosgr, cosgr2

    ! Vertical level parameters
    real :: hsg(kxp), dhs(kx), fsg(kx), dhsr(kx), fsgr(kx)

    ! Functions of latitude and longitude
    real, dimension(il) :: radang, gsin, gcos, coriol
    real, dimension(iy) :: sia, coa
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
            sia(j) = cos(3.141592654*(j - 0.25)/(il + 0.5))
            coa(j) = sqrt(1.0 - sia(j)**2.0)
            jj = il + 1 - j
            radang(j)  = -asin(sia(j))
            radang(jj) =  asin(sia(j))
            gsin(j)    = -sia(j)
            gsin(jj)   =  sia(j)
        end do

        ! Expand cosine and its reciprocal to cover both hemispheres
        do j=1,iy
            jj=il+1-j
            cosg(j)=coa(j)
            cosg(jj)=coa(j)
            cosgr(j)=1./coa(j)
            cosgr(jj)=1./coa(j)
            cosgr2(j)=1./(coa(j)*coa(j))
            cosgr2(jj)=1./(coa(j)*coa(j))
        end do

        do j = 1, il
            gcos(j) = cosg(j)
            coriol(j) = 2.*omega*gsin(j)
        end do
    end subroutine
end module