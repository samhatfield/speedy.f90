!> author: Sam Hatfield, Fred Kucharski, Franco Molteni
!  date: 29/04/2019
!  For reading and storing boundary conditions.
!
!  This module handles the reading and storing of the land-sea mask, the
!  surface geopotential (i.e. the orography), the filtered surface geopotential
!  (i.e. the smoothed orography) and the bare-land annual-mean albedo.
module boundaries
    use params

    implicit none

    private
    public initialize_boundaries, fillsf, forchk
    public fmask, phi0, phis0, alb0

    real :: fmask(ix,il) !! Original (fractional) land-sea mask

    ! Time invariant surface fields
    real :: phi0(ix,il)  !! Unfiltered surface geopotential
    real :: phis0(ix,il) !! Spectrally-filtered surface geopotential
    real :: alb0(ix,il)  !! Bare-land annual-mean albedo

contains
    !> Initialize boundary conditions (land-sea mask, surface geopotential
    !  and surface albedo).
    subroutine initialize_boundaries
        use physical_constants, only: grav
        use input_output, only: load_boundary_file

        ! Read surface geopotential (i.e. orography)
        phi0 = grav*load_boundary_file("surface.nc", "orog")

        ! Also store spectrally truncated surface geopotential
        call spectral_truncation(phi0, phis0)

        ! Read land-sea mask
        fmask = load_boundary_file("surface.nc", "lsm")

        ! Annual-mean surface albedo
        alb0 = load_boundary_file("surface.nc", "alb")
    end subroutine

    !> Check consistency of surface fields with land-sea mask and set undefined
    !  values to a constant (to avoid over/underflow).
    subroutine forchk(fmask, nf, fmin, fmax, fset, field)
        real, intent(in)    :: fmask(ix,il)    !! The fractional land-sea mask
        integer, intent(in) :: nf              !! The number of input 2D fields
        real, intent(in)    :: fmin            !! The minimum allowable value
        real, intent(in)    :: fmax            !! The maximum allowable value
        real, intent(in)    :: fset            !! Replacement for undefined values
        real, intent(inout) :: field(ix,il,nf) !! The output field

        integer :: i, j, jf, nfault

        do jf = 1, nf
            nfault = 0

            do i = 1, ix
                do j = 1, il
                    if (fmask(i,j) > 0.0) then
                        if (field(i,j,jf) < fmin .or. field(i,j,jf) > fmax) then
                            nfault = nfault + 1
                        end if
                    else
                        field(i,j,jf) = fset
                    end if
                end do
            end do
        end do
    end subroutine

    !> Compute a spectrally-filtered grid-point field.
    subroutine spectral_truncation(fg1, fg2)
        use spectral, only: grid_to_spec, spec_to_grid

        real, intent(inout) :: fg1(ix,il) !! Original grid-point field
        real, intent(inout) :: fg2(ix,il) !! Filtered grid-point field

        complex :: fsp(mx,nx)
        integer :: n, m, total_wavenumber

        fsp = grid_to_spec(fg1)

        do n = 1, nx
            do m = 1, mx
                total_wavenumber = m + n - 2
                if (total_wavenumber > trunc) fsp(m,n) = (0.0, 0.0)
            end do
        end do

        fg2 = spec_to_grid(fsp, 1)
    end subroutine

    !> Replace missing values in surface fields.
    !  @note It is assumed that non-missing values exist near the Equator.
    subroutine fillsf(sf, fmis)
        real, intent(inout) :: sf(ix,il) !! Field to replace missing values in
        real, intent(in)    :: fmis      !! Replacement for missing values

        real :: sf2(0:ix+1)

        integer :: hemisphere, j, j1, j2, j3, i, nmis
        real :: fmean

        do hemisphere = 1, 2
           if (hemisphere == 1) then
                j1 = il/2
                j2 = 1
                j3 = -1
            else
                j1 = j1+1
                j2 = il
                j3 = 1
            end if

            do j = j1, j2, j3
                sf2(1:ix) = sf(:,j)

                nmis = 0
                do i = 1, ix
                    if (sf(i,j) < fmis) then
                        nmis = nmis + 1
                        sf2(i) = 0.0
                    end if
                end do

                if (nmis < ix) fmean = sum(sf2(1:ix))/float(ix - nmis)

                do i = 1, ix
                    if (sf(i,j) < fmis) sf2(i) = fmean
                end do

                sf2(0)      = sf2(ix)
                sf2(ix+1) = sf2(1)
                do i = 1, ix
                    if (sf(i,j) < fmis) sf(i,j) = 0.5*(sf2(i-1) + sf2(i+1))
                end do
            end do
        end do
    end subroutine
end module
