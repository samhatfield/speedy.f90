!> @author
!> Sam Hatfield, AOPP, University of Oxford
!> @brief
!> A module for computing SPPT patterns to be used as multiplicative noise applied to physical tendencies
!> Stochastically Perturbed Parametrization Tendencies (SPPT) is a parametrization of model error.
!> See ECMWF Tech. Memo. #598 (Palmer et al. 2009)
module mod_sppt
    use mod_atparam
    use mod_tsteps, only: nsteps
    use mod_dyncon1, only: rearth

    implicit none

    private
    public mu, gen_sppt

    ! Array for tapering value of SPPT in the different layers of the atmosphere
    real :: mu(kx) = (/ 0, 1, 1, 1, 1, 1, 1, 0 /)

    complex :: sppt_spec(mx,nx,kx)
    logical :: first = .true.

    ! Time autocorrelation of spectral AR(1) signals
    real :: phi

    ! Decorrelation time of SPPT perturbation (in hours)
    real, parameter :: time_decorr = 6.0

    ! Correlation length scale of SPPT perturbation (in metres)
    real, parameter :: len_decorr = 500000.0

    ! Standard deviation of SPPT perturbation (in grid point space)
    real, parameter :: stddev = 0.20

    ! Total wavenumber-wise standard deviation of spectral signals
    real :: sigma(mx,nx,kx)

    contains
        !> @brief
        !> Generate grid point space SPPT pattern
        !> distribution.
        !> @return sppt_grid the generated grid point pattern
        function gen_sppt() result(sppt_grid)
            integer :: m, n, k
            real :: sppt_grid(ix,il,kx)
            complex :: eta(mx,nx,kx)
            real :: f0, randreal, randimag, twn

            ! Generate spatial amplitude pattern and time correlation
            if (first) then
                f0 = sum((/ ((2*n+1)*exp(-0.5*(len_decorr/rearth)**2*n*(n+1)),n=1,ntrun) /))
                f0 = sqrt((stddev**2*(1-phi**2))/(2*f0))

                do m = 1, mx
                    do n = 1, nx
                        ! Compute total wave number (I don't understand how this works)
                        twn=isc*(m-1)+n-1
                        
                        sigma(m,n,:) = f0 * exp(-0.25*(len_decorr/rearth)**2*twn*(twn+1))
                    end do
                end do

                phi = exp(-(24/real(nsteps))/time_decorr)
            end if

            ! Generate Gaussian noise
            do m = 1,mx
                do n = 1,nx
                    do k = 1,kx
                        randreal = randn(0.0, 1.0)
                        randimag = randn(0.0, 1.0)
                        eta(m,n,k) = cmplx(&
                            & min(10.0, abs(randreal)) * randreal/abs(randreal),&
                            & min(10.0, abs(randimag)) * randimag/abs(randimag))
                    end do
                end do
            end do

            ! If first timestep
            if (first) then
                ! First AR(1) step
                sppt_spec = (1 - phi**2)**(-0.5) * sigma * eta

                first = .false.
            else
                ! Subsequent AR(1) steps
                sppt_spec = phi*sppt_spec + sigma*eta
            end if

            ! Convert to grid point space
             do k=1,kx
                 call grid(sppt_spec(:,:,k),sppt_grid(:,:,k),1)
             end do

             ! Clip to +/- 3 standard deviations
             sppt_grid = min(3.0*stddev, abs(sppt_grid)) * sppt_grid/abs(sppt_grid)
        end function

        !> @brief
        !> Generates a random number drawn for the specified normal
        !> distribution.
        !> @param mean the mean of the distribution to draw from
        !> @param stdev the standard deviation of the distribution to draw from
        !> @return randn the generated random number
        function randn(mean, stdev)
            real, intent(in) :: mean, stdev
            real :: u, v, randn
            real :: rand(2)

            call random_number(rand)

            ! Box-Muller method
            u = (-2.0 * log(rand(1))) ** 0.5
            v =   2.0 * 6.28318530718 * rand(2)
            randn = mean + stdev * u * sin(v)
        end function
end module
