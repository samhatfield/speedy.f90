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
    use mod_spectral, only: el2

    implicit none

    private
    public mu, gen_sppt

    ! Array for tapering value of SPPT in the different layers of the atmosphere
    ! A value of 1 means the tendency is not tapered at that level
    real :: mu(kx) = (/ 1, 1, 1, 1, 1, 1, 1, 1 /)

    complex :: sppt_spec(mx,nx,kx)
    logical :: first = .true.

    ! Decorrelation time of SPPT perturbation (in hours)
    real, parameter :: time_decorr = 6.0

    ! Time autocorrelation of spectral AR(1) signals
    real :: phi = exp(-(24/real(nsteps))/time_decorr)

    ! Correlation length scale of SPPT perturbation (in metres)
    real, parameter :: len_decorr = 500000.0

    ! Standard deviation of SPPT perturbation (in grid point space)
    real, parameter :: stddev = 0.33

    ! Total wavenumber-wise standard deviation of spectral signals
    real :: sigma(mx,nx,kx)

    contains
        !> @brief
        !> Generate grid point space SPPT pattern
        !> distribution.
        !> @return sppt_grid the generated grid point pattern
        function gen_sppt() result(sppt_grid_out)
            integer :: m, n, k
            real :: sppt_grid(ix,il,kx), sppt_grid_out(ix*il,kx)
            complex :: eta(mx,nx,kx)
            real :: f0, randreal, randimag, twn

            ! Seed RNG if first use of SPPT
            if (first) call time_seed()

            ! Generate Gaussian noise
            do m = 1,mx
                do n = 1,nx
                    do k = 1,kx
                        randreal = randn(0.0, 1.0)
                        randimag = randn(0.0, 1.0)

                        ! Clip noise to +- 10 standard deviations
                        eta(m,n,k) = cmplx(&
                            & min(10.0, abs(randreal)) * sign(1.0,randreal),&
                            & min(10.0, abs(randimag)) * sign(1.0,randimag))
                    end do
                end do
            end do

            ! If first timestep
            if (first) then
                ! Generate spatial amplitude pattern and time correlation
                f0 = sum((/ ((2*n+1)*exp(-0.5*(len_decorr/rearth)**2*n*(n+1)),n=1,ntrun) /))
                f0 = sqrt((stddev**2*(1-phi**2))/(2*f0))

                do k = 1,kx
                    sigma(:,:,k) = f0 * exp(-0.25*len_decorr**2 * el2)
                end do

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
                 sppt_grid_out(:,k) = reshape(sppt_grid(:,:,k), (/ix*il/))
             end do

             ! Clip to +/- 1.0
             sppt_grid_out = min(1.0, abs(sppt_grid_out)) * sign(1.0,sppt_grid_out)
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

        !> @brief
        !> Seeds RNG from system clock.
        subroutine time_seed()
            integer :: i, n, clock
            integer, allocatable :: seed(:)
          
            call random_seed(size = n)
            allocate(seed(n))
          
            call system_clock(count=clock)
          
            seed = clock + 37 * (/ (i - 1, i = 1, n) /)
            call random_seed(put = seed)
          
            deallocate(seed)
        end subroutine
end module
