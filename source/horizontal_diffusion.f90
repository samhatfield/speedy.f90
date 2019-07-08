!> author: Sam Hatfield, Fred Kucharski, Franco Molteni
!  date: 07/05/2019
!  For performing horizontal diffusion.
module horizontal_diffusion
    use types, only: p
    use params

    implicit none

    private
    public initialize_horizontal_diffusion, do_horizontal_diffusion
    public dmp, dmpd, dmps, dmp1, dmp1d, dmp1s, tcorv, qcorv, tcorh, qcorh

    interface do_horizontal_diffusion
        module procedure do_horizontal_diffusion_2d
        module procedure do_horizontal_diffusion_3d
    end interface

    real(p) :: dmp(mx,nx)  !! Damping coefficient for temperature and vorticity (explicit)
    real(p) :: dmpd(mx,nx) !! Damping coefficient for divergence (explicit)
    real(p) :: dmps(mx,nx) !! Damping coefficient for extra diffusion in the stratosphere (explicit)

    real(p) :: dmp1(mx,nx)  !! Damping coefficient for temperature and vorticity (implicit)
    real(p) :: dmp1d(mx,nx) !! Damping coefficient for divergence (implicit)
    real(p) :: dmp1s(mx,nx) !! Damping coefficient for extra diffusion in the stratosphere
                            !! (implicit)

    real(p) :: tcorv(kx) !! Vertical component of orographic correction for temperature
    real(p) :: qcorv(kx) !! Vertical component of orographic correction for humidity

    complex(p) :: tcorh(mx,nx) !! Horizontal component of orographic correction for temperature
    complex(p) :: qcorh(mx,nx) !! Horizontal component of orographic correction for humidity

contains
    !> Initializes the arrays used for horizontal diffusion.
    subroutine initialize_horizontal_diffusion
        use dynamical_constants, only: thd, thdd, thds, gamma, hscale, hshum
        use physical_constants, only: grav, rgas
        use geometry, only: fsg

        integer :: j, k, npowhd
        real(p) :: elap, elapn, hdifd, hdiff, hdifs, qexp, rgam, rlap, twn

        ! 1. Definition of constants
        if (mod(nsteps,2) /= 0) stop ' Invalid no. of time steps'

        ! Power of Laplacian in horizontal diffusion
        npowhd = 4

        ! Coefficients for horizontal diffusion
        ! Spectral damping coefficients
        hdiff = 1./(thd *3600.)
        hdifd = 1./(thdd*3600.)
        hdifs = 1./(thds*3600.)
        rlap  = 1./float(trunc*(trunc+1))

        do j = 1, nx
            do k = 1, mx
                twn = float(k +j - 2)
                elap = (twn*(twn+1.)*rlap)
                elapn = elap**npowhd
                dmp(k,j)  = hdiff*elapn
                dmpd(k,j) = hdifd*elapn
                dmps(k,j) = hdifs*elap
            end do
            ! dmps(1,j)=0.
        end do

        ! 5.2 Orographic correction terms for temperature and humidity
        !     (vertical component)
        rgam = rgas*gamma/(1000.*grav)
        qexp = hscale/hshum

        tcorv(1)=0.
        qcorv(1)=0.
        qcorv(2)=0.

        do k = 2, kx
            tcorv(k) = fsg(k)**rgam
            if (k.gt.2) qcorv(k) = fsg(k)**qexp
        end do
    end subroutine

    !> Adds horizontal diffusion tendency of field to spectral tendency fdt
    !  using damping coefficients dmp and dmp1.
    function do_horizontal_diffusion_2d(field, fdt_in, dmp, dmp1) result(fdt_out)
        complex(p), intent(in) :: field(mx,nx), fdt_in(mx,nx)
        complex(p) :: fdt_out(mx,nx)
        real(p), intent(in) :: dmp(mx,nx), dmp1(mx,nx)

        fdt_out = (fdt_in - dmp*field)*dmp1
    end function

    !> Adds horizontal diffusion tendency of field to spectral tendency fdt
    !  at all model levels using damping coefficients dmp and dmp1.
    function do_horizontal_diffusion_3d(field, fdt_in, dmp, dmp1) result(fdt_out)
        complex(p), intent(in) :: field(mx,nx,kx), fdt_in(mx,nx,kx)
        complex(p) :: fdt_out(mx,nx,kx)
        real(p), intent(in) :: dmp(mx,nx), dmp1(mx,nx)
        integer :: k

        do k = 1, kx
            fdt_out(:,:,k) = do_horizontal_diffusion_2d(field(:,:,k), fdt_in(:,:,k), dmp, dmp1)
        end do
    end function
end module
