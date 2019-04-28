module horizontal_diffusion
    use params

    implicit none

    private
    public initialize_horizontal_diffusion, do_horizontal_diffusion
    public dmp, dmpd, dmps, dmp1, dmp1d, dmp1s, tcorv, qcorv, tcorh, qcorh

    interface do_horizontal_diffusion
        module procedure do_horizontal_diffusion_2d
        module procedure do_horizontal_diffusion_3d
    end interface

    ! Damping coef. for horizontal diffusion (explicit) (initial. in indyns)
    real, dimension(mx,nx) :: dmp, dmpd, dmps

    ! Damping coef. for horizontal diffusion (implicit) (initial. in indyns)
    real, dimension(mx,nx) :: dmp1, dmp1d, dmp1s

    ! Vertical comp. of orographic correction (initial. in INDYNS)
    real, dimension(kx) :: tcorv, qcorv

    ! Horizontal component of orographic correction (updated in FORDATE)
    complex, dimension(mx,nx) :: tcorh, qcorh

contains
    subroutine initialize_horizontal_diffusion
        use mod_dyncon0
        use physical_constants, only: grav, rgas
        use geometry, only: fsg

        integer :: j, k, npowhd
        real :: elap, elapn, hdifd, hdiff, hdifs, qexp, rgam, rlap, twn

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

    ! Add horizontal diffusion tendency of FIELD to spectral tendency FDT at NLEV
    ! levels using damping coefficients DMP and DMP1
    function do_horizontal_diffusion_2d(field, fdt_in, dmp, dmp1) result(fdt_out)
        complex, intent(in) :: field(mx,nx), fdt_in(mx,nx)
        complex :: fdt_out(mx,nx)
        real, intent(in) :: dmp(mx,nx), dmp1(mx,nx)

        fdt_out = (fdt_in - dmp*field)*dmp1
    end function

    ! Add horizontal diffusion tendency of FIELD to spectral tendency FDT at NLEV
    ! levels using damping coefficients DMP and DMP1
    function do_horizontal_diffusion_3d(field, fdt_in, dmp, dmp1) result(fdt_out)
        complex, intent(in) :: field(mx,nx,kx), fdt_in(mx,nx,kx)
        complex :: fdt_out(mx,nx,kx)
        real, intent(in) :: dmp(mx,nx), dmp1(mx,nx)
        integer :: k

        do k = 1, kx
            fdt_out(:,:,k) = do_horizontal_diffusion_2d(field(:,:,k), fdt_in(:,:,k), dmp, dmp1)
        end do
    end function

end module
