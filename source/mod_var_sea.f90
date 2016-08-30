module mod_var_sea
    use mod_atparam

    implicit none

    private
    public sstcl_ob, sicecl_ob, ticecl_ob, sstan_ob, sstcl_om, sst_am, sstan_am
    public sice_am, tice_am, sst_om, sice_om, tice_om, ssti_om, wsst_ob

    ! Daily observed climatological fields over sea
    ! Observed clim. SST
    real :: sstcl_ob(ix*il)

    ! Clim. sea ice fraction
    real :: sicecl_ob(ix*il)

    ! Clim. sea ice temperature
    real :: ticecl_ob(ix*il)

    ! Daily observed SST anomaly
    ! Observed SST anomaly
    real :: sstan_ob(ix*il)

    ! Daily climatological fields from ocean model
    ! Ocean model clim. SST
    real :: sstcl_om(ix*il)

    ! Sea sfc. fields used by atmospheric model
    ! SST (full-field)
    real :: sst_am(ix*il)

    ! SST anomaly
    real :: sstan_am(ix*il)

    ! Sea ice fraction
    real :: sice_am(ix*il)

    ! Sea ice temperature
    real :: tice_am(ix*il)

    ! Sea sfc. fields from ocean/sea-ice model
    ! Ocean model SST
    real :: sst_om(ix*il)

    ! Model sea ice fraction
    real :: sice_om(ix*il)

    ! Model sea ice temperature
    real :: tice_om(ix*il)

    ! Model SST + sea ice temp.
    real :: ssti_om(ix*il)

    ! Weight for obs. SST anomaly in coupled runs
    ! Weight mask for obs. SST
    real :: wsst_ob(ix*il)
end module
