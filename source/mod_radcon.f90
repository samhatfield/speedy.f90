module mod_radcon
    use mod_atparam

    implicit none

    private
    public solc, albsea, albice, albsn, rhcl1, rhcl2, qacl, wpcl, pmaxcl,&
        & clsmax, clsminl, gse_s0, gse_s1, albcl, albcls, epssw,  epslw, emisfc,&
        & absdry, absaer, abswv1, abswv2, abscl1, abscl2, ablwin, ablco2,&
        & ablco2_ref, ablwv1, ablwv2, ablcl1, ablcl2
    public fband
    public fsol, ozone, ozupp, zenit, stratz
    public alb_l, alb_s, albsfc, snowc
    public tau2, st4a, stratc, flux
    public qcloud, irhtop

    ! Radiation and cloud constants
    
    ! solc   = Solar constant (area averaged) in W/m^2
    ! albsea = Albedo over sea 
    ! albice = Albedo over sea ice (for ice fraction = 1)
    ! albsn  = Albedo over snow (for snow cover = 1)
    
    ! rhcl1  = relative hum. threshold corr. to cloud cover = 0
    ! rhcl2  = relative hum. corr. to cloud cover = 1
    ! qacl   = specific hum. threshold for cloud cover
    ! wpcl   = cloud c. weight for the sq. root of precip. (for p = 1 mm/day)
    ! pmaxcl = max. value of precip. (mm/day) contributing to cloud cover 
    
    ! clsmax = maximum stratiform cloud cover
    ! clsminl= minimum stratiform cloud cover over land (for RH = 1)
    ! gse_s0 = gradient of dry static energy corresp. to strat.c.c. = 0
    ! gse_s1 = gradient of dry static energy corresp. to strat.c.c. = 1
    
    ! albcl  = cloud albedo (for cloud cover = 1)
    ! albcls = stratiform cloud albedo (for st. cloud cover = 1)
    ! epssw  = fraction of incoming solar radiation absorbed by ozone
    ! epslw  = fraction of blackbody spectrum absorbed/emitted by PBL only
    ! emisfc = longwave surface emissivity
    
    !          shortwave absorptivities (for dp = 10^5 Pa) :
    ! absdry = abs. of dry air      (visible band)
    ! absaer = abs. of aerosols     (visible band)
    ! abswv1 = abs. of water vapour (visible band, for dq = 1 g/kg)
    ! abswv2 = abs. of water vapour (near IR band, for dq = 1 g/kg)
    ! abscl2 = abs. of clouds       (visible band, for dq_base = 1 g/kg)
    ! abscl1 = abs. of clouds       (visible band, maximum value)
    
    !          longwave absorptivities (per dp = 10^5 Pa) :
    ! ablwin = abs. of air in "window" band
    ! ablco2 = abs. of air in CO2 band
    ! ablwv1 = abs. of water vapour in H2O band 1 (weak),   for dq = 1 g/kg
    ! ablwv2 = abs. of water vapour in H2O band 2 (strong), for dq = 1 g/kg
    ! ablcl1 = abs. of "thick" clouds in window band (below cloud top) 
    ! ablcl2 = abs. of "thin" upper clouds in window and H2O bands

    real :: solc = 342.0

    real :: albsea = 0.07
    real :: albice = 0.60!0.75
    real :: albsn  = 0.60

    real :: rhcl1  =  0.30
    real :: rhcl2  =  1.00
    real :: qacl   =  0.20
    real :: wpcl   =  0.2
    real :: pmaxcl = 10.0

    real :: clsmax  = 0.60!0.50
    real :: clsminl = 0.15
    real :: gse_s0  = 0.25
    real :: gse_s1  = 0.40

    real :: albcl  =  0.43
    real :: albcls =  0.50

    real :: epssw  =  0.020!0.025
    real :: epslw  =  0.05
    real :: emisfc =  0.98

    real :: absdry =  0.033
    real :: absaer =  0.033
    real :: abswv1 =  0.022
    real :: abswv2 = 15.000

    real :: abscl1 =  0.015
    real :: abscl2 =  0.15

    real :: ablwin =  0.3
    real :: ablco2 =  6.0!5.0
    real :: ablwv1 =  0.7
    real :: ablwv2 = 50.0

    real :: ablcl1 = 12.0
    real :: ablcl2 =  0.6
    real :: ablco2_ref

    ! Time-invariant fields (initial. in radset)
    ! fband  = energy fraction emitted in each LW band = f(T)
    real :: fband(100:400,4)

    ! Zonally-averaged fields for SW/LW scheme (updated in sol_oz)
    ! fsol   = flux of incoming solar radiation
    ! ozone  = flux absorbed by ozone (lower stratos.)
    ! ozupp  = flux absorbed by ozone (upper stratos.)
    ! zenit  = optical depth ratio (function of solar zenith angle)
    ! stratz = stratospheric correction for polar night
    real, dimension(ix*il) :: fsol, ozone, ozupp, zenit, stratz

    ! Radiative properties of the surface (updated in fordate)
    ! alb_l  = daily-mean albedo over land (bare-land + snow)
    ! alb_s  = daily-mean albedo over sea  (open sea + sea ice)
    ! albsfc = combined surface albedo (land + sea)
    ! snowc  = effective snow cover (fraction)
    real, dimension(ix*il) :: alb_l, alb_s, albsfc, snowc

    ! Transmissivity and blackbody rad. (updated in radsw/radlw)
    ! tau2   = transmissivity of atmospheric layers
    ! st4a   = blackbody emission from full and half atmospheric levels
    ! stratc = stratospheric correction term 
    ! flux   = radiative flux in different spectral bands
    real :: tau2(ix*il,kx,4), st4a(ix*il,kx,2), stratc(ix*il,2), flux(ix*il,4)

    ! Radiative properties of clouds (updated in cloud)
    ! qcloud = Equivalent specific humidity of clouds 
    real, dimension(ix*il) :: qcloud, irhtop
end module
