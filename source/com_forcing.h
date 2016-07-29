C--
C--   /LSMASK/ land-sea masks (initial. in INFORC)
      common /LSMASK/ fmask(ix,il), fmask1(ix,il),  fmask0(ix,il),
     .                fmasko1(ix,il), fmaskl1(ix,il)
C--									
C--   /FORFIX/ Time invariant forcing fields 
C--            (initial. in INFORC, except for phis0 initial. in INVARS)
      common /FORFIX/ phi0(ix,il),
     .                phis0(ix,il), alb0(ix,il)
C--
C--   /FORDAY/ Daily forcing fields (updated in FORDATE)
      common /FORDAY/ alb1(ix,il),
     .                tmoc1(ix,il), 
     .                sst01a(ix,il),
     .                sstm1(ix,il)   

C--
C--   /FORMON/ Monthly-mean forcing fields (initial. in INFORC)
      common /FORMON/ tmp(ix,il,12),
     .                tmoc(ix,il,12),
     .                sst12a(ix,il,12)
