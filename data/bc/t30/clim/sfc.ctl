DSET   ^sfc.grd 
TITLE   constant and annual-mean climate fields from ERA           
UNDEF   -999.                                                            
OPTIONS YREV
XDEF 96 LINEAR 0.0 3.75
YDEF 48 LEVELS  -87.159 -83.479 -79.777 -76.070 -72.362 -68.652 -64.942 -61.232 -57.521 -53.810 -50.099 -46.389 -42.678 -38.967 -35.256 -31.545 -27.833 -24.122 -20.411 -16.700 -12.989  -9.278  -5.567  -1.856   1.856   5.567   9.278  12.989  16.700  20.411  24.122  27.833  31.545  35.256  38.967  42.678  46.389  50.099  53.810  57.521  61.232  64.942  68.652  72.362  76.070  79.777  83.479  87.159
ZDEF       1  LINEAR       1013     1
TDEF      1  LINEAR    jan1981     1mo
VARS       5                                                                   
orog   1  99  orographic height [m]
lsm    1  99  land-sea mask [0-1] 
alb    1  99  annual-mean albedo [0-1] 
vegh   1  99  high vegetation cover [0-1] 
vegl   1  99  low  vegetation cover [0-1] 
ENDVARS  

