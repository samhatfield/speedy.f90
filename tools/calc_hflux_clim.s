#!/bin/csh
#-------------------------------------------------------------------
# creating a heat-flux climatology from a speedy experiment
# usage : calc_hflux_clim.s  exp year1 year2   
# E.g. calc_hflux_clim.s 155 1981 1990  	
#-------------------------------------------------------------------	
set EXP=$1		
set YEAR1=$2	
set YEAR2=$3
#
# defining input and output files
#			
@ NYEAR = ($YEAR2 - $YEAR1)  +  1
@ NMONTH = $NYEAR * 12	
echo $NMONTH	
echo $NYEAR			
set OUT_FILE='hflux_'${EXP}'_'${YEAR1}'_'${YEAR2}'_clim'
set IN_FILE='lout/exp_'${EXP}'/attm'${EXP}	
#
# grads-averaging and output
#				
grads -pbc "run calc_hflux_clim.gs  ${IN_FILE} ${OUT_FILE} ${NMONTH}  b"
#

#	
# creating ctl-file
#			
echo 'DSET ^'${OUT_FILE}'.grd' > ${OUT_FILE}'.ctl'
echo 'TITLE '${YEAR1}'/'${YEAR2}' climatology heatflux from speedy experiment '${EXP}   >> ${OUT_FILE}'.ctl'
echo 'UNDEF   9.999E+19'                                                                >> ${OUT_FILE}'.ctl' 
echo '*OPTIONS   SEQUENTIAL'  	                                                        >> ${OUT_FILE}'.ctl'  
echo 'XDEF      96  LINEAR     0.000     3.750'                                         >> ${OUT_FILE}'.ctl' 
echo 'YDEF      48  LEVELS   -87.159   -83.479   -79.777   -76.070   -72.362   -68.652' >> ${OUT_FILE}'.ctl' 
echo '   -64.942   -61.232   -57.521   -53.810   -50.099   -46.389   -42.678   -38.967' >> ${OUT_FILE}'.ctl' 
echo '   -35.256   -31.545   -27.833   -24.122   -20.411   -16.700   -12.989    -9.278' >> ${OUT_FILE}'.ctl' 
echo '    -5.567    -1.856     1.856     5.567     9.278    12.989    16.700    20.411' >> ${OUT_FILE}'.ctl' 
echo '    24.122    27.833    31.545    35.256    38.967    42.678    46.389    50.099' >> ${OUT_FILE}'.ctl' 
echo '    53.810    57.521    61.232    64.942    68.652    72.362    76.070    79.777' >> ${OUT_FILE}'.ctl' 
echo '    83.479    87.159'                                                             >> ${OUT_FILE}'.ctl' 
echo 'ZDEF   1  LEVELS 1013 '                                                           >> ${OUT_FILE}'.ctl' 	
echo 'TDEF 12 LINEAR 1jan'${YEAR1}' 1mo'                                                >> ${OUT_FILE}'.ctl' 
echo 'VARS      2'                                                                      >> ${OUT_FILE}'.ctl' 
echo 'LSHF   0  99  heat flux into land sfc (dw.) [W/m^2] '                             >> ${OUT_FILE}'.ctl'  	     
echo 'SSHF   0  99 heat  flux into sea  sfc (dw.) [W/m^2] '                             >> ${OUT_FILE}'.ctl'  
echo 'ENDVARS '                                                                         >> ${OUT_FILE}'.ctl' 
#
# moving files
#			
ls -lrt ${OUT_FILE}*					
mv ${OUT_FILE}.ctl lout_priv/
mv ${OUT_FILE}.grd lout_priv/	
ls -lrt lout_priv/${OUT_FILE}.*		
#/bin/rm ${INTERPOL_DIR}/fort_code.*	
	