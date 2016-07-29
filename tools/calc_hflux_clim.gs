function writef(args)

* usage run calc_hflux_clim.gs  infile outfile nmonth b

infile=subwrd(args,1)
outfile=subwrd(args,2)
nmonth=subwrd(args,3)
batch=subwrd(args,4)
tmax=nmonth

'reinit'
'open 'infile
'set fwrite 'outfile'.grd'
'set gxout fwrite'


'set x 1 96'
'set y 1 48'
'set z 1'
'set t 1 12'
'define mlshf = ave(lshf,t+0,t='tmax',1yr)' 
'define msshf = ave(sshf,t+0,t='tmax',1yr)' 

'set t 1'
t=1 
while (t<=12)
  'set t 't
  'd mlshf'
  'd msshf'
  t = t+1
endwhile


quit


return




