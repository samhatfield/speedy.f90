# Calculates Southern-Oscillation index
# SOI = 10 * (SLP_Tahiti - SLP_Darwin) / sigma
# Where sigma is the standard deviation of the pressure difference

import iris as i
import matplotlib.pyplot as plt
import iris.plot as iplt
import iris.quickplot as qplt

plt.style.use('ggplot')

data_file = '../output/exp_103/attm103.nc'

cube = i.load(data_file, 'mean-sea-level pressure         [hPa]')[0]

tahiti = cube.extract(i.Constraint(latitude=-16.7, longitude=210.0))
darwin = cube.extract(i.Constraint(latitude=-12.989, longitude=131.25))

diff = tahiti - darwin
std = diff.collapsed('time', i.analysis.STD_DEV)
soi = 10 * diff / std
soi.rename('Southern Oscillation Index')

qplt.plot(soi)
plt.show()
