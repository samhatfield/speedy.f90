import iris as i
import matplotlib.pyplot as plt
import iris.plot as iplt
import numpy as np
from matplotlib import animation
import cartopy.crs as ccrs
from sys import argv
from numpy import linspace

color_map = plt.get_cmap('inferno')

field = i.load(argv[1] + '.nc', 'Temperature [K]')[0]
field = field.extract(i.Constraint(time=(int(argv[2])-1)*24))\
             .extract(i.Constraint(generic=850.0))

field_max = field.collapsed(['latitude', 'longitude'], i.analysis.MAX).data
field_min = field.collapsed(['latitude', 'longitude'], i.analysis.MIN).data

fig = plt.figure(figsize=(12,4))
ax = plt.axes(projection=ccrs.PlateCarree())
iplt.contourf(field, linspace(228,296,20),cmap=color_map)
plt.colorbar()
ax.coastlines()
ax.set_global()
plt.title('')

plt.savefig(argv[1] + '.pdf', bbox_inches='tight', pad_inches=0)
plt.show()
