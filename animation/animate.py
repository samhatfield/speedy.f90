import iris as i
import matplotlib.pyplot as plt
import iris.plot as iplt
from matplotlib import animation
from matplotlib.animation import ArtistAnimation
import cartopy.crs as ccrs
from sys import argv
from numpy import linspace
from os import makedirs

i.FUTURE.netcdf_promote = True

color_map = plt.get_cmap('inferno')

# Make directory to store frames, if it doesn't exit
makedirs(argv[2], exist_ok=True)

# Load file
mslp = i.load(argv[1] + '.nc', 'Temperature [K]')[0]
mslp = mslp.extract(i.Constraint(generic=925.0))

# Compute max and min
mslp_avg = mslp.collapsed(['latitude', 'longitude', 'time'], i.analysis.MEAN)
mslp_max = mslp.collapsed(['latitude', 'longitude', 'time'], i.analysis.MAX) - mslp_avg
mslp_min = mslp.collapsed(['latitude', 'longitude', 'time'], i.analysis.MIN) - mslp_avg

# Setup up figure window
fig = plt.figure(figsize=(24,8))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
ax.set_global()
plt.title('')

# Animate
for i, m in enumerate(mslp.slices_over('time')):
    print(i)
    cont = iplt.contourf(m - mslp_avg, linspace(mslp_min.data, mslp_max.data, 20), cmap=color_map)
    plt.savefig(argv[2] + '/{0:03d}.png'.format(i), bbox_inches='tight')
    cont.collections[0].remove()
