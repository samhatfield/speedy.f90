import iris as i
import matplotlib.pyplot as plt
import iris.plot as iplt
from matplotlib import animation
from matplotlib.animation import ArtistAnimation
import cartopy.crs as ccrs
from sys import argv

i.FUTURE.netcdf_promote = True

color_map = plt.get_cmap('inferno')

# Load file
mslp = i.load(argv[1] + '.nc', 'geopotential height at 500 hPa    [m]')[0]
#mslp = mslp.extract(i.Constraint(time=lambda t: t < 216))

# Compute max and min
mslp_avg = mslp.collapsed(['latitude', 'longitude', 'time'], i.analysis.MEAN)
mslp_max = mslp.collapsed(['latitude', 'longitude', 'time'], i.analysis.MAX) - mslp_avg
mslp_min = mslp.collapsed(['latitude', 'longitude', 'time'], i.analysis.MIN) - mslp_avg

# Setup up figure window
fig = plt.figure(figsize=(24,8))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
ax.set_global()
#plt.colorbar()
plt.title('')

# Animate
ims = []
for i, m in enumerate(mslp.slices_over('time')):
    print(i)
    iplt.contourf(m - mslp_avg, 20, cmap=color_map)
    plt.clim([mslp_min.data, mslp_max.data])
    plt.savefig('frames/s{0:03d}.png'.format(i), bbox_inches='tight')
#    ims.append(im.collections)
#im_ani = ArtistAnimation(fig, ims, interval=100)

#im_ani.save('im.mp4', writer='mencoder')
#plt.show()
