import iris as i
import iris.plot as iplt
import matplotlib.pyplot as plt
from numpy import linspace

i.FUTURE.netcdf_promote = True

def get_zon_vel(filename):
    orig = i.load_cube(filename, 'U-wind [m/s]')

    # Ignore first month (spin-up)
    orig = orig.extract(i.Constraint(time=lambda t: t > 124))
    return orig.collapsed(['longitude', 'time'], i.analysis.MEAN)

def plot_zon_vel(zon_vel, outname, lims=(-20,44)):
    plt.style.use('ggplot')
    color_map = plt.get_cmap('seismic')

    lim_l, lim_h = lims
    
    plt.figure()

    lat = zon_vel.coord('latitude')
    lev = zon_vel.coord('generic')
    iplt.contourf(zon_vel, linspace(lim_l,lim_h,25), coords=[lat, lev], cmap=color_map)
    plt.gca().invert_yaxis()
    plt.gca().set_yscale('log')

    cb = plt.colorbar()
    cb.set_label('Zonal-time-mean velocity (m/s)')

    plt.xlabel('Latitude')
    plt.ylabel('Pressure (hPa)')

    plt.savefig('../figs/' + outname + '.pdf', bbox_inches='tight')

exp_64bits = '../output/exp_year64bits/yyyymmddhh.nc'
exp_16bits = '../output/exp_year16bits/yyyymmddhh.nc'

plot_zon_vel(get_zon_vel(exp_64bits), '64bits_zon_mean_vel')
plot_zon_vel(get_zon_vel(exp_16bits), '16bits_zon_mean_vel')
plot_zon_vel(get_zon_vel(exp_16bits) - get_zon_vel(exp_64bits), '16bits_64bits_zon_diff', lims=(-20,20))

plt.show()
