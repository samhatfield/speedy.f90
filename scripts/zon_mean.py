import iris as i
import iris.plot as iplt
import matplotlib.pyplot as plt
from os.path import isfile
from subprocess import call

i.FUTURE.netcdf_promote = True

plt.style.use('ggplot')

def plot_zon_mean(field, title, level=-1):
    orig_f = '../output/exp_year64bits/yyyymmddhh.nc'
    new_f = '../output/exp_year16bits/yyyymmddhh.nc'
    
    orig = i.load(orig_f, field)[0]
    new  = i.load(new_f, field)[0]
    
    # First month (744 hours) is discarded as spin-up  
    if level is -1:
        orig_zon = orig.extract(i.Constraint(time=lambda t: t > 744))\
            .collapsed(['longitude', 'time'], i.analysis.MEAN)
        new_zon  = new.extract(i.Constraint(time=lambda t: t > 744))\
            .collapsed(['longitude', 'time'], i.analysis.MEAN)
    else:
        orig_zon = orig.extract(i.Constraint(generic=level))\
            .extract(i.Constraint(time=lambda t: t > 744))\
            .collapsed(['longitude', 'time'], i.analysis.MEAN)
        new_zon  = new.extract(i.Constraint(generic=level))\
            .extract(i.Constraint(time=lambda t: t > 744))\
            .collapsed(['longitude', 'time'], i.analysis.MEAN)

    fig = plt.figure()

    orig_h, = iplt.plot(orig_zon)
    new_h, = iplt.plot(new_zon)
    
    leg = plt.legend([orig_h, new_h], ['64 bits', '16 bits'], frameon=True)
    rect = leg.get_frame()
    rect.set_linewidth(0.0)
    rect.set_alpha(0.7)

    plt.xlabel('Latitude')
    plt.xlim([-90, 90])

    plt.title(field)
    plt.savefig('../figs/' + title + '.pdf')

plot_zon_mean('Temperature [K]', 'temp_64_16', 850.0)
plot_zon_mean('U-wind [m/s]', 'u_wind_64_16', 850.0)
plot_zon_mean('V-wind [m/s]', 'v_wind_64_16', 850.0)
plot_zon_mean('Specific Humidity [kg/kg]', 'humid_64_16', 850.0)

plt.show()
