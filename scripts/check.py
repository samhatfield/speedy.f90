import iris as i
import iris.plot as iplt
import matplotlib.pyplot as plt
from os.path import isfile
from subprocess import call

plt.style.use('ggplot')

def plot_zon_mean(field, title, level=-1):
    orig_f = '../output/exp_100/attm100.nc'
    new_f = '../output/exp_101/attm101.nc'
    
    orig = i.load(orig_f, field)[0]
    new  = i.load(new_f, field)[0]
    
    if level is -1:
        orig_zon = orig.collapsed(['longitude', 'time'], i.analysis.MEAN)
        new_zon  = new.collapsed(['longitude', 'time'], i.analysis.MEAN)
    else:
        orig_zon = orig.extract(i.Constraint(generic=level)).collapsed(['longitude', 'time'], i.analysis.MEAN)
        new_zon  = new.extract(i.Constraint(generic=level)).collapsed(['longitude', 'time'], i.analysis.MEAN)

    fig = plt.figure()

    orig_h, = iplt.plot(orig_zon)
    new_h, = iplt.plot(new_zon)
    
    leg = plt.legend([orig_h, new_h], ['Original', 'New'], frameon=True)
    rect = leg.get_frame()
    rect.set_linewidth(0.0)
    rect.set_alpha(0.7)

    plt.xlabel('Latitude')
    plt.xlim([-90, 90])

    plt.title(field)
    plt.savefig('../figs/' + title + '.pdf')

call(['sed', '-i', 's/ 365_day_calendar//g', '../output/exp_101/attm101.ctl'])
call(['cdo', '-f', 'nc', 'import_binary', '../output/exp_101/attm101.ctl', '../output/exp_101/attm101.nc'])

plot_zon_mean('geopotential height               [m]', 'geopotential', 925.0)
plot_zon_mean('streamfunction           [10^6 m^2/s]', 'streamfunction', 925.0)
plot_zon_mean('abs. temperature               [degK]', 'temperature', 925.0)
plot_zon_mean('relative humidity                 [%]', 'rhumidity', 925.0)
plot_zon_mean('land-surface temp.             [degK]', 'land_temp')

plt.show()
