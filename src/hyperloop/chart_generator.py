from hyperloop_sim import HyperloopSim
from openmdao.units.units import convert_units as cu
from matplotlib import pyplot, rcParams
import numpy as np
from time import time
import sys
from os import devnull

def plot(p, x_array, x_varname, y_varnames, x_label, y_label,
        title='HyperloopSim', postprocess_funcs=tuple(),
        show=True, filename=''):
    '''
    Runs an OpenMDAO problem for multiple values of x and plots the specified
    results.
    
    Parameters
    ----------
    p : openmdao.core.problem.problem
    x_array : numpy.array
        X values to sample.
    x_varname : str
        OpenMDAO variable name that can be set using p[`x_varname`].
    y_varnames : list of str
        OpenMDAO variable names that can be read using p[`y_varnames`[i]].
    x_label, y_label : str
        Axis labels for graph.
    title : str
        Title of graph.
    postprocess_funcs : tuple of functions
        For each y variable, if the corresponding function in the tuple exists
        and is not equal to None, the function will be called with the y value
        as the only parameter. The returned value will be stored and plotted.
    show : bool
        Whether or not to open plot in a window.
    filename : str
        If specified, the location to save a PNG of the plot.
    '''
    progress_width = 50
    y_arrays = list([] for varname in y_varnames)
    opt_num = 0
    print 'Running optimizations...'
    print ''
    start_time = time()
    for val in x_array:
        elapsed_s = time() - start_time
        progress = float(opt_num) / len(x_array)
        remaining = elapsed_s / progress - elapsed_s if progress > 0 else float('nan')
        sys.stdout.write('[%s] %.2f minutes remaining    \r' %
            ('#' * int(progress * progress_width) + '-' * (progress_width - int(progress * progress_width)),
            float(remaining) / 60))
        sys.stdout.flush()
        p[x_varname] = val
        sys.stdout = open(devnull, 'w')
        p.run()
        sys.stdout.close()
        sys.stdout = sys.__stdout__
        for i in range(len(y_varnames)):
            out = p[y_varnames[i]]
            if len(postprocess_funcs) > i and postprocess_funcs[i] != None:
                out = postprocess_funcs[i](out)
            y_arrays[i].append(out)
        opt_num += 1
    sys.stdout.write('[%s] %.2f minutes elapsed    \r' %
        ('#' * progress_width, float(elapsed_s) / 60))
    sys.stdout.flush()
    print ''
    print ''
    colors = ('r', 'g', 'b', 'c', 'm', 'y', 'k')

    #rcParams['interactive'] = True
    #rcParams['font.family'] = 'Lato'
    #rcParams['font.weight'] = 700
    f = pyplot.figure()
    ax = f.add_subplot(1, 1, 1, frame_on=True, xlabel=x_label,
            ylabel=y_label, title=title)
    ax.tick_params(axis='both', which='major', labelsize=12)
    f.set_size_inches(8, 5)

    for i in range(len(y_varnames)):
        label = y_varnames[i] if (len(postprocess_funcs) < i or
                postprocess_funcs[i] == None) else '* ' + y_varnames[i]
        ax.plot(x_array, y_arrays[i], '-', label=label,
                lw=2, c=colors[i % len(colors)], alpha=0.6)
    ax.legend(loc='best')
    if filename != '':
        f.savefig(filename, dpi=130)
    if show:
        f.show()

class PostProcess():
    @staticmethod
    def lbm2kg(val):
        return cu(val, 'lbm', 'kg')
    @staticmethod
    def degR2degC(val):
        return cu(val, 'degR', 'degC')
    @staticmethod
    def psi2Pa(val):
        return cu(val, 'psi', 'Pa')
    @staticmethod
    def ft2m(val):
        return cu(val, 'ft', 'm')
    @staticmethod
    def ft_s2mph(val):
        return cu(val, 'ft/s', 'mi/h')
    @staticmethod
    def sq_in2sq_m(val):
        return cu(val, 'inch**2', 'm**2')
    @staticmethod
    def invert(val):
        return -val

if __name__ == "__main__":
    p = HyperloopSim.p_factory(inlet_area=0.4735, cross_section=0.8538)

    # Mass flow
    plot(p, x_array=np.arange(0.025, 0.6, 0.025),
        x_varname='pod_MN',
        y_varnames=('tube_flow.W', 'bypass_W', 'compression_system.inlet.Fl_I:stat:W'),
        x_label='Travel Mach',
        y_label='Mass flow (kg/s)',
        title='OpenMDAO: UW Pod, Rev: 10 Nov 2015',
        postprocess_funcs=(None, None, PostProcess.lbm2kg),
        show=True, filename='/Users/brent/Desktop/PyCycle/auto/Mass_Flow.png')
    #p['percent_to_bypass'] = 0.01
    #p['pod_MN'] = 0.3
