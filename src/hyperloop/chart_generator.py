from hyperloop_sim import HyperloopSim
from openmdao.units.units import convert_units as cu
import pylab
import numpy as np
from time import time
import sys
from os import devnull

def plot(p, x_array, x_varname, y_varnames, x_label, y_label,
        title='HyperloopSim', postprocess_funcs=tuple()):
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
    for i in range(len(y_varnames)):
        pylab.plot(x_array, y_arrays[i], '-', label=y_varnames[i], lw=2,
            c=colors[i % len(colors)], alpha=0.6)
    pylab.tick_params(axis='both', which='major', labelsize=12)
    pylab.xlabel(x_label, fontsize=12)
    pylab.ylabel(y_label, fontsize=12)
    pylab.title(title, fontsize=14)
    pylab.legend(loc='best')
    pylab.show()

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
    p = HyperloopSim.p_factory()
    p['compression_system.comp1.PR_design'] = 1.0
    p['comp1_exit_MN'] = 0.9
    #p['percent_to_bypass'] = 0.01
    #p['pod_MN'] = 0.3
    plot(p, x_array=np.arange(0.1, 0.6, 0.25),
        x_varname='pod_MN',
        y_varnames=('tube_flow.W', 'bypass_W', 'compression_system.inlet.Fl_I:stat:W'),
        x_label='Travel Mach',
        y_label='Mass flow (kg/s)',
        title='OpenMDAO: UW Pod, Rev: 1 Nov 2015',
        postprocess_funcs=(None, None, PostProcess.lbm2kg))