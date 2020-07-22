# Copyright Flo unso

import numpy as np
import matplotlib.pyplot as plt
import lmfit
import pandas as pd
from collections import OrderedDict


def set_peaks(m, n=-1, timeout=-1, xlim=None, **plot_kwargs): #,xlim=(0,10)
    """
    m : measurement object
    
    returns: array of tuples with [(v0, a0), ...] for make_model()
    ______________________________________________________________
    should be used in jupyter notebook like this:
    
    %matplotlib
    peaks = set_peaks(m)
    %matplotlib inline
    make_model(peaks)    
    """

    m.plot_spec([0], xlim=xlim, **plot_kwargs)
    plt.tight_layout()
    plt.legend([])
    plt.text(0.95,0.95,'left mouse button: set coordinate \nright mouse button: delete last coordinate\nmiddle mouse button: exit\n',
         horizontalalignment='right',verticalalignment='top', transform=plt.gca().transAxes)
    arr = plt.ginput(n=n, timeout=timeout)

    return arr

def make_model(peaks, print_params=True):
    '''
    returns a SpecModel()-object with parameters generated from set_peaks()-array
    '''
    
    model = SpecModel(len(peaks))
    model.set_initial_values(['v%i'%i for i in range(len(peaks))], [i[0] for i in peaks])
    model.set_initial_values(['a%i'%i for i in range(len(peaks))], [i[1]/100 for i in peaks])
    if print_params:
        model.params.pretty_print()
    return model


def reduce_fitted_phases(Measurementobject, SpecModel):
    
    m = Measurementobject
    for k in ['ph%i'%i for i in range(SpecModel.n)]:
        m.eNMRraw[k+'reduced'] = m.eNMRraw[k]*m.d/m.delta/m.Delta/m.g/m.gamma
        m.eNMRraw[k+'reduced'] -= m.eNMRraw.loc[m.eNMRraw['U / [V]']==0, k+'reduced'][0]


def fit_Measurement(obj_M, obj_S, plot=False, peak_deconvolution=False, savepath=None, fixed_parameters=None):
    '''
    obj_M: object of the class eNMR_Measurement
    obj_S: object of the class SpecModel
    fixed_parameters: List of parameters to be fixed after the first fit
    '''
    
    fp = []
    
    
    if fixed_parameters is None:
        i = 0
        
    elif fixed_parameters is not None:
        #finds the matching parameter keys
        if len(fixed_parameters[0]) == 2:
            fp = fixed_parameters
        else:
            for i in range(len(fixed_parameters)):
                for p in obj_S.params:
                    if p[0] == fixed_parameters[i]:
                        fp.append(p)
                    
        i = 1 #counter set to one for the rest of the spectra to be fitted
        fig = obj_S.fit(obj_M.ppm, obj_M.data[0], plot=plot, peak_deconvolution=peak_deconvolution)
        ph_res = obj_S.get_phasedata()
        
        for par in ph_res.keys():
            #obj_M.eNMRraw.set_value(row, par, ph_res[par])
            obj_M.eNMRraw.at[0, par] = ph_res[par]
            
        for p in fp: # fixes all variables listed in fixed_parameters
            obj_S.params[p].set(vary=False)
            print('%s not varied!'%p)
            
        if (plot is True) and (savepath is not None):
            fig.savefig(savepath+'%.1f'%obj_M.eNMRraw.loc[0, obj_M._x_axis]+'.png', dpi=300)
    
    print('start fitting')
    
    for row in range(i, obj_M.data[:,0].size):
        fig = obj_S.fit(obj_M.ppm, obj_M.data[row], plot=plot, peak_deconvolution=peak_deconvolution)
        ph_res = obj_S.get_phasedata()
        
        for par in ph_res.keys():
            #obj_M.eNMRraw.set_value(row, par, ph_res[par])
            obj_M.eNMRraw.at[row, par] = ph_res[par]
            
        if (plot is True) and (savepath is not None):
            fig.savefig(savepath+'%.1f'%obj_M.eNMRraw.loc[row, obj_M._x_axis]+'.png', dpi=300)
            
    for p in fp: # reset all vary-Values
        obj_S.params[p].set(vary=True)
    
    
    print('fitting finished')

def drop_errors(df):
    '''
    drops all columns that which keys end with _err --> created from the fitting model
    '''

    df.keys()
    #drop the vc, Voltag, Gradient and outlier columns
    try:
        sel = df.drop(['vc', 'U / [V]', 'g in T/m', 'outlier'], axis=1)
    except KeyError:
        try:
            sel = df.drop(['vd', 'U / [V]', 'g in T/m', 'outlier'], axis=1)
        except:
            print('no vd or vc found, no standard parameters dropped')
            sel = df
            #sel = df.drop(['vd', 'U / [V]', 'g in T/m', 'outlier'], axis=1)
    # drop all error-columns-boolean
    _bool = np.array([k[-4:] != '_err' for k in sel.keys()])
    # elect only the non-error columns
    sel = sel[np.array(sel.keys()[_bool])]
    return sel

def plot_correlations(df, method='pearson', without_errors=True, textcolor="#222222", **fig_kwargs):
    """
    correlation coefficients plot (heatmap) for any pandas DataFrame
    
    method:
        pearson, kendall, spearman
    """
    if without_errors:
        df = drop_errors(df)
    corr = df.corr(method=method)
    columns = corr.keys()
    indices = np.array(corr.index)

    fig, ax = plt.subplots(**fig_kwargs)
    im = ax.imshow(corr, cmap='Spectral', vmin=-1, vmax=1)

    # We want to show all ticks...
    ax.set_xticks(np.arange(len(columns)))
    ax.set_yticks(np.arange(len(indices)))
    # ... and label them with the respective list entries
    ax.set_xticklabels(columns)
    ax.set_yticklabels(indices)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    for i in range(len(indices)):
        for j in range(len(columns)):
            text = ax.text(j, i, '%.3f'%corr.iloc[i, j],
                           ha="center", va="center", color=textcolor)
    plt.colorbar(im)
    fig.tight_layout()
    
    return fig, corr


def lorentz_real(x, x0, _lambda, ph, amplitude):
    """
    calculates the real part of the spectrum
    here, the conversion from rad to ° happens
    """
    def dispersion(x, x0, _lambda):
        return -(x-x0)/(_lambda**2+(x-x0)**2)

    def absorption(x, x0, _lambda):
        return _lambda/(_lambda**2+(x-x0)**2)
    #Umrechnung von ° in rad
    ph = (ph/360*2*np.pi)
    amplitude, _lambda = abs(amplitude), abs(_lambda)

    return amplitude*(absorption(x, x0, _lambda)*np.cos(ph)-dispersion(x, x0, _lambda)*np.sin(ph))

def lorentz_imag(x, x0, _lambda, ph, amplitude):
    """
    calculates the real part of the spectrum
    here, the conversion from rad to ° happens
    """
    def dispersion(x, x0, _lambda):
        return -(x-x0)/(_lambda**2+(x-x0)**2)

    def absorption(x, x0, _lambda):
        return _lambda/(_lambda**2+(x-x0)**2)
    
    #Umrechnung von ° in rad
    ph = (ph/360*2*np.pi)
    amplitude, _lambda = abs(amplitude), abs(_lambda)

    return amplitude*(dispersion(x, x0, _lambda)*np.cos(ph)+absorption(x, x0, _lambda)*np.sin(ph))

def gauss(x, x0, s, amp):
    return amp*np.exp(-(x-x0)**2/(2*s**2))/(np.sqrt(np.pi*2)*s)

#def makefunc_Lorentz(n = 1):
    #'''
    #returns a combination of n lorentz-real as a lambda function
    #'''
    #s = 'lambda x, baseline'
    #for i in range(n):
        #s += ', v%i, l%i, ph%i, a%i'%(i, i, i, i)
    
    #s+=': baseline'
    #for i in range(n):
        #s += ' + lorentz_real(x, v%i, l%i, ph%i, a%i)'%(i, i, i, i)
    
    #func = eval(s)
    ##func.__name__="Flo's cooles Spektrum"
    #func.__name__="Lorentz Peak Superposition"
    #return func

def makefunc_Lorentz_cmplx(n = 1):
    '''
    returns a combination of n lorentz-real as a lambda function
    '''
    s = 'lambda x, baseline'
    for i in range(n):
        s += ', v%i, l%i, ph%i, a%i'%(i, i, i, i)
    
    s+=': baseline'
    for i in range(n):
        s += ' + lorentz_real(x, v%i, l%i, ph%i, a%i)'%(i, i, i, i)
        s += ' + 1j*lorentz_imag(x, v%i, l%i, ph%i, a%i)'%(i, i, i, i)
    
    func = eval(s)
    #func.__name__="Flo's cooles Spektrum"
    func.__name__="Lorentz Peak Superposition"
    return func

def makefunc_Voigt(n = 1):
    '''
    returns a combination of n Voigt-real as a lambda function
    '''
    s = 'lambda x, baseline'
    for i in range(n):
        s += ', v%i, l%i, ph%i, a%i, s%i'%(i, i, i, i, i)
    
    s+=': baseline'
    for i in range(n):
        s += ' + lorentz_real(x, v%i, l%i, ph%i, a%i)*gauss(x, v%i, s%i, 1)'%(i, i, i, i, i, i)
        s += ' + 1j*lorentz_imag(x, v%i, l%i, ph%i, a%i)*gauss(x, v%i, s%i, 1)'%(i, i, i, i, i, i)
    
    func = eval(s)
    #func.__name__="Flo's cooles Spektrum"
    func.__name__="Voigt Peak Superposition"
    return func

class SpecModel(object):
    '''
    This class creates a lmfit Model from s function of n_peaks lorentz-distributions
    parameter explanation
        an = amplitude n
        ln = peak broadness of n
        phn = phase angle n
        vn = chemical shift v of n
        baseline = basline of the whole spectrum
        model = 'Lorentz' or 'Voigt'
    '''
    
    def __init__(self, n_peaks, verbose=False, model='Lorentz'):
        
        self.mkey = model
        _func = {'Lorentz': makefunc_Lorentz_cmplx,
                 'Voigt' : makefunc_Voigt}[model]
        
        self.n = n_peaks
        self.model = lmfit.Model(_func(n_peaks))
        self.params = self.model.make_params(verbose=False)
        
        # This loop initializes all parameters by 1 and sets standard boundaries depending on the respective parameterset
        for k in self.params.keys():
            self.set_initial_values(k, 1)
            if k[0] == 'a':
                self.set_boundaries(k, 0, np.inf)
                self.set_initial_values(k, 1e6)
            elif k[0] == 'l':
                self.set_boundaries(k, 0, np.inf)
                self.set_initial_values(k, 1e-2)
            elif k[0] == 'p':
                self.set_boundaries(k, -180, 180)
            elif k[0] =='s':
                self.set_boundaries(k, 0, np.inf)
                self.set_initial_values(k, 1e-2)
            else:
                pass

        if verbose:
            self.params.pretty_print()

    def set_mathematical_constraints(self, expr=None, reset=True):
        """
        expr: mathematical expression without whitespace
        
        reset: will reset all mathematical constraints
        
        example 1:
            ph0 should be always equal to ph1
            set_mathematical_constraint('ph0=ph1')
        example 2:
            ph0 should be always ph1 - 90 degree
            set_mathematical_constraint('ph0=ph1-90')
        """
        
        if (expr == None) and (reset):
            print('model without constraints')
            return
        
        elif reset:
            for i in self.params:
                self.params[i].expr = None
            print('constraints were reset before new assignment')
        
        if type(expr) == list:
            for e in expr:
                a = e.split('=')
                self.params[a[0]].expr = a[1]
        else:
            a = expr.split('=')
            self.params[a[0]].expr = a[1]
        
         
        
    def set_initial_values(self, par, val):
        '''
        takes single parameter par or list of parameters
        '''
        #if type(par) == list:
            ##i = 0
            #for i in range(self.n):
                #self.params[par[i]].value = val[i]
                ##i+=1
        if type(par) == list:
            i = 0
            for p in par:
                self.params[p].value = val[i]
                i+=1
        else:
            self.params[par].value = val
        
    def set_boundaries(self, par, min, max):
        self.params[par].min, self.params[par].max = min,max
        #print(self.params[par])
    
    def calc_single_peak(self, x, params=None, peak=0):
        '''
        calculates the single peak of a parameter set
        '''
        if params is None:
            params = self.params
        
        def make_singlepeak_params(params, peak=0):
            '''
            creates a dictionary containing all parameters of the respective peak named as 
            0-Peak to pass it to the SpecModel evaluation function of a single peak
            '''
            dic = {k[:-1]+'0': params[k]for k in params.keys() if k[-1] == str(peak)}
            dic['baseline'] = params['baseline']
            return dic

        single_peak = SpecModel(1, model=self.mkey)
        s0par = make_singlepeak_params(params, peak)
        return single_peak.model.eval(s0par, x=x)
        
    def plot_init_spec(self, xdata, single_peaks=True, fig=None):
        '''
        plots the spectrum for the initialized parameter set
        '''
        if fig is None:
            fig, ax = plt.subplots()
        else:
            ax = fig.gca()
        
        if single_peaks:
            for n in range(self.n):
                ax.plot(xdata, self.calc_single_peak(np.array(xdata), params=self.params, peak=n).real, '--', label='peak '+str(n))
        ax.plot(np.array(xdata), self.model.eval(x=np.array(xdata), params=self.params).real, c='k', label='start parameters')
        ax.legend()
        
        shiftlist = list(filter(lambda x: x[0]=='v', self.params.keys()))
        for i in range(len(shiftlist)):   
            x = self.params['v%i'%i].value
            y = self.params['a%i'%i].value*100
            s = 'P%i'%i
            ax.annotate(s, (x,y), fontsize=16)
            
        return fig
    
    def fit(self, xdata, ydata, plot=False, peak_deconvolution=False):
        
        self.result = self.model.fit(ydata, x=xdata, params=self.params)
        
        fig = None
        if plot:
            fig = self.result.plot()[0]
        
        if peak_deconvolution and plot:
            ax = fig.gca()
            for n in range(self.n):
                ax.plot(xdata, self.calc_single_peak(np.array(xdata), params=self.result.params, peak=n), '--', label='peak '+str(n))
        #if plot:# =='fit':
            #fig = self.result.plot_fit()
        return fig
    
    def get_phasedata(self):
        dic = {}
        #for i in range(self.n):
            #dic['ph%i'%i] = self.result.best_values['ph%i'%i]
            #dic['ph%i_err'%i] = self.result.params['ph%i'%i].stderr
            
        for k in self.params.keys():
            dic[k] = self.result.best_values[k]
            dic[k+'_err'] = self.result.params[k].stderr

        return dic
    
    def report(self):
        print(self.result.fit_report())
        
    def reassign_parameter(self, par, best_values_dic, vary=False):
        '''
        reassigns the initial value for the parameter by taking the value from the best-fit dictionary (usually obtained form a test-fit)
        
        vary: sets if the parameter should be varied in further fittings or not.
        '''
        self.params[par].set(best_values_dic[par], vary=vary)
        print(par, '=', best_values_dic[par], ';', 'vary =', vary)


class SpecSim(object):
    '''
    This class lets you simulate peaks with given properties in order compare e.g. MOSY and fitting results to yours if Peaks are overlapping.
    '''
    def __init__(self, n_Peaks, ppm=None, xlim=None, n_points=2**11, model="Lorentz", fit_params=None):
        """
        model: 'Lorentz' or 'Voigt'
        """
        
        _func = {'Lorentz': makefunc_Lorentz_cmplx,
            'Voigt' : makefunc_Voigt}[model]
        self.data = None
        self.func = _func(n_Peaks)
        self.pkeys = self.func.__code__.co_varnames
        self.pdic = OrderedDict({k: 0.1 for k in self.pkeys})
        if fit_params is not None:
            try:
                for p in fit_params.keys():
                    self.pdic[p] = fit_params[p].value
            except:
                print('This fit_params are no lmfit-Parameters object')

        print('func_params', self.pkeys)
        
        if ppm is not None:
            self.ppm = ppm
        elif xlim is not None:
            self.ppm = np.linspace(*xlim, n_points)
        else:
            self.ppm = None
        self.pdic['x'] = self.ppm
        
        self.eNMRraw = pd.DataFrame()
        self.eNMRraw['U / [V]'] = None
        self.params = None

        
    def set_params(self, par, val):
        '''
        sets the parameter par as value value in the parameters dictionary self.pdic
        
        par:
            string or list of strings representing the respective keyword for the parameters
        var: single value or list of values. This can be anything, also numpy arrays.
        '''
        if type(par) == list:
            i=0
            for p in par:
                self.pdic[p] = val[i]
                i += 1
        else:
            self.pdic[par] = val
            
    def add_noise(self, scale=1):
        '''
        adds normal distributed noise to the spectrum
        '''
        self.data += scale*np.random.standard_normal(self.data.shape)
        

    def calc_spec(self, noise=0):
        '''
        calculates the spectrum from the parameters dictionary self.pdic
        '''
        self.data = self.func(*self.pdic.values())
        self.add_noise(noise)
        return self.data
    
    def calc_spec_series(self, params, vlist=None, noise=0):
        '''
        calculates a series of spectra with the given lists.
        if par is a list, then vlist needs to be 2-dimensional and rectangular.
        '''
        if vlist is None:
            par = params.keys()
            vlist = params.values()
            #self.Ulim = params.spec_par['Ulim']

        else:
            vlist = np.array(vlist)
        
        if (vlist.dtype == 'int64') or (vlist.dtype == 'float64'):
            if type(par) != list:
                self.set_params(par, vlist[0])
                self.data = np.array([self.calc_spec(noise)])
                for i in range(1, len(vlist)):
                    self.set_params(par, vlist[i])                
                    self.data = np.append(self.data, np.array([self.calc_spec(noise)]), axis=0)

            elif type(par) == list:
                # set all first parameter values
                for i in range(len(par)):
                    self.set_params(par[i], vlist[i][0])
                # create first row of dataset so rest can be appended in the following loop
                self.data = np.array([self.calc_spec(noise)])
                for j in range(1, len(vlist[0])):
                    for i in range(len(par)):
                        self.set_params(par[i], vlist[i][j])
                    self.data = np.append(self.data, np.array([self.calc_spec(noise)]), axis=0)
        else:
            print('error: vlist is not in a rectangular shape')

        #Instanzvariablen Updaten
        self.params = params
        self.eNMRraw['U / [V]'] = params.U

        _d = pd.DataFrame(self.data)
        self.eNMRraw['data'] = None
        for i in range(len(self.data[:,0])):
            self.eNMRraw.at[i, 'data'] = _d[i]#test[i]
        
        for p in params.keys():
            self.eNMRraw[p] = params[p]
        
        return self.data
    
    def plot(self, text_pos=None):
        fig, ax = plt.subplots()
        if self.data.ndim == 1:
            ax.plot(self.ppm, self.data.real)
        else:    
            for n in range(len(self.data[:,0])):
                ax.plot(self.ppm, self.data[n].real, label=n)
        ax.legend()
        if text_pos is not None:
            ax.text(*text_pos, self.params.spec_par)
        return fig
    
    
class ParameterVariation(OrderedDict):
    '''
    OrderedDictionary with extra methods for the creation of datasets
    '''
    def __init__(self):
        self.spec_par = {'Delta': None,
                         'delta': None,
                         'g': None,
                         'Ulim':None,
                         'NUC':None,
                         }
        self.U = None
        pass
    def values(self):
        return np.array(list(super().values()))
    
    def keys(self):
        return list(super().keys())
    
    def set_par(self, par, val):
        if type(par) != list:
            self[par] = val
        else:
            for i, p in enumerate(par):
                self[p] = val[i]
    
    def set_ph_from_mobility(self, mobilities, n, spec_par=None):
        '''
        mobilities:
            must be an array. mobilities are sorted according
<Parameter 'l2', 1, bounds=[0:inf]>
<Parameter 'ph2', 1, bounds=[-180:180]>
<Parameter 'a2', 1, bounds=[0:inf]>
Name         Value      Min      Max   Stderr     Vary     Expr Brute_Step
a0           4e+05        0      inf to the peaks
        
        spec_par: dictionary with the spectroscopic parameters in SI units:
            Delta in s
            delta in s
            NUC: 1H, 7Li or 19F
            g in T/m
            Ulim: Voltage limits in V
            d: electrode distance
        '''
        if spec_par is None:
            spec_par = self.spec_par
        
        #gamma in MHz/T
        #gamma = {"7Li": 5957443266,
                #"1H": 15326621020,
                #"19F": 14415618125}[spec_par['NUC']]
        # gamma in rad/Ts
        gamma = {'1H':26.7513e7,
                '7Li': 10.3962e7,
                '19F': 25.1662e7}[spec_par['NUC']]
        # Umrechnung von rad in °
        gamma = gamma/2/np.pi*360
        print(gamma)
        
        U = np.linspace(*spec_par['Ulim'], n)
        print(U)
        g = spec_par['g']
        d = spec_par['d']
        Delta = spec_par['Delta']
        delta = spec_par['delta']
        
        print(g, d, Delta, delta)
        for i, mu in enumerate(mobilities):
            self['ph%i'%i] = gamma*Delta*delta*g*(U/d)*mu
        print(self)
        self.U = U
        self.spec_par = spec_par
        pass




