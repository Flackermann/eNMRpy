"""
Here, you can find the SpecModel Class for creating a fit model to analyse the eNMR spectra.
You can also use the wrapper function set_peaks() to create a fit model from a matplotlib-based GUI input

You can also Simulate Spectra with Phase shifts.
"""

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
    
    #use built-in method to create pyplot-figure
    m.plot_spec([0], xlim=xlim, **plot_kwargs)
    # adjust the layout for the GUI to minimize bezels
    plt.tight_layout() 
    # insert hints how to correctly use the GUI input
    plt.text(0.95,0.95,'left mouse button: set coordinate \nright mouse button: delete last coordinate\nmiddle mouse button: exit\n',
         horizontalalignment='right',verticalalignment='top', transform=plt.gca().transAxes)
    # activates the input-function delivered by matplotlib
    arr = plt.ginput(n=n, timeout=timeout)
    #returns an array with the selected coordinates
    return arr


def peakpicker(x, y=None, inverse_order=False, width=10, threshold=1e5):
    """
    tool to automatically determine peak positions in an x-y spectrum
    
    x: ppm/frequency-array
        x [optional]: Measurement-object of which the first slice is taken
    y: intensity
    
    inverse_order:
        inverts the order of peaks in the array
        this could be helpful when passing it to make_model
    
    width: scan width. should be even numbered
    
    
    returns: 2D-np.array with [[ppm,intensity],[...]] which can be passed to make_model
    """
    
    #res: step size for scan
    res=1

    if type(x) != np.ndarray:
        y = x.data[0]
        x = x.ppm

    
    if width%2 != 0:
        raise ValueError('width should be an even integer')
#         return
    
    def advindx(arr, min, max):
        """more or less introduces repeating boundaries for the array"""
        if (min < 0) or (max > len(arr)):
            a = arr[min:]
            b = arr[0:max]
            return np.append(a, b)
        else:
            return arr[min:max]
    
    def check_peak(y_out, width=width):
        '''checks for consecutive increase and decrease of values'''
        
        arr = y_out
        
        sel_arr = np.array([])
        for i in range(width):
            # if adjacent values increase, the difference is positive.
            # if adjacent values increase, the difference is negative.
            sel_arr = np.append(sel_arr, (arr[i+1].real - arr[i].real))
        
        # check for increase or decrease
        bool_arr = sel_arr > 0
        
        # if in the first part all values increase and the second part all values decrease
        if (all(bool_arr[:width//2]) == True) and all(not x for x in bool_arr[width//2:]):
            return True
        else:
            return False
    
    peaks = []
    for i in range(0,len(y), res):
        x_out = advindx(x,i-width,i+width)
        y_out = advindx(y,i-width,i+width)

        check = check_peak(y_out, width)
    
        if check and (y[i-width//2].real > threshold):
            # i-width//2 ensures that the right coordinates are given for the respective peak
            peaks.append((x[i-width//2], y[i-width//2].real))
    
    if inverse_order:
        return np.array(peaks)[::-1]
    
    elif not inverse_order:
        return np.array(peaks)


def make_model(peaks, print_params=True):
    '''
    returns a SpecModel()-object with parameters generated from set_peaks()-array
    '''
    
    model = SpecModel(len(peaks))
    # sets the chemical shift vX of each peak
    model.set_initial_values(['v%i'%i for i in range(len(peaks))], [i[0] for i in peaks])
    # sets the estimated amplitude as the peak hight divided by 100
    model.set_initial_values(['a%i'%i for i in range(len(peaks))], [i[1]/100 for i in peaks])
    #shows the resulting parameters in output line
    if print_params:
        model.params.pretty_print()
    return model

# maybe should be handled as static method
def reduce_fitted_phases(Measurementobject, SpecModel):
    """
    Static method for the calculation of reduced phase shift values, where the slope is equal to the
    electrophoretic mobility when plotted against the applied voltage.
    
    This is especially helpful for the comparison of phase shifts obtained for different nuclei
    or under different experimental conditions.
    """
    m = Measurementobject
    for k in ['ph%i'%i for i in range(SpecModel.n)]:
        m.eNMRraw[k+'reduced'] = m.eNMRraw[k]*m.d/m.delta/m.Delta/m.g/m.gamma
        m.eNMRraw[k+'reduced'] -= m.eNMRraw.loc[m.eNMRraw['U / [V]']==0, k+'reduced'][0]

def fit_Measurement(obj_M, obj_S, fixed_parameters=None, plot=False, savepath=None, **plot_kwargs):
    '''
    function to fit a the series of voltage dependent spectra contained in a typical eNMR measurement
    
    obj_M: object of the class eNMR_Measurement
    obj_S: object of the class SpecModel
    fixed_parameters: List of parameters to be fixed after the first fit
    
    **plot_kwargs are passed to SpecModel.fit:
        peak_deconvolution=False, parse_complex='abs','real, or 'imag'
    
    '''
    
    i=0

    if fixed_parameters is not None:
        fp = [] # working list of parameter-objects
        
        # iterates the list of parameters
        for k in obj_S.params.keys(): 
            # if the parameter name p[0] matches any string in fixed_parameters
            if any(k == np.array(fixed_parameters)): 
                # append the parameter to the working list
                fp.append(k)                
        
        
        fig = obj_S.fit(obj_M.ppm, obj_M.data[0], **plot_kwargs)
        print('row 0 fitted including fixed_parameters being varied')
        ph_res = obj_S.get_result_values()
        
        for par in ph_res.keys():
            # saves the results from row 0 in the eNMRraw-DataFrame
            obj_M.eNMRraw.at[0, par] = ph_res[par]
        
        
        for p in fp: # fixes all variables listed in fixed_parameters
            obj_S.params[p].set(obj_S.result.params[p].value)
            obj_S.params[p].set(vary=False)
            print('%s not varied!'%p)
            
        if (plot is True) and (savepath is not None):
            fig.savefig(savepath+'%.1f'%obj_M.eNMRraw.loc[0, obj_M._x_axis]+'.png', dpi=300)
        
        i = 1 #counter set to one for the rest of the spectra to be fitted

    print('start fitting from row %i'%i)
    
    for row in range(i, obj_M.data[:,0].size):
        fig = obj_S.fit(obj_M.ppm, obj_M.data[row], plot=plot, **plot_kwargs)
        ph_res = obj_S.get_result_values()
        
        for par in ph_res.keys():
            #obj_M.eNMRraw.set_value(row, par, ph_res[par])
            obj_M.eNMRraw.at[row, par] = ph_res[par]
            
        if (plot is True) and (savepath is not None):
            fig.savefig(savepath+'%.1f'%obj_M.eNMRraw.loc[row, obj_M._x_axis]+'.png', dpi=300)
            
    #for p in fp: # reset all vary-Values
        #obj_S.params[p].set(vary=True)
    
    
    print('fitting finished')

def drop_errors(df):
    '''
    drops all columns that which keys end with _err --> created from the fitting model
    This function is used in the plot_correlations function
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

def plot_correlations_heatmap(df, method='pearson', without_errors=True, textcolor="#222222", **fig_kwargs):
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
    # transforming the phase angle from degree ° in rad
    ph = (ph/360*2*np.pi)
    amplitude, _lambda = abs(amplitude), abs(_lambda)

    return amplitude*(absorption(x, x0, _lambda)*np.cos(ph)-dispersion(x, x0, _lambda)*np.sin(ph))

def lorentz_imag(x, x0, _lambda, ph, amplitude):
    """
    calculates the imaginary part of the spectrum
    here, the conversion from rad to ° happens
    """
    def dispersion(x, x0, _lambda):
        return -(x-x0)/(_lambda**2+(x-x0)**2)

    def absorption(x, x0, _lambda):
        return _lambda/(_lambda**2+(x-x0)**2)
    
    # transforming the phase angle from degree ° in rad
    ph = (ph/360*2*np.pi)
    amplitude, _lambda = abs(amplitude), abs(_lambda)

    return amplitude*(dispersion(x, x0, _lambda)*np.cos(ph)+absorption(x, x0, _lambda)*np.sin(ph))

def gauss(x, x0, s, amp=1):
    """
    this is just the function of a gaußian distribution with x0 as the center value, s the standard deviation and
    amp as an amplitude which would be 1 in the case of a normal distribution.
    """
    return amp*np.exp(-(x-x0)**2/(2*s**2))/(np.sqrt(np.pi*2)*s)

def makefunc_Lorentz_cmplx(n = 1):
    '''
    returns a function describing the complex lorentzians with real and imaginary part
    n: number of lorentzians contained in the function
    '''
    s = 'lambda x, baseline'
    for i in range(n):
        s += ', v%i, l%i, ph%i, a%i'%(i, i, i, i)
    
    s+=': baseline'
    for i in range(n):
        s += ' + lorentz_real(x, v%i, l%i, ph%i, a%i)'%(i, i, i, i)
        s += ' + 1j*lorentz_imag(x, v%i, l%i, ph%i, a%i)'%(i, i, i, i)
    
    func = eval(s)
    func.__name__="Lorentz Peak Superposition"
    return func

def makefunc_Voigt(n = 1):
    '''
    returns a combination of n Voigt-real as a lambda function
    
    similar to makefunc_Lorentz_cmplx, but with voigt profiles
    '''
    s = 'lambda x, baseline'
    for i in range(n):
        s += ', v%i, l%i, ph%i, a%i, s%i'%(i, i, i, i, i)
    
    s+=': baseline'
    for i in range(n):
        s += ' + lorentz_real(x, v%i, l%i, ph%i, a%i)*gauss(x, v%i, s%i, 1)'%(i, i, i, i, i, i)
        s += ' + 1j*lorentz_imag(x, v%i, l%i, ph%i, a%i)*gauss(x, v%i, s%i, 1)'%(i, i, i, i, i, i)
    
    func = eval(s)
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
        
        # number of peaks in the fitmodel
        self.n = n_peaks
        
        # the chosen model type is parsed with the number of peaks as a the respective function
        # to the lmfit.Model-wrapper
        self.model = lmfit.Model(_func(n_peaks))
        # parameter set is create
        self.params = self.model.make_params(verbose=False)
        
        self.result = None
        
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
        
        # should this be depreciated?
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
            for i in self.params:
                self.params[i].expr = None
                self.params[i].vary = True
            print('model without constraints')
            return
        
        elif reset:
            for i in self.params:
                self.params[i].expr = None
                self.params[i].vary = True
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
        to set the initial value as val, which can be a list as well
        '''

        if type(par) == list:
            i = 0
            for p in par:
                self.params[p].value = val[i]
                i+=1
        else:
            self.params[par].value = val
        
    def set_boundaries(self, par, min=None, max=None):
        """
        explicit name for setting min and max of said parameter
        
        par: parameter key for self.params
        """
        self.params[par].min, self.params[par].max = min,max

    
    def calc_single_peak(self, x, params=None, peak=0):
        '''
        calculates the single peak of a parameter set params
        
        mainly used for plotting purposes to visualize the individual peaks fitted in a superposition
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
        
    def plot_init_spec(self, xdata, single_peaks=False, fig=None):
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
        
        # get the chemical shifts of all peaks one may also use self.n at this point
        # shiftlist = list(filter(lambda x: x[0]=='v', self.params.keys()))
        for i in range(self.n):
            x = self.params['v%i'%i].value
            # the factor 100 comes from the division in the make_model function, where the amplitude was estimated as y/100
            y = self.params['a%i'%i].value*100 
            # string for the annotation
            s = 'P%i'%i
            ax.annotate(s, (x,y), fontsize=16)
            
        return fig
    
    def fit(self, xdata, ydata, plot=False, peak_deconvolution=False, parse_complex='real', figsize=None):
        """
        method to fit a single spectrum consisting of xdata and ydata
        with the previously defined model and parameters
        
        stores the result in self.result
        
        parse_complex: representation of the complex data when plotting
            "real", "imag", or "abs"
        
        returns a matplotlib figure if plot=True
        """
        
        # use the fit method built in the lmfit model object lmfit
        self.result = self.model.fit(ydata, x=xdata, params=self.params)
        
        fig = None
        if plot:
            fig = self.result.plot(parse_complex=parse_complex)[0]
            if figsize != None:
                fig = plt.gcf()
                fig.set_size_inches(figsize)
            
        if peak_deconvolution and plot:
            ax = fig.gca()
            for n in range(self.n):
                ax.plot(xdata, self.calc_single_peak(np.array(xdata), params=self.result.params, peak=n), '--', label='peak '+str(n))

        return fig
    
    #formerly: def get_phasedata(self):
    def get_result_values(self):
        """
        method to return parameters and their errors resulting from the fit (self.result)
        """
        dic = {}
        #for i in range(self.n):
            #dic['ph%i'%i] = self.result.best_values['ph%i'%i]
            #dic['ph%i_err'%i] = self.result.params['ph%i'%i].stderr
            
        for k in self.params.keys():
            dic[k] = self.result.best_values[k]
            dic[k+'_err'] = self.result.params[k].stderr

        return dic
    
    def report(self):
        """shortcut to print self.result.fit_report()"""
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
        
        U = np.linspace(*spec_par['Ulim'], n)
        print(U)
        g = spec_par['g']
        d = spec_par['d']
        Delta = spec_par['Delta']
        delta = spec_par['delta']
        
        for i, mu in enumerate(mobilities):
            self['ph%i'%i] = gamma*Delta*delta*g*(U/d)*mu

        self.U = U
        self.spec_par = spec_par




