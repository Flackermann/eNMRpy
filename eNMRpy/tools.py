from eNMRpy.Measurement.eNMR_Methods import _eNMR_Methods
import pandas as pd
import numpy as np
from io import StringIO

class Load_eNMRpy_data(_eNMR_Methods):
    '''
    this class is used to load .eNMRpy-files which were saved previously
    
    returns: eNMR-Measurement-object which can be used as usual.
    '''
    def __init__(self, path):
        
        if path.split('.')[-1] != 'eNMRpy':
            raise NameError('the given path does not correspond to a .eNMRpy file!')
            return
        
        f = open(path)
        
        dic = eval(f.read().replace('array',''))
            
        for k in dic:
            if '_type_pd.DataFrame' in k:
                setattr(self, k.replace('_type_pd.DataFrame',''), pd.read_json(dic[k]))
            else:
                setattr(self, k, dic[k])
        
        self.data = np.loadtxt(StringIO(self.data), dtype=complex)
        self.ppm = np.loadtxt(StringIO(self.ppm), dtype=float)
        
        # derived instance variables
        self.data_orig = self.data
        self._ppmscale = np.linspace(self._ppm_l, self._ppm_r, self.n_zf_F2)  # np.size(self.data[0,:]))
        #self.ppm = self._ppmscale
        self.fid = self.data
    
    def plot_fid(self):
        raise ValueError('this method is not available in loaded .eNMRpy file! Sorry! please consider the original data')
    
    def proc(self):
        raise ValueError('this method is not available in loaded .eNMRpy file! Sorry! please consider the original data')

def open_measurement(_cls, path):
    '''
    opens .eNMRpy file in path as an instance of _cls
    '''
    f = open(path)
    s = eval(f.read())
    f.close()
        
    obj =  _cls(**(s))
    obj.__dict__.update(s)

    return obj




def relegend(fig, new_labels, **kwargs):
    '''
    Takes a figure with a legend, gives them new labels (as a list) 
    
    **kwargs: ncol, loc etc.
        ncol: number of columns
        loc: location of the legend (see matplotlib documentation)
    '''
    ax = fig.gca()
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, new_labels, **kwargs) 

def calibrate_x_axis(fig, val, majorticks=1, new_xlim=None):
    """
    this function takes a pyplot figure and shifts the x-axis values by val.
    majorticks defines the distance between the major ticklabels.
    """
    ax = fig.gca()
    xlim = ax.get_xlim()
    
    if xlim[0] > xlim[-1]:
        x1, x2 = ax.get_xlim()[::-1]
    else:
        x1, x2 = ax.get_xlim()
        
    ax.set_xticks(np.arange(x1, x2+1, majorticks)-val%1)
    ax.set_xticklabels(['%.1f'%f for f in ax.get_xticks()+val])
    if new_xlim is None:
        ax.set_xlim(xlim)
    else:
        ax.set_xlim(new_xlim)
        
    return spec

def phc_normalized(phc0, phc1):
    """
    nifty function to correct the zero order phase correction
    when phasing 1st order distortions. 
    phc0 -= phc1/(np.pi/2)
    
    returns: (phc0, phc1)
    
    best use in conjunctoin with self.proc()
   
    """
    
    #phc0 -= phc1/(np.pi/2)
    phc0 -= phc1/(np.pi)
    return phc0, phc1
