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
