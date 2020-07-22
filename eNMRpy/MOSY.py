import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

class MOSY(object):
    """
    takes an eNMR Measurement obj and yields as MOSY object for MOSY processing and depiction.
    the Measurement object needs to be processed (fouriertransformed in F2) before passing it to this class
    """
    def __init__(self, obj):
        self.obj = obj
        self.eNMRraw = obj.eNMRraw
        self.eNMRraw['data'] = None
        self.eNMRraw['fid'] = None
        for row in range(len(obj.data[:,0])):
            self.eNMRraw.set_value(row, 'data', pd.Series(obj.data[row,:]))
            #self.eNMRraw.set_value(row, 'fid', pd.Series(obj.fid[row,:]))
        data_sorted = np.array(obj.eNMRraw['data'].tolist())
        self.data = data_sorted.astype('complex').copy()
        self.TD1 = len(self.data[:,0])
        self.mscale = None
        self.mX = None
        self.mY = None
#         self.mZ = None
        
    def zerofilling_old(self, n=128, dimension=0):
        """
        n: number of total datapoints along the F1-dimension (Voltage)
        dimension 0 = F1
        dimension 1 = F2
        """
        if dimension == 0:
            self.data = np.concatenate((self.data, 
                                [[0 for x in range(len(self.data[0,:]))] for y in range(n-len(self.data[:,0]))]))
        if dimension == 1:
            self.data = np.concatenate((self.data, 
                                [[0 for x in range(n-len(self.data[0,:]))] for y in range(len(self.data[:,0]))]), axis=1)
        print('got it!')
        
    def zerofilling(self, n=128, dimension=0):
        """
        n: number of total datapoints along the F1-dimension (Voltage)
        dimension 0 = F1
        dimension 1 = F2
        """
        print('zero filling started!')
        
        if dimension == 0:
            old_shape = self.data.shape
            new_shape = old_shape
            new_shape = n, old_shape[1]
            print(old_shape, new_shape)
            data_new = np.zeros(new_shape, dtype=complex)
            for i in range(old_shape[0]):
                data_new[i,:] = self.data[i,:]
            self.data = data_new
            del(data_new)
            
        if dimension == 1:
        
            old_shape = self.data.shape
            new_shape = old_shape
            new_shape = old_shape[0], n
            print(old_shape, new_shape)
            data_new = np.zeros(new_shape, dtype=complex)
            for i in range(old_shape[1]):
                data_new[:,i] = self.data[:,i]
            self.data = data_new
            del(data_new)
        
        print('zero filling finished!')

    def fft_F1 (self):
        """
        Fourier Transformation from nmrglue.proc_base.fft)
        along the F1 dimension
        """
        
        from nmrglue.proc_base import fft
        
        #data_temp = np.zeros(self.data.shape)
        for n in range(len(self.data[0,:])):
            self.data[:,n] = fft(self.data[:,n])
        print("done")
    
    def calc_MOSY(self, u_max = None, n_zf=2**12, mobility_scale=True, include_0V=True, electrode_distance=2.2e-2, old_zf=False):
        
        '''
        Calculates the MOSY according to the States-Haberkorn-Ruben method described for MOSY applications by AUTOR NENNEN!
        
        u_max:
            maximum Voltage to be used for processing.
        n_zf:
            number of zerofilling for the F1 dimension
        mobility_scale:
            prints the converted mobilty y-scale if True
        include_0V:
            include the 0V spectrum or not.
        '''
        x = self.obj._x_axis
        
        if u_max is None:
            u_max = 10000
        
        if include_0V:
            positive = self.eNMRraw[(self.eNMRraw[x] >= 0)
                                    &(self.eNMRraw[x] <= u_max)].sort_values(x)
            negative = self.eNMRraw[(self.eNMRraw[x] <= 0)
                                    &(self.eNMRraw[x] >= -u_max)].sort_values(x, ascending=False)
        else:
            positive = self.eNMRraw[self.eNMRraw[x] >0].sort_values(x)
            negative = self.eNMRraw[self.eNMRraw[x] <0].sort_values(x, ascending=False)

        
        #dataset conversion for the SHR-method
        #if include_0V:
            #positive = self.eNMRraw[(self.eNMRraw['U / [V]'] >= 0)
                                    #&(self.eNMRraw['U / [V]'] <= u_max)].sort_values('U / [V]')
            #negative = self.eNMRraw[(self.eNMRraw['U / [V]'] <= 0)
                                    #&(self.eNMRraw['U / [V]'] >= -u_max)].sort_values('U / [V]', ascending=False)
        #else:
            #positive = self.eNMRraw[self.eNMRraw['U / [V]'] >0].sort_values('U / [V]')
            #negative = self.eNMRraw[self.eNMRraw['U / [V]'] <0].sort_values('U / [V]', ascending=False)

        SHR_real = np.zeros((len(positive['data']), len(positive['data'].iloc[0])))
        SHR_imag = np.zeros((len(positive['data']), len(positive['data'].iloc[0])))

        for n in range(1, len(positive['data'])):
            SHR_real[n,:] = positive['data'].iloc[n].real + negative['data'].iloc[n].real
            SHR_imag[n,:] = positive['data'].iloc[n].imag - negative['data'].iloc[n].imag

        SHR = SHR_real + SHR_imag*1j
        del(SHR_real)
        del(SHR_imag)
        self.data = SHR
        del(SHR)
        if old_zf:
            self.zerofilling_old(n=n_zf)
        else:
            self.zerofilling(n=n_zf)
        self.fft_F1()
     
        X = np.array([self.obj.ppm for y in range(len(self.data[:,0]))])
        Y = np.array([[y for x in range(len(self.data[0,:]))] for y in range(len(self.data[:,0]))])
        
        Y = ((0.5*n_zf-Y)*self.TD1/n_zf*360)/self.TD1 # conversion to ° phasechange per Uink
        Y = Y/self.obj.uInk  # ° phasechange per Volt
        
        if mobility_scale:
            #Y = (Y*self.obj.d)/(self.obj.gamma*self.obj.delta*self.obj.Delta*self.obj.g) # old version with self.d
            Y = (Y*electrode_distance)/(self.obj.gamma*self.obj.delta*self.obj.Delta*self.obj.g)
        
        self.mX = X
        del(X)
        self.mY = -Y
        del(Y)
    
    
    def plot_MOSY(self, xlim=None, ylim=(-1e-8, 1e-8), yscale='linear',
                  save=False, savepath='',
                  tight_layout=False, dpi=300, y_autoscale=True,
                  h_ratio=7, w_ratio=8, figsize=(5,5),
                  latex_backend=False, hspace=0, wspace=0,
                  **kwargs):
        
        from matplotlib import gridspec
        """
        yscale:
            .. ACCEPTS: [ 'linear' | 'log' | 'symlog' | 'logit' | ... ]
        **kwargs:
            are passed to the .contour()
            for example:
                levels: list of values to be drawn as height lines
                
        :returns: figure
        """
        fig= plt.figure(figsize=figsize)
        
        gs = gridspec.GridSpec(2, 2, width_ratios=[1, w_ratio], height_ratios=[1,h_ratio]) 

        mitte = fig.add_subplot(gs[3])
        oben = fig.add_subplot(gs[1], sharex=mitte, frameon=False)
        links = fig.add_subplot(gs[2], sharey=mitte, frameon=False)

        #fig.subplots_adjust(hspace=0, wspace=0)

        mitte.yaxis.tick_right()
        mitte.contour(self.mX, self.mY, self.data.real, **kwargs)
        mitte.set_ylim(ylim)
        mitte.set_xlim(xlim)
        links.set_xticks([])

        # calculation and plotting of the projections taking the maxima
        links.plot(np.amax(self.data, axis=1), self.mY[:,0],  'k')
        links.set_xlim(links.get_xlim()[::-1])
        oben.plot(self.mX[0], np.amax(self.data, axis=0), 'k')
        
        #axis format
        oben.tick_params(bottom='off')
        links.tick_params(left='off')
        oben.set_yticks([])
        plt.setp(oben.get_xticklabels(), visible=False);
        plt.setp(links.get_yticklabels(), visible=False);
        plt.setp(links.get_yaxis(), visible=False);
        if latex_backend:
            mitte.set_ylabel(r'$\mu\;/\;(\textrm{m}^2 \textrm{V}^{-1} \textrm{s}^{-1})$')
        else:
            mitte.set_ylabel(r'$\mu$ / (m$^2$V$^{-1}$s$^{-1})$')
        mitte.yaxis.set_label_position("right")
        mitte.set_xlabel(r'$\delta$ / ppm')
        mitte.set_yscale(yscale)
        mitte.ticklabel_format(axis='y', style='sci', scilimits=(0,2))
        mitte.grid()
        
        fig.subplots_adjust(hspace=hspace, wspace=wspace)
        
        def autoscale_y(bx, margin=0.1):
            """This function rescales the y-axis based on the data that is visible given the current xlim of the axis.
            ax -- a matplotlib axes object
            margin -- the fraction of the total height of the y-data to pad the upper and lower ylims"""

            def get_bottom_top(line):
                xd = line.get_xdata()
                yd = line.get_ydata()
                lo,hi = bx.get_xlim()
                try:
                    y_displayed = yd[((xd>lo) & (xd<hi))]
                    h = np.max(y_displayed) - np.min(y_displayed)
                except:
                    y_displayed = yd[((xd<lo) & (xd>hi))]
                    h = np.max(y_displayed) - np.min(y_displayed)
                
                bot = np.min(y_displayed)#-margin*h
                top = np.max(y_displayed)#+margin*h
                print(bot, top, h)
                return bot,top

            lines = bx.get_lines()
            bot,top = np.inf, -np.inf

            for line in lines:
                new_bot, new_top = get_bottom_top(line)
                if new_bot < bot: bot = new_bot
                if new_top > top: top = new_top

            bx.set_ylim(bot,top)
            
        if y_autoscale:
            autoscale_y(oben, xlim)
        
        if tight_layout:
            fig.tight_layout()
        
        if save:
            fig.savefig(savepath, dpi=dpi)

        return fig
    
    def plot_slices_F1(self, ppm, xlim=None, scaling=None, normalize=False, vline_0=False, annotate=False, annotate_pos=None, legend_loc=None, latex_backend=False, colors=None, figsize=None):
        """
        plots slices in F1-direction

        ppm:
            list of chemical shift-slices to be plotted
        xmin, xmax:
            limits for the x-axis
        scaling:
            list of scaling factors for spectra to adjust intensities
        normalize:
            normalization of each peak to 1
        vline_0:
            prints a vertical line at µ = 0
        annotate:
            automatically annotates the maxima of the slices and displays the respective mobility
        colors:
            list of colors/colorcodes as strings
        
        :returns: figure
        """

        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)
        # getting the standard color cycle
        prop_cycle = plt.rcParams['axes.prop_cycle']
        if colors is None:
            colors = prop_cycle.by_key()['color']

        if scaling is None:
            scaling = np.ones((len(ppm)))
        i=0
        # setting the initial textposition depending on the number of slices
        xpostext = 0.5-((len(ppm)-1)*0.1-0.1)
        if xpostext < 0.2:
            xpostext = 0.2
        ypostext = 0.2
        
        for n in ppm:
            pos = np.searchsorted(self.obj.ppm, n)
            x = self.mY[:,pos]
            if normalize:
                y = scaling[i]*self.data[:,pos].real/self.data[:,pos].max()
                ax.plot(x, y, label=r'$\delta=%.2f$'%n, c=colors[i])
                ax.set_ylabel('normalized intensity')
            else:
                y = scaling[i]*self.data[:,pos].real
                ax.plot(x, y, label=r'$\delta=%.2f$'%n, c=colors[i])
                ax.set_ylabel('intensity / a.u.')
            
            if annotate:
                xmax = x[np.argmax(y)]
                ymax = y.max()
                text = "$\mu$=%.2E"%(xmax)
                bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
                arrowprops=dict(arrowstyle="->", color=colors[i])#, connectionstyle="angle,angleA=0,angleB=60")
                kw = dict(xycoords='data',textcoords="axes fraction",
                arrowprops=arrowprops, bbox=bbox_props, ha="right", va="top")
                if type(annotate_pos) == list:
                    xpostext, ypostext = annotate_pos[i]
                elif type(annotate_pos) == tuple:
                    xpostext += annotate[0]
                    ypostext += annotate[1]
                ax.annotate(text, xy=(xmax, ymax), xytext=(xpostext,ypostext), **kw)
                #placing of the textboxes
                if xpostext >= 1:
                    xpostext = 0.2
                    ypostext += 0.1
                else:
                    xpostext += 0.2
            i += 1

        if vline_0:
            ax.vlines(0, *ax.get_ylim(), linestyles='dotted')
        if latex_backend:
            ax.set_xlabel('$\mu$\;/\;(m$^2$V$^{-1}$s$^{-1})$')
        else:
            ax.set_xlabel('$\mu$ / (m$^2$V$^{-1}$s$^{-1})$')
        ax.set_xlim(xlim)
        ax.legend(loc=legend_loc)
        ax.ticklabel_format(style='sci', scilimits=(0,3))
        return fig
    
    def export_data(self, path):
        """
        saves the X, Y and Z-data to 3 textfiles in order to be used with other programs
        """
        pX = path+'_X'
        pY = path+'_Y'
        pZ = path+'_Z'
        np.savetxt(pX, self.mX, delimiter=',')
        np.savetxt(pY, self.mY, delimiter=',')
        np.savetxt(pZ, self.data, delimiter=',')
 
