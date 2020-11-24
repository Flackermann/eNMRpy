from .base import Measurement
from sklearn.linear_model import huber as hub
import numpy as np
import nmrglue as ng
import matplotlib.pyplot as plt
import pandas as pd
import lmfit as lf

class _eNMR_Methods(Measurement):
    """
    This is the subclass of Masurement() containing all methods
    
    path:
        relative or absolute path to the measurements folder
    measurement:
        the to the experiment corresponding EXPNO
    alias:
        Here you can place an individual name relevant for plotting. If None, the path is taken instead.
    """

    def __init__(self, path, expno, alias=None,  lineb=5):
        super().__init__(path, expno, alias=None,  lineb=5)
        
        self._x_axis = {"U": "U / [V]",
                   "G": "g in T/m",
                   "I": "I / mA",
                   "RI": "RI / V"
                   }[self.dependency.upper()]
    
    def __repr__(self):
        return '''%s, expno %s, Delta = %.1fms, ppm range: %.1f to %.1f
    delta= %.1fms, g= %.3f T/m, e-distance=%.0fmm'''%(
        self.nuc, self.expno, self.Delta*1000, 
        self.ppm[0], self.ppm[-1],
        self.delta*1000, self.g, self.d*1000 # outputs the gradient time and diffusion time in ms
        )
    
    def __getitem__(self, key):
        if type(key) == str:
            return self.eNMRraw[key]
        elif type(key) == int:
            return (self.ppm, self.data[key])                

    def autophase_phase_analysis(self,
                  returns=False,
                  method="acme",
                  progress=False,
                  period_compensation=True,
                  normalize=True):
        """
        analyzes the phase of the spectral data and returns phased data

        at this point no data cropping --> full spectrum is analyzed

        returns: if True, returns the raw phased Data. If false, returns nothing

        method: chooses the method for the phase correction
            "acme": standard method using entropy minimization
            "difference": --> _ps_abs_difference_score()
                minimizes the linear difference to the first spectrum by fitting p0
            "sqdifference": --> _ps_sq_difference_score()
                minimizes the square difference to the first spectrum by fitting p0

        progress:
            prints the progress. In case you are dealing with large datasets

        period_compensation:
            corrects for the 360 degree error that occasionaly occurs since a phase of 0 and 360
            is equivalent and indistinguishable

        normalize:
            normalizes all phase data to 0 degrees for 0V. This compensates for different reference 
            phases during the aquisition and phase analysis.
        """
        if self.dependency.upper() == 'G':
            normalize = False

        data_ph = np.array([])
        for n in range(len(self.data[:, 0])):
            if progress:
                print("row %i / %i" % (n+1, len(self.data[:, 0])))

            #####################################
            # using the data with the algorithm
            _val = self._autops(self.data[n, :], _fnc=method)
            #####################################

            data_ph = np.append(data_ph, _val)
                
            if method == 'acme':
                self.eNMRraw.loc[n, "ph0acme"] = _val[0]  
                if progress:
                    clear_output() # keeps the output clean. remove for more info on iteration steps etc.
            elif method == 'difference':
                self.eNMRraw.loc[n, "ph0diff"] = _val[0]  
                if progress:
                    clear_output() # keeps the output clean. remove for more info on iteration steps etc.
            else:
                self.eNMRraw.loc[n, "ph0sqdiff"] = _val[0]  
                if progress:
                    clear_output() # keeps the output clean. remove for more info on iteration steps etc.
       
        corr_enmr = self.eNMRraw.sort_values(self._x_axis)
        
        if method == 'acme':
            if period_compensation:
                for m in range(corr_enmr["ph0acme"].size):
                    if corr_enmr["ph0acme"].iloc[m]-corr_enmr["ph0acme"].iloc[abs(m-1)] < -300:
                        corr_enmr["ph0acme"].iloc[m] += 360
                    elif corr_enmr["ph0acme"].iloc[m]-corr_enmr["ph0acme"].iloc[abs(m-1)] > 300:
                        corr_enmr["ph0acme"].iloc[m] -= 360
                self.eNMRraw = corr_enmr
        
            if normalize:
                #self.U_0 = corr_enmr[corr_enmr["U / [V]"] == 0]["ph0acme"][0]
                self.U_0 = corr_enmr[corr_enmr[self._x_axis] == 0]["ph0acme"].iloc[0]
                for m in range(corr_enmr["ph0acme"].size):
                    corr_enmr["ph0acme"].iloc[m] -= self.U_0
                self.eNMRraw = corr_enmr

            self.eNMRraw["ph0acmereduced"] = corr_enmr['ph0acme']*self.d/(self.eNMRraw['g in T/m'][0]*self.Delta*self.delta*self.gamma)
        elif method == 'difference':
            if period_compensation:
                for m in range(corr_enmr["ph0diff"].size):
                    if corr_enmr["ph0diff"].iloc[m]-corr_enmr["ph0diff"].iloc[abs(m-1)] < -300:
                        corr_enmr["ph0diff"].iloc[m] += 360
                    elif corr_enmr["ph0diff"].iloc[m]-corr_enmr["ph0diff"].iloc[abs(m-1)] > 300:
                        corr_enmr["ph0diff"].iloc[m] -= 360
                self.eNMRraw = corr_enmr
        
            if normalize:
                self.U_0 = corr_enmr[corr_enmr[self._x_axis] == 0]["ph0diff"][0]
                for m in range(corr_enmr["ph0diff"].size):
                    corr_enmr["ph0diff"].iloc[m] -= self.U_0
                self.eNMRraw = corr_enmr

            self.eNMRraw["ph0diffreduced"] = corr_enmr['ph0diff']*self.d/(self.eNMRraw['g in T/m'][0]*self.Delta*self.delta*self.gamma)
        else:
            if period_compensation:
                for m in range(corr_enmr["ph0sqdiff"].size):
                    if corr_enmr["ph0sqdiff"].iloc[m]-corr_enmr["ph0sqdiff"].iloc[abs(m-1)] < -300:
                        corr_enmr["ph0sqdiff"].iloc[m] += 360
                    elif corr_enmr["ph0sqdiff"].iloc[m]-corr_enmr["ph0sqdiff"].iloc[abs(m-1)] > 300:
                        corr_enmr["ph0sqdiff"].iloc[m] -= 360
                self.eNMRraw = corr_enmr
        
            if normalize:
                self.U_0 = corr_enmr[corr_enmr[self._x_axis] == 0]["ph0sqdiff"][0]
                for m in range(corr_enmr["ph0sqdiff"].size):
                    corr_enmr["ph0sqdiff"].iloc[m] -= self.U_0
                self.eNMRraw = corr_enmr

            self.eNMRraw["ph0sqdiffreduced"] = corr_enmr['ph0sqdiff']*self.d/(self.eNMRraw['g in T/m'][0]*self.Delta*self.delta/1000*self.gamma)

        print("done, all went well")
        
        if returns is True:
            return data_ph  # , data_spec

    def _autops(self, data, _fnc="acme", p0=0.0, p1=None):  # , p1=0.0):
        """
        FS: modified for eNMR purpose. optimizes only p0
        ----------
        Automatic linear phase correction

        Parameters
        ----------
        data : ndarray
            Array of NMR data.
        _fnc : str or function
            Algorithm to use for phase scoring. Built in functions can be
            specified by one of the following strings: "acme", "peak_minima"
        p0 : float
            Initial zero order phase in degrees.
        p1 : float
            Initial first order phase in degrees.

        Returns
        -------
        ndata : ndarray
            Phased NMR data.

        """
        
        from  scipy.optimize import fmin
        
        self._fnc = _fnc
        self._p0 = p0
        self._p1 = p1
        
        if not callable(_fnc):
            self._fnc = {
                # 'peak_minima': _ps_peak_minima_score,
                'acme': self._ps_acme_score,
                'difference': self._ps_abs_difference_score,
                'sqdifference': self._ps_sq_difference_score
            }[self._fnc]
        
        # this is the behavior for the analysis in the eNMR context
        if p1 == None:
            # self.opt ist die optimierte Phase für das Spektrum.
            self.opt = self._p0  # [self._p0, self.p1]
            self.opt = fmin(self._fnc, x0=self.opt, args=(data, ), disp=0)

            # self.phasedspc = ng.proc_base.ps(data, p0=self.opt[0], p1=self.opt[1])

            return self.opt  # self.phasedspc, self.opt
        
        # this is used to autophase the spectrum and get the values zero and first order
        elif p1 != None:
            # self.opt ist die optimierte Phase für das Spektrum.
            self.opt = (self._p0, self._p1)
            self.opt = fmin(self._fnc, x0=self.opt, args=(data, ), disp=0)

            #self.phasedspc = ng.proc_base.ps(data, p0=self.opt[0], p1=self.opt[1])
            return self.opt  # self.phasedspc, self.opt

    @staticmethod
    def _ps_acme_score(ph, data):
        """
        FS: modified for eNMR purpose. optimizes only p0
        ------------
        Phase correction using ACME algorithm by Chen Li et al.
        Journal of Magnetic Resonance 158 (2002) 164-168

        using only the first derivative for the entropy

        Parameters
        ----------
        pd : tuple
            Current p0 and p1 values
        data : ndarray
            Array of NMR data.

        Returns
        -------
        score : float
            Value of the objective function (phase score)

        """
        stepsize = 1
        

        if (len(ph) == 1):
            s0 = ng.proc_base.ps(data, p0=-ph, p1=0)  # , p1=phc1 --> p1=0 lets the algorithm optimize only p0
        
        else:
            s0 = ng.proc_base.ps(data, p0=-ph[0], p1=ph[1])  # optimizing for ph0 and ph1
        
        data = np.real(s0)
        

        # Calculation of first derivatives --> hier wird das absolute Spektrum erzeugt
        ds1 = np.abs((data[1:]-data[:-1]) / (stepsize*2))
        p1 = ds1 / np.sum(ds1)

        # Calculation of entropy
        p1[p1 == 0] = 1  # was macht diese Zeile? das Verstehe ich noch nicht richtig

        h1 = -p1 * np.log(p1)
        h1s = np.sum(h1)

        # Calculation of penalty
        pfun = 0.0
        as_ = data - np.abs(data)
        sumas = np.sum(as_)

        if sumas < 0:
            pfun = pfun + np.sum((as_/2) ** 2)

        p = 1000 * pfun

        return h1s + p
    
    def _ps_abs_difference_score(self, ph, data):
        """
        FS: modified for eNMR purpose. optimizes only p0
        ------------
        Parameters
        ----------
        pd : tuple
            Current p0 and p1 values
        data : ndarray
            Array of NMR data.

        Returns
        -------
        score : float
            Value of the objective function (phase score)
        """
       
        s0 = ng.proc_base.ps(data, p0=-ph, p1=0)  # , p1=phc1 --> p1=0 lets the algorithm optimize only p0
        phasedspec = np.real(s0)

        penalty = np.sum(np.abs(self.data[0, :].real - phasedspec))
        
        return penalty
    
    def _ps_sq_difference_score(self, ph, data):
        """
        FS: modified for eNMR purpose. optimizes only p0
        ------------
        Parameters
        ----------
        pd : tuple
            Current p0 and p1 values
        data : ndarray
            Array of NMR data.

        Returns
        -------
        score : float
            Value of the objective function (phase score)
        """
       
        s0 = ng.proc_base.ps(data, p0=-ph, p1=0)  # , p1=phc1 --> p1=0 lets the algorithm optimize only p0
        phasedspec = np.real(s0)

        penalty = np.sum(np.square(self.data[0,:].real - phasedspec))  # *1/n-1 wäre die Varianz,
        # ist hier egal, da alle Spektren die gleiche Anzahl an Punkten haben.
        
        return penalty
    
    def analyze_intensity(self, data='cropped', ph_var='ph0acme', normalize=True, ylim=None):
        """
        uses the phase data information to rephase and integrate all spectra and plot a comparison
        stores the intensity information in the measurement folder
        
        data:
            'orig': self.data_orig
            'cropped': self.data

        returns: fig, intensity_data
        """
        _data = {'orig': self.data_orig,
                'cropped': self.data}[data]

        # the correction factor for the normalization is added again.
        self.phased = [ng.proc_base.ps(_data[n, :], p0=-(self.eNMRraw.loc[n, ph_var]))# + self.U_0))
                  for n in range(len(_data[:, 0]))]

        intensity = np.array([self.phased[n].real.sum() for n in range(len(_data[:, 0]))])
        
        if normalize:
            intensity /= intensity[0]

        u = [self.eNMRraw.loc[i, self._x_axis] for i, n in enumerate(self.eNMRraw[self._x_axis])]

        intensity_data = pd.DataFrame()
        intensity_data["U"] = u
        intensity_data['intensity'] = intensity
        intensity_data['ph'] = self.eNMRraw[ph_var]# + self.U_0

        self.intensity_data = intensity_data

        fig = plt.figure(figsize=(8, 6))
        _ax = plt.subplot(221)
        _ax.scatter(intensity_data['U'], intensity_data['intensity'], c='k')
        _ax.set_ylabel('intensity / a.u.')
        if normalize and (ylim is None):
            _ax.set_ylim(0,1.05)
        elif normalize:
            _ax.set_ylim(*ylim)
        _bx = plt.subplot(222, sharey=_ax)
        _bx.plot(intensity_data['intensity'], 'ok')

        _cx = plt.subplot(223, sharex=_ax)
        _cx.scatter(intensity_data['U'], intensity_data['ph'], c='k')
        
        if self.dependency.upper() == "U":
            _cx.set_xlabel('$U$ / V')
        elif self.dependency.upper() == "G":
            _cx.set_xlabel("$g$ / $($T$\cdot$m$^{-1})$")
        elif self.dependency.upper() == 'I':
            _cx.set_xlabel("$I$ / mA")
        elif self.dependency.upper() == 'RI':
            _cx.set_xlabel("$(R \cdot I)$ / V")

        _cx.set_ylabel('$\t{\Delta}\phi$ / °')

        _dx = plt.subplot(224, sharex=_bx)
        _dx.plot(intensity_data['ph'], 'ok')
        _dx.set_xlabel('vc')

        fig.savefig(self.path+'intensity_plot_'+self.expno+".pdf")

        intensity_data.to_csv(self.path+'intensity_data_'+self.expno+".csv")

        return fig, intensity_data

    def linreg(self, ulim=None, y_column='ph0'):
        """
        standard linear regression method based on the least-square method

        ulim: 
            tuple defining the voltage limits for the regression e.g. ulim = (-100, 100)
        y_column:
            column(keyword) to be analyzed from the eNMRraw dataset

        stores results in lin_res_dic[y_column]
        :returns: nothing
        """

        # select x-axis
        #self._x_axis = {"U": "U / [V]",
                #"G": "g in T/m"
                #}[self.dependency.upper()]

        # convert data
        _eNMRreg = self.eNMRraw[[self._x_axis, y_column]].sort_values(self._x_axis)

        # setting the axis for regression
        if ulim is None:
            umin = min(self.eNMRraw[self._x_axis])
        else:
            umin = ulim[0]

        if ulim is None:
            umax = max(self.eNMRraw[self._x_axis])
        else:
            umax = ulim[1]

        _nparray = np.array(_eNMRreg[(self.eNMRraw[self._x_axis] <= umax)
                                                == (self.eNMRraw[self._x_axis] >= umin)])
        _X_train, _Y_train = _nparray[:, 0], _nparray[:, 1]


        def lineareq(x, m, b):
            return x*m+b
        

        # regression object
        linmodel = lf.Model(lineareq)
        linparams = linmodel.make_params()
        linparams['b'].set(0)
        linparams['m'].set(1)
        result = linmodel.fit(_Y_train, x=_X_train, params=linparams)

        # linear parameters
        m = result.best_values['m']  # slope
        b = result.best_values['b'] # y(0)
        _Y_pred = result.best_fit
    
        # calculation of the slope deviation
        _sig_m_a = np.sqrt(np.sum((_Y_train-_Y_pred)**2)/(np.size(_Y_train)-2))
        _sig_m_b = np.sqrt(np.sum((_X_train-_X_train.mean())**2))
        sig_m = _sig_m_a/_sig_m_b

        # debug
        #print(self.sig_m)

        # R^2    
        r_square = 1 - result.residual.var() / np.var(_Y_train)
        
        self.lin_res_dic[y_column] = {'b': b,
                                'm': m,
                                'sig_m': sig_m,
                                'r_square': r_square,
                                'x': np.array(_X_train.tolist()).ravel(),
                                'y': np.array(_Y_train.tolist()).ravel(),
                                'y_fitted': _Y_pred.ravel(),
                                }
        return

    # should be replaced by a more recent function since it will be deprecated
    def lin_huber(self, epsilon=3, ulim=None, y_column='ph0'):
        """
        robust linear regression method from scikit-learn module based on the least-square method with an additional threshhold (epsilon) for outlying datapoints
        outlying datapoints are marked as red datapoints
        
        epsilon:
            threshhold > 1
        ulim: 
            tuple defining the voltage limits for the regression e.g. ulim = (-100, 100)
        y_column:
            column(keyword) to be analyzed from the eNMRraw dataset
        
        stores results in lin_res_dic[y_column]
        :returns: nothing
        """

        # select x-axis
        #self._x_axis = {"U": "U / [V]",
                   #"G": "g in T/m"
                   #}[self.dependency.upper()]

        # convert data
        _eNMRreg = self.eNMRraw[[self._x_axis, y_column]].sort_values(self._x_axis)
        
        # setting the axis for regression
        if ulim is None:
            umin = min(self.eNMRraw[self._x_axis])
        else:
            umin = ulim[0]
            
        if ulim is None:
            umax = max(self.eNMRraw[self._x_axis])
        else:
            umax = ulim[1]

        _npMatrix = np.matrix(_eNMRreg[(self.eNMRraw[self._x_axis] <= umax)
                                                 == (self.eNMRraw[self._x_axis] >= umin)])

        _X_train, _Y_train = _npMatrix[:, 0], _npMatrix[:, 1]
        
        # regression object
        huber = hub.HuberRegressor(epsilon=epsilon)
        huber.fit(_X_train, _Y_train)
        
        # linear parameters
        m = huber.coef_  # slope
        b = huber.intercept_  # y(0)
        _y_pred = huber.predict(_X_train)
        _y_pred = _y_pred.reshape(np.size(_X_train), 1)
        
        # drop the outliers
        _outX_train = np.array(_X_train[[n == False for n in huber.outliers_]])
        _outY_train = np.array(_Y_train[[n == False for n in huber.outliers_]])
        _outY_pred = np.array(_y_pred[[n == False for n in huber.outliers_]])
        
        # mark outliers in dataset
        # self._inliers = [n is not True for n in self.huber.outliers_]

        self.eNMRraw["outlier"] = True

        for n in range(len(_npMatrix[:, 0])):
            self.eNMRraw.loc[self.eNMRraw[self._x_axis] == _npMatrix[n, 0], "outlier"] = huber.outliers_[n]

        # calculation of the slope deviation
        _sig_m_a = np.sqrt(np.sum((_outY_train-_outY_pred)**2)/(np.size(_outY_train)-2))
        _sig_m_b = np.sqrt(np.sum((_outX_train-_outX_train.mean())**2))
        sig_m = _sig_m_a/_sig_m_b

        # debug
        #print(self.sig_m)

        # R^2
        r_square = huber.score(_outX_train, _outY_train)
        
        self.lin_res_dic[y_column] = {'b': b,
                                   'm': m,
                                   'r_square': r_square,
                                   'x': np.array(_X_train.tolist()).ravel(),
                                   'y': _Y_train,
                                   'y_fitted': _y_pred.ravel(),
                                   'sig_m': sig_m}

    def lin_display(self, ylim=None, show_slope_deviation=True, n_sigma_displayed=1, dpi=500, y_column='ph0', textpos=(0.5,0.15), extra_note=''):
        """
        displays the linear huber regression
        If there is an alias available from measurement object, it will replace the path in the title
        
        ylim:
            set limits of the y-axis
        
        show_slope_deviation:
            display the standard deviation
        n_sigma_displayed:
            multiplicator of the displayed standard deviation
        y_column:
            column(keyword) to be analyzed from the eNMRraw dataset
        textpos:
            tuple for the textposition
        extra_note:
            added to the text
        dpi:
            adjusts the output dpi to the required value
        
        :returns: figure
        """
        
        textx, texty = textpos
        #_x_axis = {"U":"U / [V]", "G":"g in T/m"}
        #self._x_axis = {"U":"U / [V]", "G":"g in T/m"}[self.dependency]

        print("formula: y = {0}x + {1}".format(self.lin_res_dic[y_column]['m']
                                                ,self.lin_res_dic[y_column]['b']))
        
        # create figure
        fig_enmr = plt.figure()

        # sublot phase data
        _ax = fig_enmr.add_subplot(111)
        
        # color format for outliers
        try:
            colors = ["r" if n else "k" for n in self.eNMRraw.sort_values(self._x_axis)['outlier']]
        except KeyError:
            colors = ["k" for n in self.eNMRraw[y_column]]
            print('no outliers registered')
            pass
        
        #_ax.scatter(x=np.ravel(self._eNMRreg[self._x_axis]),
                    #y=np.ravel(self._eNMRreg[y_column]),
        _ax.scatter(x=np.ravel(self.eNMRraw[self._x_axis]),
                    y=np.ravel(self.eNMRraw[y_column]),
                    marker="o",
                    c=colors)
        _ax.set_ylim(ylim)

        # format the data for plotting
        _xdata = np.ravel(self.lin_res_dic[y_column]['x'])
        _ydata = np.ravel(self.lin_res_dic[y_column]['y_fitted'])
        
        # Plot the regression
        _ax.plot(_xdata, _ydata, "r-")

        if show_slope_deviation:
            #_ax.fill_between(_xdata, _xdata*(self.m+n_sigma_displayed*self.sig_m)+self.b,
                                     #_xdata*(self.m-n_sigma_displayed*self.sig_m)+self.b,
                                     #alpha=0.5,
                                     #facecolor="blue")
            _ax.fill_between(_xdata, _xdata*(self.lin_res_dic[y_column]['m']+n_sigma_displayed*self.lin_res_dic[y_column]['sig_m'])+self.lin_res_dic[y_column]['b'],
                                     _xdata*(self.lin_res_dic[y_column]['m']-n_sigma_displayed*self.lin_res_dic[y_column]['sig_m'])+self.lin_res_dic[y_column]['b'],
                                     alpha=0.5,
                                     facecolor="blue")

        # make title
        if self.alias is None:
            title_printed = r'%s'%((self.path+self.expno).split("/")[-2]+", EXPNO: "+self.expno+extra_note)
#             plt.title('LiTFSIDOLDME') # debugging
        else:
            title_printed = self.alias+", EXPNO: "+self.expno+extra_note
#             plt.title(self.alias+", EXPNO: "+self.expno+extra_note)
#             plt.title('test2') # for debugging purposes
        plt.title(title_printed.replace('_',r' '))
    
        if self.dependency.upper() == "U":
            plt.xlabel("$U$ / V")
        elif self.dependency.upper() == "G":
            plt.xlabel("$g$ / $($T$\cdot$m$^{-1})$")
        elif self.dependency.upper() == 'I':
            plt.xlabel("$I$ / mA")
        elif self.dependency.upper() == 'RI':
            plt.xlabel("$(R \cdot I)$ / V")

        plt.ylabel("$\Delta\phi$ / °")
        
        # plotting the Textbox
        plt.text(textx, texty,
                 "y = %.4f $\cdot$ x + %4.2f\n$R^2$=%4.3f; $\sigma_m=$%4.4f"%(self.lin_res_dic[y_column]['m'],
                                                                              self.lin_res_dic[y_column]['b'],
                                                                              self.lin_res_dic[y_column]['r_square'],
                                                                              self.lin_res_dic[y_column]['sig_m']),
                 fontsize=14,
                 bbox={'facecolor':'white', 'alpha':0.7,'pad':10},
                 horizontalalignment='center',
                 verticalalignment='center',
                 transform=_ax.transAxes)

        return fig_enmr

    def lin_results_display(self, cols, regression=True, normalize=True, colors=None, markers=None, fig=None, x_legend=1.1, y_legend=1.0, ncol_legend=2, ymin=None, ymax=None, figsize=None):
        '''
        displays the phase results including the regression. Can take a previous figure from a different eNMR measurement in order to compare results
        the relegend() function can be used on these graphs to reprint the legend() for publishing
        
        cols:
            list of keywords for dataselection from the eNMRraw table
        colors:
            list of colors or a colormap --> see matplotlib documentation
        
        :returns: figure
        '''
        
        #self._x_axis = {"U":"U / [V]", "G":"g in T/m"}[self.dependency.upper()]
        
        if type(cols) != list:
            raise TypeError('cols should be a list')
        
        if fig is None:
            fig = plt.figure(figsize=figsize)
        else:
            fig = fig

        ax = fig.add_subplot(111)

        if colors is None:
            prop_cycle = plt.rcParams['axes.prop_cycle']
            colors = prop_cycle.by_key()['color']
        else:
            colors = colors
        
        if markers is None:
            markers = ['o', '^', 's', '+', '*', 'D', 'v', 'x', '<']
        else:
            markers = markers

        message = {
            'ph0acme': 'Autophase with entropy minimization',
            'ph0diff': 'Autophase with linear difference minimization',
            'ph0sqdiff': 'Autophase with square difference minimization'
            }
        
        if type(cols) == list:
            for i, col in enumerate(cols):
                if normalize:
                    corr = self.eNMRraw.loc[self.eNMRraw['U / [V]']==0, col].iloc[0]
                else:
                    corr = 0
                    
                try:
                    ax.errorbar(self.eNMRraw[self._x_axis], self.eNMRraw[col]-corr, yerr=self.eNMRraw['%s_err'%col], fmt='o', label=message[col], c=colors[i], marker=markers[i])
                except KeyError:
                    if col == 'ph0acme':
                        ax.errorbar(self.eNMRraw[self._x_axis], self.eNMRraw[col]-corr, fmt='o', label=message[col], c=colors[i], marker=markers[i])
                    else:
                        ax.errorbar(self.eNMRraw[self._x_axis], self.eNMRraw[col]-corr, yerr=self.eNMRraw['%s_err'%col], fmt='o', label='fitted data %s'%col, c=colors[i], marker=markers[i])
                except ValueError:
                    ax.errorbar(self.eNMRraw[self._x_axis], self.eNMRraw[col]-corr, fmt='o', label=message[col], c=colors[i], marker=markers[i])
                if regression:
                    ax.plot(self.lin_res_dic[col]['x'], self.lin_res_dic[col]['y_fitted']-corr, '--', label='%s lin regression'%col, c=colors[i], marker=None)

        else:
            if normalize:
                corr = self.eNMRraw.loc[self.eNMRraw['U / [V]']==0, cols][0]
            else:
                corr = 0
                    
            try:
                ax.errorbar(self.eNMRraw[self._x_axis], self.eNMRraw[cols]-corr, yerr=self.eNMRraw['%s_err'%cols], fmt='o', label='Phasefitting of Peak %s' %cols, c=colors[0], marker=markers)
                if regression:
                    ax.plot(self.lin_res_dic[cols]['x'], self.lin_res_dic[cols]['y']-corr, '--', label='%s lin regression'%cols, c=colors[0], marker=markers)
            except KeyError:
                if cols=='ph0acme':
                    ax.errorbar(self.eNMRraw[self._x_axis], self.eNMRraw[cols]-corr, fmt='o', label=message[cols], c=colors[0], marker=markers)
                else:
                    ax.errorbar(self.eNMRraw[self._x_axis], self.eNMRraw[col]-corr, yerr=self.eNMRraw['%s_err'%col], fmt='o', label='fitted data %s'%cols, c=colors[i], marker=markers)
            except ValueError:
                ax.errorbar(self.eNMRraw[self._x_axis], self.eNMRraw[cols]-corr, fmt='o', label=message[cols], c=colors[0], marker=markers)
                if regression:
                    ax.plot(self.lin_res_dic[cols]['x'], self.lin_res_dic[cols]['y_fitted']-corr, '--', label='%s lin regression'%cols, c=colors[0], marker=None)

        ax.legend(bbox_to_anchor=(x_legend, y_legend), ncol=ncol_legend)
        
        xlabel = {'U': r'$U$ / V',
                  'G': r'$g$ / (T\,m$^{-1}$)',
                  'I': '$I$ / mA',
                  'RI': '$(R \cdot I)$ / V'
                      }[self.dependency.upper()]
        
        ax.set_xlabel(xlabel)
        if normalize:
            ax.set_ylabel(r'$\phi-\phi_0$ / °')
        else:
            ax.set_ylabel(r'$\phi$ / °')
        return fig

    def plot_phcorrected_spectra(self, phasekey='ph0acme', xlim=None, save=False, savepath=None, show=True, ppm=True, orig_data=True, x_legend=1.1, y_legend= 1.0, ncol_legend=2):
        """
        plots phase corrected rows stacking the spectra to visualize the analysis result

        ppm: When True, plots the x-axis with ppm scale, xmin and xmax then are percentages
            if False, plots only the number of datapoints, xmin and xmax then are absolute values

        save: saves the spectrum in the data folder

        show: displays the spectrum {plt.show()} or not
        
        returns:
            phased data within given range
        """
        if savepath is None:
            path = self.path+"layered_spectra_"+self.expno+".png"
        else:
            path = savepath
        if orig_data:
            phasedspc = [ng.proc_base.ps(self.data_orig[n, :], p0=-(self.eNMRraw[phasekey][n]))# obsolet: +self.U_0))  # the correction factor for the normalization is added again.
                          for n in range(len(self.data_orig[:, 0]))]  # p1=self.opt[1])
        else:
            phasedspc = [ng.proc_base.ps(self.data[n, :], p0=-(self.eNMRraw[phasekey][n]))# obsolet: +self.U_0))  # the correction factor for the normalization is added again.
                          for n in range(len(self.data[:, 0]))]  # p1=self.opt[1])
            
        fig, ax = plt.subplots()
        
        if not ppm:
            for n in range(len(self.data[:, 0])):
                ax.plot(phasedspc[n].real)
            ax.set_xlim(xlim)
            ax.set_xlabel("data points")
            ax.set_ylabel("intensity (a.u.)")

        else:
            if orig_data:
                for n in range(len(self.data_orig[:, 0])):
                    ax.plot(self.ppm, phasedspc[n].real, label="row %i" %n)
            else:
                for n in range(len(self.data[:, 0])):
                    ax.plot(self.ppm, phasedspc[n].real, label="row %i" %n)
            ax.set_xlim(xlim)
            ax.set_xlabel("$\delta$ / ppm")
            ax.set_ylabel("intensity / a.u.")
            ax.legend(bbox_to_anchor=(x_legend, y_legend), ncol=ncol_legend)

        ax.ticklabel_format(style='sci', scilimits=(0,0), axis='y')

        if save:
            fig.savefig(path, dpi=500)
        if show:
            plt.show()
        phasedspc = np.asarray(phasedspc)
        
        return fig, ax
        
    def output_phase_data(self, path=None):
        """
        saves the raw phase data in the measurement-folder as a csv
        
        if path is None: saves in self.path+"phase_data_"+self.expno+".csv"
        """
        if path != None:
            self.eNMRraw.to_csv(path, sep=" ")
            return
        else:
            self.eNMRraw.to_csv(self.path+"phase_data_"+self.expno+".csv", sep=" ")
            return

    def output_results(self, path=None):
        """
        saves the mobility result data in the measurement-folder similar to obj.output_data()
        """
        from warnings import warn
        warn('this function was renamed output_properties_csv')
        
        results_output = pd.Series([self.nuc,
                                    self.mu[0],
                                    self.sig_m*self.mu[0],
                                    self.d,
                                    self.g,
                                    self.delta,
                                    self.Delta,
                                    self.uInk],
                                   index=["nucleus",
                                          "µ",
                                          "sig_µ",
                                          "d / m",
                                          "g / (T/m)",
                                          "delta / s",
                                          "Delta / s",
                                          "Uink / V"],
                                   name=self.dateipfad)
        if path is None:
            results_output.to_csv(self.path+"mobility_data_"+self.expno+".csv")
        elif path is not None:
            results_output.to_csv(path+"mobility_data_"+self.expno+".csv")
        else:
            print('ooops!')
            
    def output_properties_csv(self, path=None):
        """
        saves the mobility result data in the measurement-folder similar to obj.output_phase_data()
        """
        self.output_results(self, path)
    
    def output_all_results(self, path=None, data=False):
        """
        saves the mobility result data in the measurement-folder similar to obj.output_data()
        """
        
        from pandas import ExcelWriter
        from openpyxl import load_workbook
        
        try:
            book = load_workbook(path)
            writer = pd.ExcelWriter(path, engine='openpyxl') 
            writer.book = book
            writer.sheets = dict((ws.title, ws) for ws in book.worksheets)
        except:
            writer = ExcelWriter(path, engine='xlsxwriter')

        for key in self.lin_res_dic.keys():
            #self.mobility(key)
            results_output = pd.Series([self.nuc,
                                        self.expno,
                                        self.lin_res_dic[key]['mu'][0],
                                        self.lin_res_dic[key]['mu_err'][0],
                                        self.d,
                                        self.g,
                                        self.delta,
                                        self.Delta,
                                        self.uInk],
                                    index=["nucleus",
                                           "expno",
                                            "µ",
                                            "sig_µ",
                                            "d / m",
                                            "g / (T/m)",
                                            "delta / s",
                                            "Delta / s",
                                            "Uink / V"],
                                    name=self.dateipfad)
            results_output.to_excel(writer, sheet_name=self.nuc+'expno'+str(self.expno)+'_'+key)
        
        if data:
            try:
                self.eNMRraw.drop(['data', 'fid'], 1).to_excel(writer, sheet_name=self.nuc+'expno'+str(self.expno)+'_data')
            except KeyError:
                self.eNMRraw.to_excel(writer, sheet_name=self.nuc+'expno'+str(self.expno)+'_data')
            except:
                print('oops, an unexpected error occured')
            writer.save()
        return
    
    def mobility(self, y_column='ph0acme', electrode_distance=None, verbose=True):
        """
        calculates and returns (mobility, deviation) from the regression data
        """
        
        #self.lin_res_dic = {y_column: {'b': self.b,
                            #'m': self.m,
                            #'r^2': self.r_square,
                            #'y_reg': self._y_pred,
                            #'sig_m': self.sig_m}}
        if electrode_distance is None:
            d = self.d
        else:
            d = electrode_distance
        
        if self.dependency.upper() == "G":
            g = self.uInk
        elif self.dependency.upper() == "U":
            g = self.g
        elif self.dependency.upper() == "I":
            g = self.g
        elif self.dependency.upper() == "RI":
            g = self.g
        else:
            print("no dependency was set")
        
        if y_column is None:
            self.mu = (self.m*d)/(self.gamma*self.delta*self.Delta*g)
            #return self.mu, self.mu*(self.sig_m/self.m)

        else:
            m = self.lin_res_dic[y_column]['m']
            sig_m = self.lin_res_dic[y_column]['sig_m']
            mu = (m*d)/(self.gamma*self.delta*self.Delta*g)
            self.lin_res_dic[y_column]['mu']= mu
            self.lin_res_dic[y_column]['mu_err']= mu*(sig_m/m)
            #return self.mu, self.mu*(sig_m/m)
            #self.sig_m, self.m = sig_m, m
        
        if verbose:
            print ('%.2E (m^2/Vs)'%self.lin_res_dic[y_column]['mu'],'+- %.2E'%(self.lin_res_dic[y_column]['mu']*(self.lin_res_dic[y_column]['sig_m']/self.lin_res_dic[y_column]['m'])))
        
        return self.lin_res_dic[y_column]['mu'], self.lin_res_dic[y_column]['mu']*(self.lin_res_dic[y_column]['sig_m']/self.lin_res_dic[y_column]['m'])

