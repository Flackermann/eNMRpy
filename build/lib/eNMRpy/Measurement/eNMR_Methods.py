from .base import Measurement
from sklearn.linear_model import huber as hub
import numpy as np
import nmrglue as ng
import matplotlib.pyplot as plt

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
    # def get_xaxis(self):
    #
    #     if self.dependency is not None:
    #         return {
    #             "G": self.eNMRraw["g in T/m"],
    #             "U": self.eNMRraw["U / [V]"]
    #         }[self.dependency.upper()]
    #     else:
    #         print("an error occured, the dependency is None")
    def __init__(self, path, measurement, alias=None,  linebroadening=5):
        super().__init__(path, measurement, alias=None,  linebroadening=5)
        
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
        self.delta*1000, self.g, self.d*1000
        )
    
    def __getitem__(self, key):
        if type(key) == str:
            return self.eNMRraw[key]
        elif type(key) == int:
            return (self.ppm, self.data[key])
        
    #def __add__(self, other):
        
    
    def spam_and_eggs(self):
        """
        For Class inheritance test purposes
        """
        #self.x_var = {None:None,
                #"U": "U / [V]",
                #"G": "g in T/m"
        #}
        print(self.x_var)
        print("SPAM and EGGS!")
        

    def autophase(self,
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
                self.eNMRraw.loc[n, "ph0acme"] = _val[0]  # *-1
                if progress:
                    clear_output() # keeps the output clean. remove for more info on iteration steps etc.
            elif method == 'difference':
                self.eNMRraw.loc[n, "ph0diff"] = _val[0]  # *-1
                if progress:
                    clear_output() # keeps the output clean. remove for more info on iteration steps etc.
            else:
                self.eNMRraw.loc[n, "ph0sqdiff"] = _val[0]  # *-1
                if progress:
                    clear_output() # keeps the output clean. remove for more info on iteration steps etc.
       
        #corr_enmr = self.eNMRraw.sort_values("U / [V]")
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

            #self.eNMRraw["ph0acmereduced"] = corr_enmr['ph0acme']/(self.eNMRraw['g in T/m'][0]*self.Delta*self.delta*self.gamma*self.d)
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

            #self.eNMRraw["ph0diffreduced"] = corr_enmr['ph0diff']/(self.eNMRraw['g in T/m'][0]*self.Delta*self.delta*self.gamma*self.d)
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

    def _autops(self, data, _fnc="acme", p0=0.0):  # , p1=0.0):
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
        # self.p1 = p1
        
        if not callable(_fnc):
            self._fnc = {
                # 'peak_minima': _ps_peak_minima_score,
                'acme': self._ps_acme_score,
                'difference': self._ps_abs_difference_score,
                'sqdifference': self._ps_sq_difference_score
            }[self._fnc]

        # self.opt ist die optimierte Phase für das Spektrum.
        self.opt = self._p0  # [self._p0, self.p1]
        self.opt = fmin(self._fnc, x0=self.opt, args=(data, ), disp=0)

        # self.phasedspc = ng.proc_base.ps(data, p0=self.opt[0], p1=self.opt[1])

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
        
        # s0 --> initial spectrum?
        s0 = ng.proc_base.ps(data, p0=-ph, p1=0)  # , p1=phc1 --> p1=0 lets the algorithm optimize only p0
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

        return: fig
        """
        _data = {'orig': self.data_orig,
                'cropped': self.data}[data]

        # x_axis = {}[self.dependency]

        # the correction factor for the normalization is added again.
        self.phased = [ng.proc_base.ps(_data[n, :], p0=-(self.eNMRraw.loc[n, ph_var]))# + self.U_0))
                  for n in range(len(_data[:, 0]))]

        intensity = np.array([self.phased[n].real.sum() for n in range(len(_data[:, 0]))])
        
        if normalize:
            intensity /= intensity[0]
                
        #u = [self.eNMRraw.loc[i, 'U / [V]'] for i, n in enumerate(self.eNMRraw['U / [V]'])]
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
        _cx.set_xlabel('$U$ / V')
        _cx.set_ylabel('$\t{\Delta}\phi$ / °')

        _dx = plt.subplot(224, sharex=_bx)
        _dx.plot(intensity_data['ph'], 'ok')
        _dx.set_xlabel('vc')

        fig.savefig(self.path+'intensity_plot_'+self.expno+".pdf")

        intensity_data.to_csv(self.path+'intensity_data_'+self.expno+".csv")

        return fig#, intensity_data

    def lin_huber(self, epsilon=1.35, ulim=None, y_column='ph0acme'):
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
        self._eNMRreg = self.eNMRraw[[self._x_axis, y_column]].sort_values(self._x_axis)
        
        # setting the axis for regression
        if ulim is None:
            self.umin = min(self.eNMRraw[self._x_axis])
        else:
            self.umin = ulim[0]
            
        if ulim is None:
            self.umax = max(self.eNMRraw[self._x_axis])
        else:
            self.umax = ulim[1]

        self._npMatrix = np.matrix(self._eNMRreg[(self.eNMRraw[self._x_axis] <= self.umax)
                                                 == (self.eNMRraw[self._x_axis] >= self.umin)])

        self._X_train, self._Y_train = self._npMatrix[:, 0], self._npMatrix[:, 1]
        
        # regression object
        self.huber = hub.HuberRegressor(epsilon=epsilon)
        self.huber.fit(self._X_train, self._Y_train)
        
        # linear parameters
        self.m = self.huber.coef_  # slope
        self.b = self.huber.intercept_  # y(0)
        self._y_pred = self.huber.predict(self._X_train)
        self._y_pred = self._y_pred.reshape(np.size(self._X_train), 1)
        
        # drop the outliers
        self._outX_train = np.array(self._X_train[[n == False for n in self.huber.outliers_]])
        self._outY_train = np.array(self._Y_train[[n == False for n in self.huber.outliers_]])
        self._outY_pred = np.array(self._y_pred[[n == False for n in self.huber.outliers_]])
        
        # mark outliers in dataset
        # self._inliers = [n is not True for n in self.huber.outliers_]

        self.eNMRraw["outlier"] = True

        for n in range(len(self._npMatrix[:, 0])):
            self.eNMRraw.loc[self.eNMRraw[self._x_axis] == self._npMatrix[n, 0], "outlier"] = self.huber.outliers_[n]

        # calculation of the slope deviation
        _sig_m_a = np.sqrt(np.sum((self._outY_train-self._outY_pred)**2)/(np.size(self._outY_train)-2))
        _sig_m_b = np.sqrt(np.sum((self._outX_train-self._outX_train.mean())**2))
        self.sig_m = _sig_m_a/_sig_m_b

        # debug
        #print(self.sig_m)

        # R^2
        self.r_square = self.huber.score(self._outX_train, self._outY_train)
        
        self.lin_res_dic[y_column] = {'b': self.b,
                                   'm': self.m,
                                   'r^2': self.r_square,
                                   'x': np.array(self._X_train.tolist()).ravel(),
                                   'y': self._y_pred.ravel(),
                                   'sig_m': self.sig_m}

    def lin_display(self, ylim=None, show_slope_deviation=True, n_sigma_displayed=1, dpi=500, y_column='ph0acme', textpos=(0.5,0.15), extra_note=''):
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

        print("formula: y = {0}x + {1}".format(self.m,self.b))
        
        # create figure
        fig_enmr = plt.figure()

        # sublot phase data
        _ax = fig_enmr.add_subplot(111)
        
        # color format for outliers
        colors = ["r" if n else "k" for n in self.eNMRraw.sort_values(self._x_axis)['outlier']]
        
        _ax.scatter(x=np.ravel(self._eNMRreg[self._x_axis]),
                    y=np.ravel(self._eNMRreg[y_column]),
                    marker="o",
                    c=colors)
        _ax.set_ylim(ylim)

        # format the data for plotting
        _xdata = np.ravel(self._X_train)
        _ydata = np.ravel(self._y_pred)
        
        # Plot the regression
        _ax.plot(_xdata, _ydata, "r-")

        if show_slope_deviation:
            _ax.fill_between(_xdata, _xdata*(self.m+n_sigma_displayed*self.sig_m)+self.b,
                                     _xdata*(self.m-n_sigma_displayed*self.sig_m)+self.b,
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
                 "y = %.4f $\cdot$ x + %4.2f\n$R^2$=%4.3f; $\sigma_m=$%4.4f"%(self.m,
                                                                              self.b,
                                                                              self.r_square,
                                                                              self.sig_m),
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
                    ax.plot(self.lin_res_dic[col]['x'], self.lin_res_dic[col]['y']-corr, '--', label='%s lin regression'%col, c=colors[i], marker=None)
        #if type(cols) == list:
            #for i, col in enumerate(cols):
                #try:
                    #ax.errorbar(self._x_axis, col, yerr='%s_err'%col, fmt='o', data=self.eNMRraw, label=message[col], c=colors[i], marker=markers[i])
                #except KeyError:
                    #ax.errorbar(self._x_axis, col, yerr='%s_err'%col, fmt='o', data=self.eNMRraw, label='fitted data %s'%col, c=colors[i], marker=markers[i])
                #except ValueError:
                    #ax.errorbar(self._x_axis, col, fmt='o', data=self.eNMRraw, label=message[col], c=colors[i], marker=markers[i])
                #if regression:
                    #ax.plot(self.lin_res_dic[col]['x'], self.lin_res_dic[col]['y'], '--', label='%s lin regression'%col, c=colors[i], marker=None)
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
                    ax.plot(self.lin_res_dic[cols]['x'], self.lin_res_dic[cols]['y']-corr, '--', label='%s lin regression'%cols, c=colors[0], marker=None)
        #else:
            #try:
                #ax.errorbar(self._x_axis, cols, yerr='%s_err'%cols, fmt='o', data=self.eNMRraw, label='Phasefitting of Peak %s' %cols, c=colors[0], marker=markers)
                #if regression:
                    #ax.plot(self.lin_res_dic[cols]['x'], self.lin_res_dic[cols]['y'], '--', label='%s lin regression'%cols, c=colors[0], marker=markers)
            #except KeyError:
                    #ax.errorbar(self._x_axis, col, yerr='%s_err'%col, fmt='o', data=self.eNMRraw, label='fitted data %s'%cols, c=colors[i], marker=markers)
            #except ValueError:
                #ax.errorbar(self._x_axis, cols, fmt='o', data=self.eNMRraw, label=message[cols], c=colors[0], marker=markers)
                #if regression:
                    #ax.plot(self.lin_res_dic[cols]['x'], self.lin_res_dic[cols]['y'], '--', label='%s lin regression'%cols, c=colors[0], marker=None)
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

    def plot_spec_phcorr(self, phasekey='ph0acme', xlim=None, save=False, savepath=None, show=True, ppm=True, orig_data=True, x_legend=1.1, y_legend= 1.0, ncol_legend=2):
        """
        plots phase corrected rows

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
            
        #xmin = 0 if xmin is None else xmin
        fig, ax = plt.subplots()
        
        if not ppm:
            #xmax = self.TD if xmax is None else xmax
            #xmax = len(self.data[0,:]) if xmax is None else xmax
            for n in range(len(self.data[:, 0])):
                ax.plot(phasedspc[n].real)
            ax.set_xlim(xlim)
            ax.set_xlabel("data points")
            ax.set_ylabel("intensity (a.u.)")

        else:
            #xmax = self.ppm[0] if xmax is None else xmax
            #xmin = self.ppm[-1] if xmin is None else xmin
            #ixmax, ixmin = np.where(self.ppm >= xmin)[0][0], np.where(self.ppm >= xmax)[0][1]
#             irange = np.where(self._ppmscale <= xmin)[0][0], np.where(self._ppmscale <= xmax)[0][1]

            if orig_data:
                for n in range(len(self.data_orig[:, 0])):
                    # plt.plot(self._ppmscale[ixmin:ixmax], phasedspc[n].real)
                    ax.plot(self.ppm, phasedspc[n].real, label="row %i" %n)
            else:
                for n in range(len(self.data[:, 0])):
                    # plt.plot(self._ppmscale[ixmin:ixmax], phasedspc[n].real)
                    ax.plot(self.ppm, phasedspc[n].real, label="row %i" %n)
            # plt.axis(xmin=self._ppm_l-self.dic["acqus"]["SW"]*xmin/100,
            #          xmax=self._ppm_r+self.dic["acqus"]["SW"]*(1-xmax/100))
            ax.set_xlim(xlim)
            ax.set_xlabel("$\delta$ / ppm")
            ax.set_ylabel("intensity / a.u.")
            ax.legend(bbox_to_anchor=(x_legend, y_legend), ncol=ncol_legend)

        #ax = plt.gca()
        #ax.set_xlim(ax.get_xlim()[::-1])  # inverts the x-axis
        ax.ticklabel_format(style='sci', scilimits=(0,0), axis='y')

        if save:
            fig.savefig(path, dpi=500)
        if show:
            plt.show()
        #print("ixmax:", ixmax)
        #print("ixmin:", ixmin)
#         return ixmax, ixmin
        phasedspc = np.asarray(phasedspc)
        
        return fig, ax
        #if orig_data:
            #return phasedspc[:,ixmin:ixmax]
        #else:
            #return phasedspc
        # ps = ng.proc_base.ps

    def plot_spec_comparison_to_0(self, row, xmax=None, xmin=None, ppm=True):
        """
        plots row 0 and row n in the range of xmax and xmin
        """
        _max = None if xmax is None else xmax
        _min = None if xmin is None else xmin

        fig = plt.figure()
        if ppm:
            _max = self.ppm[0] if xmax is None else xmax
            _min = self.ppm[-1] if xmin is None else xmin
            
            plt.plot(self.ppm, self.data[0, ::1].real, label='row ' + str(0))
            plt.plot(self.ppm, self.data[row, ::1].real, label='row ' + str(row))
            plt.legend()
            plt.axis(xmax=_max, xmin=_min)
            plt.title("this is row %i at %.0f V" % (row, self.eNMRraw.loc[row, self._x_axis]))
            plt.xlabel("ppm")
            plt.ylabel("intensity / a.u.")
        if not ppm:
            plt.plot(self.data[0, ::1].real, label='row '+str(0))
            plt.plot(self.data[row, ::1].real, label='row '+str(row))
            plt.legend()
            plt.axis(xmax=_max, xmin=_min)
            plt.title("this is row %i at %.0f V" % (row, self.eNMRraw.loc[row, self._x_axis]))
            plt.xlabel("datapoints")#.sum()
            plt.ylabel("intensity / a.u.")

        plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

        return fig

    def output_data(self):
        """
        saves the raw phase data in the measurement-folder as a csv
        """
        self.eNMRraw.to_csv(self.path+"phase_data_"+self.expno+".csv", sep=" ")


    def output_results(self, path=None):
        """
        saves the mobility result data in the measurement-folder similar to obj.output_data()
        """
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
            self.mu = (m*d)/(self.gamma*self.delta*self.Delta*g)
            self.lin_res_dic[y_column]['mu']=self.mu
            self.lin_res_dic[y_column]['mu_err']=self.mu*(sig_m/m)
            #return self.mu, self.mu*(sig_m/m)
            self.sig_m, self.m = sig_m, m
        
        if verbose:
            print ('%.2E (m^2/Vs)'%self.mu[0],'+- %.2E'%(self.mu*(self.sig_m/self.m)))
        
        return self.mu, self.mu*(self.sig_m/self.m)


    def analyzePhasecorrection(self, linebroadening=10, lin_threshhold=2.5,
                           graph_save=False, savepath=None, method="acme", xmin=None,
                           xmax=None, cropmode="absolute", progress=True, umin=None, umax=None, output_path=None):

        """
        standard phase correction analyzing routine for eNMR
        
        linebroadening:
            sets the linebroadening of the fourier transformation
        lin_threshhold:
            sets the threshhold for outlier detection of the linear regression
        graph_save:
            True: saves the graph as a png in the measurment folder
        savepath:
            None: Takes the standard-values from the measurement
            'string': full path to customized graph
                    works with .png, .pdf and .eps suffixes
        method: chooses between phase correction algorithms
			"acme": standard method using entropy minimization
			"difference": --> _ps_abs_difference_score()
                minimizes the linear difference to the first spectrum by fitting p0
            "sqdifference": --> _ps_sq_difference_score()
                minimizes the square difference to the first spectrum by fitting p0
        xmin, xmax:
            min and maximum x-values to crop the data for processing (and display)
            can take relative or absolute values depending on the cropmode.
        
        cropmode: changes the x-scale unit
            "percent": value from 0% to 100% of the respective x-axis-length --> does not fail
            "absolute": takes the absolute values --> may fail
        progress: This shows the spectral row that is processed in this moment. May be disabled in order to be able to stop clearing the output.
        """
        if output_path is None:
            output_path = 'expno_%i_results.xlsx'%self.expno
        
        self.proc(linebroadening=linebroadening, xmin=xmin, xmax=xmax, cropmode=cropmode)
        self.autophase(analyze_only=False, method=method, progress=progress)
        self.lin_huber(epsilon=lin_threshhold, umin=umin, umax=umax)
        # obj.lin() #---> dies ist die standardn, least squares method.
        self.lin_display(save=graph_save, dpi=330, savepath=savepath)
        self.mobility()
        self.output_all_results(output_path)
 
