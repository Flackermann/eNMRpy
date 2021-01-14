import matplotlib.pyplot as plt
import nmrglue as ng
from re import findall
import numpy as np
import pandas as pd
from io import StringIO

class Measurement(object):
    """
    This is the base class for the analysis of eNMR raw spectra obtained on Bruker spectrometers with given settings.
    if creating an object from an eNMR measurement creates an error, Uink might have failed to be set and needs to be set manually instead.
    
    path:
        relative or absolute path to the measurements folder
    measurement:
        the to the experiment corresponding EXPNO
    alias:
        Here you can place an individual name relevant for plotting. If None, the path is taken instead.
        
    lineb:
        linebroadening
    """

    def __init__ (self, path, expno, alias=None,  lineb=.3, n_zf_F2=2**14):
        self.path = path
        self.expno = str(expno)
        self.dateipfad = self.path+self.expno
        self.alias = alias
        self.name = self.dateipfad.split(sep='/')[-2]+'/'+self.expno
        self.lineb = lineb
        
        self.comment_string = None
        
        # takes the title page to extract the volt increment
        try:
            title = open(self.dateipfad+"/pdata/1/title").read()
            # reformatting the title_page removing excess of newlines (\n)

        except UnicodeDecodeError:
            title = open(self.dateipfad+"/pdata/1/title", encoding='ISO-8859-15').read()
            
        title = findall('.+', title)              
        self.title_page = self.name+'\n'

        for n in title:
            self.title_page += n+'\n'        

        # read in the bruker formatted data
        self.dic, self.data = ng.bruker.read(self.dateipfad)
        self.pdata = ng.bruker.read_procs_file(self.dateipfad+"/pdata/1")
        # original dataset
        self.data_orig = self.data

        # Berechnung der limits der ppm-skala
        self._ppm_r = self.pdata["procs"]["OFFSET"]
        self._ppm_l = -(self.dic["acqus"]["SW"]-self.pdata["procs"]["OFFSET"])
        self._ppmscale = np.linspace(self._ppm_l, self._ppm_r, n_zf_F2)  # np.size(self.data[0,:]))
        self.ppm = self._ppmscale

        # Bestimmung der dictionaries für die NMR-Daten aus den Messdateien
        # --> wichtig für die Entfernung des digitalen Filters
        #self.udic = ng.bruker.guess_udic(self.dic, self.data)
        #self._uc = ng.fileiobase.uc_from_udic(self.udic)

        # getting some other variables
        self.TD = self.dic["acqus"]["TD"]
        self.fid = self.data.copy()
        self.n_zf_F2 = n_zf_F2

        # the gamma_values in rad/Ts
        gamma_values = {'1H':26.7513e7,
                        '7Li': 10.3962e7,
                        '19F': 25.1662e7}

        self.gamma = gamma_values[self.dic["acqus"]["NUC1"]]
        # conversion from rad in °
        self.gamma = self.gamma/2/np.pi*360
        
        self.nuc = self.dic["acqus"]["NUC1"]
        
        # initialize dictionary for linear regression results
        self.lin_res_dic = {}
        #self.dependency = None

        # END __INIT__
    
    def ph_ref(self, ref):
        '''
        reference your phase shifts on the shifts of one peak.
        ref: string, name of the peak which phase is used for reference (usually solvent).
        '''

        ref_orig_ph = self.eNMRraw[ref].copy
        for i in filter(lambda x: x[:2] =='ph' and x[-1]!='r', self.eNMRraw.keys()):
            self.eNMRraw[i] = self.eNMRraw[i] - ref_orig_ph
    
    
    def comment(self, s, append_to_title=False):
        '''
        write a comment as a multi-line-string or normal string
        '''
        self.comment_string = s
        
        if append_to_title:
            #self._title_page_orig = self.title_page
            self.title_page += self.comment_string

    def calibrate_ppm(self, ppmshift):
        """
        calibrate ppm-scale by adding the desired value.

        :param ppmshift: value added to ppm-scale
        :return: nothing
        """
        self.ppm = self._ppmscale + ppmshift

    def proc(self, linebroadening=None, phc0=0, phc1=0, zfpoints=None, xmin=None, xmax=None, cropmode="percent"):
        """
        processes the spectral data by:
            - removing digital filter
            - create separate fid data
            - linebroadening on spectral data
            - zero filling
            - fourier transformation
            
        crops the data after fft set to xmin and xmax on the x-axis and returns the value
        when xmin and xmax or not both None
        
        xmin, xmax:
            min and maximum x-values to crop the data for processing (and display)
            can take relative or absolute values depending on the cropmode.
        
        cropmode: changes the x-scale unit
            "percent": value from 0% to 100% of the respective x-axis-length --> does not fail
            "absolute": takes the absolute values --> may fail            
        
        :return: nothing
        """
        
        
        if linebroadening is None:
            linebroadening = self.lineb
        
        if zfpoints is not None:
            zfp = zfpoints
        else:
            zfp = self.n_zf_F2
        _lineb = linebroadening
        
        # remove the digital filter
        self.data = ng.bruker.remove_digital_filter(self.dic, self.data)

        # process the spectrum
        # linebroadening
        self.data = ng.proc_base.em(self.data, lb=_lineb/self.dic["acqus"]["SW_h"])

        # zero fill to 32768 points
        try:
            self.data = ng.proc_base.zf_size(self.data, zfp)
        except ValueError:
            zfp = 2**15
            self.data = ng.proc_base.zf_size(self.data, zfp)
        # Fourier transform
        self.data = ng.proc_base.fft(self.data)
        
        # Phasecorrection
        self.data = ng.proc_autophase.ps(self.data, phc0, phc1)
        #correct ppm_scale
        self._ppmscale = np.linspace(self._ppm_l, self._ppm_r, zfp)  # np.size(self.data[0,:]))
        self.ppm = self._ppmscale

        self.data_orig = self.data

        if (xmin is not None) or (xmax is not None):
            self.data = self.set_spectral_region(xmin=xmin, xmax=xmax, mode=cropmode)
    

    def plot_fid(self, xmax=None, xmin=0, step=1):
        """
        plots every n-th(step) fid and scales the time axis with xmax and xmin
        
        :returns: figure
        """
        
        _xmax = self.TD/2 if xmax is None else xmax
        _xmin = xmin

        fig, ax = plt.subplots()
        
        for n in range(len(self.data[::step, 0])):
            ax.plot(self.fid[n, ::1].real)

        ax.set_xlim((_xmin, _xmax))
        ax.set_xlabel("data points")
        ax.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
        ax.legend([x*step for x in range(len(self.data[::step, 0]))],
                   ncol=3,
                   title="row",
                   loc=1)
        return fig

    def plot_spec(self, row, xlim=None, figsize=None, invert_xaxis=True, sharey=True):#, ppm=True):
        """
        plots row 0 and row n in the range of xmax and xmin
        
        :returns: figure
        """
        
        _max = None if xlim is None else xlim[0]
        _min = None if xlim is None else xlim[1]
        
        if type(xlim) != list:
            fig = plt.figure(figsize=figsize)
            ax = fig.add_subplot(111)
        elif type(xlim) == list:
            fig, ax = plt.subplots(ncols=len(xlim), nrows=1, figsize=figsize, sharey=sharey)
        
        _min = self.ppm[0] if xlim is None else xlim[1]
        _max = self.ppm[-1] if xlim is None else xlim[0]
        
        if type(xlim) != list:
            if type(row) == list:
                for r in row:
                    ax.plot(self.ppm, self.data[r, ::1].real, label='row %i'%r)
            else:
                ax.plot(self.ppm, self.data[row, ::1].real, label='row %i'%row)
            ax.legend()
            ax.set_xlim(xlim)
            #ax.set_title("this is row %i at %.0f V" % (row, self.eNMRraw.loc[row, "U / [V]"]))
            ax.set_xlabel("$\delta$ / ppm")
            ax.set_ylabel("intensity / a.u.")

            ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            if invert_xaxis:
                xlimits = ax.get_xlim()
                ax.set_xlim(xlimits[::-1])
        
        elif type(xlim) == list:
            for axis, xlim in zip(ax,xlim):
                if type(row) ==list:
                    for r in row:
                        axis.plot(self.ppm, self.data[r, ::1].real, label='row %i'%r)
                else:
                    axis.plot(self.ppm, self.data[row, ::1].real, label='row %i'%row)
                
                axis.set_xlim(xlim)
                #ax.set_title("this is row %i at %.0f V" % (row, self.eNMRraw.loc[row, "U / [V]"]))
                axis.set_xlabel("$\delta$ / ppm")
                axis.set_ylabel("intensity / a.u.")

                axis.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
                if invert_xaxis:
                    xlimits = axis.get_xlim()
                    axis.set_xlim(xlimits[::-1])
            #fig.legend()
        return fig
    
    def plot_spec_1d(self, xlim=None, figsize=None, invert_xaxis=True):#, ppm=True):
        """
        plots row 0 and row n in the range of xmax and xmin
        
        :returns: figure
        """
        _max = None if xlim is None else xlim[0]
        _min = None if xlim is None else xlim[1]

        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)
        
        _min = self.ppm[0] if xlim is None else xlim[1]
        _max = self.ppm[-1] if xlim is None else xlim[0]
        
        ax.plot(self.ppm, self.data[::1].real)
        
        #ax.legend()
        ax.set_xlim(xlim)
        #ax.set_title("this is row %i at %.0f V" % (row, self.eNMRraw.loc[row, "U / [V]"]))
        ax.set_xlabel("$\delta$ / ppm")
        ax.set_ylabel("intensity / a.u.")

        ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        if invert_xaxis:
            xlimits = ax.get_xlim()
            ax.set_xlim(xlimits[::-1])
        return fig

    def set_spectral_region(self, xmin, xmax, mode='absolute', ppm=True, original_data=True):
        """
        crops the data of the object set to xmin and xmax on the x-axis and returns the value
        
        mode: changes the x-scale unit
            "percent": value from 0% to 100% of the respective x-axis-length --> does not fail
            "absolute": takes the absolute values --> may fail
        
        :returns: cropped_data, cropped_ppmscale
        """
        
        if (xmin is not None) or (xmax is not None):

            if mode == "percent":
                ppmmin = xmin/100*self.data[0, :].size
                ppmmax = xmax/100*self.data[0, :].size
                if original_data:
                    self.data = self.data_orig[:, int(ppmmin):int(ppmmax)]
                else:
                    self.data = self.data[:, int(ppmmin):int(ppmmax)]
                self.ppm = self._ppmscale[int(ppmmin):int(ppmmax)]

            elif mode == "absolute":
                if ppm:
                    #print('xmin: %i' % xmin)
                    #print('xmax: %i' % xmax)
                    xmax = np.where(self._ppmscale <= xmax)[0][-1]  # if xmax is not None else -1
                    xmin = np.where(self._ppmscale >= xmin)[0][0]  # if xmin is not None else 0
                    #print('xmin: %i' % xmin)
                    #print('xmax: %i' % xmax)
                    #print(self._ppmscale)
                
                if original_data:
                    self.data = self.data_orig[:, xmin:xmax]
                else:
                    self.data = self.data[:, xmin:xmax]
                    
                self.ppm = self._ppmscale[xmin:xmax]
            else:
                raise ValueError('oops, you mistyped the mode')
            
        return self.data, self.ppm # self._ppmscale[xmin:xmax]
 

    def save_eNMRpy(self, path):
        '''
        saves the instance variables of the object in a .eNMRpy file
        '''
        
        out_dic = {"import class": self.__class__.__name__,}
        
        # reduced number of instance variables for 1D-Measurements which may be potentially saved
        if self.__class__.__name__ == "Measurement":
            instance_variables_keys = ['path', 'expno', 'dateipfad', 'alias', 'name', 'lineb', 'title_page', 'dic', 'data', 'ppm',
                                'pdata', '_ppm_r', '_ppm_l', 'TD', 'n_zf_F2', 'gamma', 'nuc']
        
        # full range of instance variables for an eNMR-Measurement
        else:
            instance_variables_keys = ['path', 'expno', 'dateipfad', 'alias', 'name', 'lineb', 'title_page', 'dic', 'data', 'ppm',
                                'pdata', '_ppm_r', '_ppm_l', 'TD', 'n_zf_F2', 'gamma', 'nuc', 'lin_res_dic', '_x_axis',
                                'dependency', 'cell_resistance', 'Delta', 'delta', 'eNMRraw', 'difflist', 'd', 'g']
        
        # check if the given path ends with .eNMRpy
        if path.split('.')[-1] != 'eNMRpy':
            # if not, change path end to .eNMRpy
            path += '.eNMRpy'
        
        # write selected instance variables in output-dictionary
        for k in instance_variables_keys:
            # if instance variable is np-array
            if type(self.__dict__[k]) == np.ndarray:
                # write array as re-readable string
                s = StringIO()
                np.savetxt(s, self.__dict__[k])
                out_dic[k] = s.getvalue()
                pass
            # write pandas DataFrames in json-format
            elif type(self.__dict__[k]) == pd.DataFrame:
                out_dic[k+'_type_pd.DataFrame'] = self.__dict__[k].to_json()
            # write everything else
            else:
                out_dic[k] = self.__dict__[k]
                
        f = open(path,'w')
        f.write(str(out_dic))
        f.close()
        return
    
    def compare_dicts(self, other):
        """
        compares the "acqus" parameters in self.dic und returns a dataframe of deviating entries
        """
        if type(self) != type(other):
            print('beware! types of both measurement objects are not equal!')
                
        klist = []
        for k in self.dic['acqus']:
            if not self.dic['acqus'][k]==other.dic['acqus'][k]:
    #             print(self.dic['acqus'][k], other.dic['acqus'][k])
                klist.append(k)
        selflist = [self.dic['acqus'][k] for k in klist]
        otherlist = [other.dic['acqus'][k] for k in klist]
        
        
        dfdiff = pd.DataFrame(index=klist)
        dfdiff['self'] = selflist
        dfdiff['other'] = otherlist
        return dfdiff
    
