from .eNMR_Methods import _eNMR_Methods
import matplotlib.pyplot as plt

class Juergen1(_eNMR_Methods):
    '''
    This is the subsubclass of Masurement() and subclass of eNMR_Methods specialised to process data obtained from the experimental Schönhoff set-up
    
    path:
        relative or absolute path to the measurements folder
    measurement:
        the to the experiment corresponding EXPNO
    alias:
        Here you can place an individual name relevant for plotting. If None, the path is taken instead.    
    Uink:
        voltage increment. Usually extracted from the title file if defined with e.g. "Uink = 10V"
        If Uink cannot be found or is wrong it can be entered manually during the data import.
        The voltage list is calculated from the voltage increment and the vc list when the incrementation loop is used in the pulse program
    dependency:
        'U': voltage dependent eNMR measurement
        'G': fieldgradient dependent eNMR measurement

    linebroadening:
        setting a standard-value for the linebroadening.
    '''
    def __init__(self, path, expno, Uink=None, dependency="U", alias=None, linebroadening=0.5, electrode_distance=2.2e-2):
        Measurement.__init__(self, path, expno, linebroadening=linebroadening, alias=alias)
        self.dependency = dependency.upper()
        
        self._x_axis = {"U": "U / [V]",
                   "G": "g in T/m",
                   "I": "I / mA",
                   'RI': 'RI / V'
                   }[self.dependency.upper()]
        
        #self._x_axis = {"G": "g in T/m",
                        #"U": "U / [V]"}[self.dependency.upper()]
        
        if dependency.upper() == "U":
            try:
                # takes the title page to extract the volt increment
                title = open(self.dateipfad+"/pdata/1/title").read()
                # gets the voltage increment using a regular expression
                #uimport = findall('[U|u]in[k|c]\s*=?\s*\d+', title)[0]
                uimport = findall('[U|u]in[k|c]\s*=+\s*\d+', title)[0]
                self.uInk = int(findall('\d+', uimport)[0])
            except ValueError:
                print('no volt increment found\nyou may want to put it in manually')
                self.uInk = Uink
            except IndexError:
                print('No Uink found! May not be an eNMR experiment.')
                self.uInk = Uink
                
        elif dependency.upper() == "G":
            try:
                # takes the title page to extract the volt increment
                title = open(self.dateipfad+"/pdata/1/title").read()
                # gets the voltage increment using a regular expression
                uimport = findall('[U|u]\s*=?\s*\d+', title)[0]
                self.uInk = int(findall('\d+', uimport)[0])

            except ValueError:
                print('no volt increment found\nyou may want to put it in manually')
                self.uInk = Uink
            except IndexError:
                print('No Uink found! May not be an eNMR experiment.')
                self.uInk = Uink # Uinktext

        if self.dependency.upper() == "U":
            try:
                self.vcList = pd.read_csv(self.dateipfad+"/vclist",
                                          names=["vc"]).loc[:len(self.data[:, 0])-1]
            except:
                print("There is a Problem with the VC-list or you performed a gradient dependent measurement")
        elif self.dependency.upper() == "G":
            self.vcList = pd.DataFrame(np.ones((len(self.data[:, 0]), 1)),
                                       columns=["vc"])
        else:
            print("The dependency is not properly selected, try again!")

        self.difflist = pd.read_csv(self.dateipfad+"/difflist",
                                    names=["g in T/m"])*0.01
        
        if Uink is not None:
            self.uInk = Uink
            
        self.vcList["U / [V]"] = [self.vcList["vc"][n]/2*self.uInk if self.vcList["vc"][n] % 2 == 0
                                  else (self.vcList["vc"][n]+1)/2*self.uInk*-1
                                  for n in range(len(self.data[:, 0]))]
        
        # try to open phase data, otherwise create new
        try:
            self.eNMRraw = pd.read_csv(self.path+"phase_data_"+self.expno+".csv",
                                       index_col=0, sep=" ")
            # --> update voltage list
            self.eNMRraw["U / [V]"] = self.vcList["U / [V]"]
        except:
            print("eNMRraw was missing and is generated")
            self.vcList["ph0"] = np.zeros(len(self.data.real[:, 0]))
            self.eNMRraw = self.vcList
        finally:
            self.eNMRraw["g in T/m"] = self.difflist
        
        self.p1 = self.dic["acqus"]["P"][1]
        self.d1 = self.dic["acqus"]["D"][1]
        
        try:
            # import of diffusion parameters for newer Spectrometers
            import xml.etree.ElementTree as etree
            diffpar = etree.parse(self.dateipfad+'/diff.xml')
            root = diffpar.getroot()
            self.Delta = float(root.findall('DELTA')[0].text)*1e-3
            self.delta = float(root.findall('delta')[0].text)*1e-3  # it should be read as in microseconds at this point due to bruker syntax
            print('The diffusion parameters were read from the respectie .XML!')
        except:
            # determination of the diffusion parameters for Emma
            self._d2 = self.dic["acqus"]["D"][2]
            self._d5 = self.dic["acqus"]["D"][5]
            self._d9 = self.dic["acqus"]["D"][9]
            self._d11 = self.dic["acqus"]["D"][11]
            self._p19, self._p18, self._p17 = self.dic["acqus"]["P"][19],\
                                            self.dic["acqus"]["P"][18],\
                                            self.dic["acqus"]["P"][17]
            print('That did not work. Your data is from an old spectrometer!')
            # calculating usable parameters
            self.delta = self._p17+self._p18
            self._Delta_1 = 0.001*(self._p17*2+self._p18)+(self._d2+self._d9+self._d5+self._d11)*1000+0.001*self.p1+self._d11
            self._Delta_2 = 0.001*(self._p17*2+self._p18)+(self._d2+self._d9+self._d5+self._d11)*1000+0.001*self.p1*2
            self._spoiler = (self._d11+self._p17+self._p19+self._p17)*0.001+self._d2*1000
            self.Delta = self._Delta_1+self._Delta_2+2*self._spoiler
            self.Delta *=1e-3
            self.delta *=1e-6
            

        # Elektrodenabstand in m
        self.d = electrode_distance
        self.g = self.eNMRraw["g in T/m"][0]
    
    def plot_spec(self, row, xlim=None, figsize=None, invert_xaxis=True, sharey=True):#, ppm=True):
        """
        plots row 0 and row n in the range of xmax and xmin
        
        :returns: figure
        """

        _max = None if xlim is None else xlim[0]
        _min = None if xlim is None else xlim[1]
        
        if type(xlim) is not list:
            fig = plt.figure(figsize=figsize)
            ax = fig.add_subplot(111)
        elif type(xlim) == list:
            fig, ax = plt.subplots(ncols=len(xlim), nrows=1, figsize=figsize, sharey=sharey)
        
        _min = self.ppm[0] if xlim is None else xlim[1]
        _max = self.ppm[-1] if xlim is None else xlim[0]
        
        def plot(r, axes=ax):
            
            unit = {'U': 'V', 'G': 'T/m', 'I': 'mA', 'RI': 'V'}[self.dependency.upper()]
            
            axes.plot(self.ppm, self.data[r, ::1].real, label='row %i, %i %s'%(r, self.eNMRraw[self._x_axis].iloc[r], unit))
        
        if type(xlim) is not list:
            if type(row) ==list:
                for r in row:
                    plot(r)
            else:
                plot(row)
        
            ax.set_xlim(xlim)
            #ax.set_title("this is row %i at %.0f V" % (row, self.eNMRraw.loc[row, "U / [V]"]))
            ax.set_xlabel("$\delta$ / ppm")
            ax.set_ylabel("intensity / a.u.")

            ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            if invert_xaxis:
                xlimits = ax.get_xlim()
                ax.set_xlim(xlimits[::-1])
            
            ax.legend()
        
        elif type(xlim) == list:
            for axis, xlim in zip(ax,xlim):
                if type(row) ==list:
                    for r in row:
                        plot(r, axis)
                else:
                    plot(row, axis)
                
                axis.set_xlim(xlim)
                #ax.set_title("this is row %i at %.0f V" % (row, self.eNMRraw.loc[row, "U / [V]"]))
                axis.set_xlabel("$\delta$ / ppm")
                axis.set_ylabel("intensity / a.u.")
                axis.legend()
                axis.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
                if invert_xaxis:
                    xlimits = axis.get_xlim()
                    axis.set_xlim(xlimits[::-1])
            ax[-1].legend()
            #fig.legend()
        return fig
    
 
