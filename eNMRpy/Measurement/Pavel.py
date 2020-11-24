from .eNMR_Methods import _eNMR_Methods
import pandas as pd
import xml.etree.ElementTree as etree

class Pavel(_eNMR_Methods):
    '''
    This is the subsubclass of Masurement() and subclass of eNMR_Methods specialised to process data obtained from the experimental Swedish from Pavel set-up
    the voltage list is valculated from the vd-values

    path:
        relative or absolute path to the measurements folder
    expno:
        the to the experiment number corresponding EXPNO
    dependency:
        'U': voltage dependent eNMR measurement
        'G': fieldgradient dependent eNMR measurement
             -->reads only the first value in the vd-List to calculate the set voltage

    alias:
        Here you can place an individual name relevant for plotting. If None, the path is taken instead.
    lineb:
        setting a standard-value for the linebroadening.
    d:
        electrode_distance
        
    Uconst:
        constant voltage value set for gradient dependent measurements (dependency="G")
    '''
    def __init__(self, path, expno, dependency='U', alias=None, lineb=.3, d=2.2e-2, cell_resistance=None, Uconst=None):
        
        self.dependency = dependency
        self.cell_resistance = cell_resistance
        
        super().__init__(path, expno, lineb=lineb, alias=alias)
        
        #self._x_axis = {"G": 'g in T/m', "U": 'U / [V]'}[dependency.upper()]
        # import the diffusion parameters
        diffpar = etree.parse(self.dateipfad+'/diff.xml')
        root = diffpar.getroot()
        self.Delta = float(root.findall('DELTA')[0].text)*1e-3
        self.delta = float(root.findall('delta')[0].text)*1e-3  #in Seconds
        print('The diffusion parameters were read from the respective .XML!')
        
        try:
            self.vdList = pd.read_csv(self.dateipfad+"/vdlist",
                                      names=["vd"]).loc[:len(self.data[:, 0])-1]
        except IndexError:
            raise IndexError('Your data is %i-dimensional instead of 2-dimensional. This may be an issue with Topspin. Most likely your measurement was aborted.'%self.data.ndim)
        except:
            raise FileNotFoundError('no VD-List found!')
            #print('no vdList found, generated ones list instead')
            #self.vdList = pd.DataFrame(np.ones((len(self.data[:, 0]), 1)),
                                       #columns=["vd"])
        self.eNMRraw = self.vdList

        #self.vdList["U / [V]"] = hier die Konversion von vdlist zu Spannungsliste
        try:
            self.difflist = pd.read_csv(self.dateipfad+"/gradlist",
                                    names=["g in T/m"])*0.01
        except:
            print('gradlist not found. difflist imported instead')
            self.difflist = pd.read_csv(self.dateipfad+"/difflist",
                                    names=["g in T/m"])*0.01
        self.eNMRraw["g in T/m"] = self.difflist
        
        self.d = d
        self.g = self.eNMRraw["g in T/m"][0]
        

        if self.dependency.upper() == "G":
            self.uInk = None
            # was introduced for vdlists with only a single entry but works regardless of the number of entries if dependency=='G', since the voltage is kept constant
            #if self.vdList.shape[0] == 1: 
            self.eNMRraw = self.difflist
            self.eNMRraw['vd'] = float(self.vdList.iloc[0, 0][:-1])
            self.eNMRraw['U / [V]'] = [
            0 if (self.eNMRraw.loc[i,'vd'] <= 0.6)
            else 
            n 
            for i, n in enumerate(self.eNMRraw['vd']*5)]
            
        else:
            # converts the vd-List
            for i, n in enumerate(self.eNMRraw['vd']):
                self.eNMRraw.loc[i, 'vd_temp'] = float(n[:-1])
            # calculates the applied Voltages        
        
        if (self.dependency.upper() == "U") and (self.dic['acqus']['PULPROG'] == 'py_ENMR_DSTE_3.2_AL'): # checks additionally for the used PULPROG
            self.eNMRraw[self._x_axis] = [
                0 if (self.eNMRraw.loc[i,'vd_temp'] <= 0.6)
                else 
                n if i%2==0 
                else 
                n*-1 # voltage reversal inherent to the PULPROG
                for i, n in enumerate(self.eNMRraw['vd_temp']*5)]
        
            self.uInk = self.eNMRraw['U / [V]'][0] - self.eNMRraw['U / [V]'][1] 
            if self.uInk == 0:
                self.uInk = self.eNMRraw['U / [V]'][0] - self.eNMRraw['U / [V]'][2] 
            if self.uInk < 0:
                self.uInk *= -1
        
        elif self.dependency.upper() == "I":
            self.uInk = None
            self.eNMRraw[self._x_axis] = [
                    0 if (self.eNMRraw.loc[i,'vd_temp'] <= 0.6)
                    else 
                    n if i%2==0 
                    else 
                    n*-1
                    for i, n in enumerate(self.eNMRraw['vd_temp'])
                ]
        
        elif self.dependency.upper() == "RI":
            self.uInk = None
            self.eNMRraw[self._x_axis] = [
                    0 if (self.eNMRraw.loc[i,'vd_temp'] <= 0.6)
                    else 
                    n if i%2==0 
                    else 
                    n*-1
                    for i, n in enumerate(self.eNMRraw['vd_temp'])
                ]
            self.uInk = self.eNMRraw['RI / V'][0] - self.eNMRraw['RI / V'][1] 
            if self.uInk == 0:
                self.uInk = self.eNMRraw['RI / V'][0] - self.eNMRraw['RI / V'][2] 
            if self.uInk < 0:
                self.uInk *= -1
            # calculation of the Voltage from cell resistance and Current /1000 because of mA
            self.eNMRraw[self._x_axis] *= self.cell_resistance/1000

         
    def plot_spec(self, row, xlim=None, figsize=None, invert_xaxis=True, sharey=True):#, ppm=True):
        from .Juergen1 import Juergen1 as eNMR_Measurement
        return eNMR_Measurement.plot_spec(self, row, xlim, figsize, invert_xaxis, sharey)#, ppm=True):
