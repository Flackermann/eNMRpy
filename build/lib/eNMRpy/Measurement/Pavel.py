from .eNMR_Methods import _eNMR_Methods
import pandas as pd


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

    alias:
        Here you can place an individual name relevant for plotting. If None, the path is taken instead.
    linebroadening:
        setting a standard-value for the linebroadening.
    '''
    def __init__(self, path, expno, dependency='U', alias=None, linebroadening=5, electrode_distance=2.2e-2, cell_resistance=None):
        
        self.dependency = dependency
        self.cell_resistance = cell_resistance
        
        super().__init__(path, expno, linebroadening=linebroadening, alias=alias)
        
        #self._x_axis = {"G": 'g in T/m', "U": 'U / [V]'}[dependency.upper()]
        # import the diffusion parameters
        import xml.etree.ElementTree as etree
        diffpar = etree.parse(self.dateipfad+'/diff.xml')
        root = diffpar.getroot()
        self.Delta = float(root.findall('DELTA')[0].text)*1e-3
        self.delta = float(root.findall('delta')[0].text)*1e-3  #in Seconds
        print('The diffusion parameters were read from the respectie .XML!')
        
        try:
            self.vdList = pd.read_csv(self.dateipfad+"/vdlist",
                                      names=["vd"]).loc[:len(self.data[:, 0])-1]
        except:
            print('no vdList found, generated ones list instead')
            self.vdList = pd.DataFrame(np.ones((len(self.data[:, 0]), 1)),
                                       columns=["vd"])
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
        
        self.d = electrode_distance
        self.g = self.eNMRraw["g in T/m"][0]
        

        # converts the vd-List
        for i, n in enumerate(self.eNMRraw['vd']):
            self.eNMRraw.loc[i, 'vd_temp'] = float(n[:-1])
        # calculates the applied Voltages

        if self.dependency.upper() == "U":
            self.eNMRraw[self._x_axis] = [
                0 if (self.eNMRraw.loc[i,'vd_temp'] <= 0.6)
                else 
                n if i%2==0 
                else 
                n*-1
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

#class Idependent(Pavel):
    #def __init__(self, path, expno, dependency='U', alias=None, linebroadening=5, electrode_distance=2.2e-2):
        #Pavel.__init__(path, expno, dependency='U', alias=None, linebroadening=5, electrode_distance=2.2e-2)
        
        #self.eNMR['I'] = None#vd=Funktion
        
        #self.eNMR.drop('U / [V]', inplace=True)
 
