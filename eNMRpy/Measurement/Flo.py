from .eNMR_Methods import _eNMR_Methods
import pandas as pd
import xml.etree.ElementTree as etree
import re
import numpy as np


class Flo(_eNMR_Methods):
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
    lineb:
        setting a standard-value for the linebroadening.
    d:
        electrode_distance
    '''
    def __init__(self, path, expno, dependency='U', alias=None, lineb=.3, d=2.2e-2, cell_resistance=None):
        
        self.dependency = dependency
        self.cell_resistance = cell_resistance
        
        super().__init__(path, expno, lineb=lineb, alias=alias)
        
        self._x_axis = {"G": 'g in T/m', "U": 'U / [V]', "I": "I / mA"}[dependency.upper()]
        # import the diffusion parameters
        diffpar = etree.parse(self.dateipfad+'/diff.xml')
        root = diffpar.getroot()
        self.Delta = float(root.findall('DELTA')[0].text)*1e-3
        self.delta = float(root.findall('delta')[0].text)*1e-3  #in Seconds
        print('The diffusion parameters were read from the respective .XML!')

        if path[-1] == '/':
            pulseprogram = open(path+str(expno)+'/pulseprogram')
        else:
            pulseprogram = open(path+'/'+str(expno)+'/pulseprogram')
            
        self.pulseprogram = pulseprogram.read()

        bitregex = r"(define list.*bit.*|define list.*pol.*)"
        vlistregex = r".*Voltage\sset\sList.*= \n;.*]"
        ilistregex = r".*Current\sset\sList.*= \n;.*]"
        
        # list, reading all lines with bit lists in the pulse program
        rawlist = re.findall(bitregex, self.pulseprogram)
        rawvlist = re.findall(vlistregex, self.pulseprogram) 
        
        # check vor const U or I mode
        if len(rawvlist) == 0:
            rawvlist = re.findall(ilistregex, self.pulseprogram) 
            self.dependency = 'I'
            self._x_axis = "I / mA"
            print('dependency changed to "I" --> const current mode')
        
        
        rawvlist = rawvlist[0].split('= \n;')
        vlist = eval(rawvlist[1])

        # array of integers generated from the rawlist. Indexing like [bit, voltagestep]
        bitarray = np.array([[int(i) for i in re.findall('{.*}', rawlist[j])[0][1:-1].split(' ')] for j in range(len(rawlist))])


        def byte_to_int(bitarray, row):
            #converting the array into the correct string
            bitstring = str(bitarray[:,row])[1:-1].replace(' ', '')
            
            # check the last bit for the polarity
            if bitstring[-1] == '1':
                polarity = 1
            elif bitstring[-1] == '0':
                polarity = -1
            # transformation of the bitstring minus polarity, which is the last bit, with a base of 2
            intvar = int(bitstring[:-1][::-1], 2)
            
            if self.dependency.upper() == 'U':
                return round(polarity*intvar*200/255, 2)
            
            if self.dependency.upper() == 'I':
                return round(polarity*intvar*50/255, 2)
            


        ulist = [byte_to_int(bitarray, i) for i in range(len(bitarray[0]))]
        
        if ulist == vlist:
            pass
        else:
            raise ValueError('The decoded voltage list does not match the endcoding voltage list! Revisit your pulse program!\n {} \n{}'.format(ulist, vlist))
        
        
        self.eNMRraw = pd.DataFrame(ulist, columns=[self._x_axis])
        
        
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
        

        ## converts the vd-List
        #for i, n in enumerate(self.eNMRraw['vd']):
            #self.eNMRraw.loc[i, 'vd_temp'] = float(n[:-1])
        ## calculates the applied Voltages

        #if self.dependency.upper() == "U":
            #self.eNMRraw[self._x_axis] = [
                #0 if (self.eNMRraw.loc[i,'vd_temp'] <= 0.6)
                #else 
                #n if i%2==0 
                #else 
                #n*-1
                #for i, n in enumerate(self.eNMRraw['vd_temp']*5)]
        
            #self.uInk = self.eNMRraw['U / [V]'][0] - self.eNMRraw['U / [V]'][1] 
            #if self.uInk == 0:
                #self.uInk = self.eNMRraw['U / [V]'][0] - self.eNMRraw['U / [V]'][2] 
            #if self.uInk < 0:
                #self.uInk *= -1
        
        #elif self.dependency.upper() == "I":
            #self.uInk = None
            #self.eNMRraw[self._x_axis] = [
                    #0 if (self.eNMRraw.loc[i,'vd_temp'] <= 0.6)
                    #else 
                    #n if i%2==0 
                    #else 
                    #n*-1
                    #for i, n in enumerate(self.eNMRraw['vd_temp'])
                #]
        
        #elif self.dependency.upper() == "RI":
            #self.uInk = None
            #self.eNMRraw[self._x_axis] = [
                    #0 if (self.eNMRraw.loc[i,'vd_temp'] <= 0.6)
                    #else 
                    #n if i%2==0 
                    #else 
                    #n*-1
                    #for i, n in enumerate(self.eNMRraw['vd_temp'])
                #]
            #self.uInk = self.eNMRraw['RI / V'][0] - self.eNMRraw['RI / V'][1] 
            #if self.uInk == 0:
                #self.uInk = self.eNMRraw['RI / V'][0] - self.eNMRraw['RI / V'][2] 
            #if self.uInk < 0:
                #self.uInk *= -1
            ## calculation of the Voltage from cell resistance and Current /1000 because of mA
            #self.eNMRraw[self._x_axis] *= self.cell_resistance/1000
            
    def plot_spec(self, row, xlim=None, figsize=None, invert_xaxis=True, sharey=True):#, ppm=True):
        from .Juergen1 import Juergen1 as eNMR_Measurement
        return eNMR_Measurement.plot_spec(self, row, xlim, figsize, invert_xaxis, sharey)#, ppm=True):
