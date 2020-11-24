## coding: utf-8
# This is the eNMR class. THE library for the evaluation of bruker-eNMR-spectra based on the VC-PowerSource of the
# Schönhoff working group.
# It works with the volt-increment method which calculates the respective voltage with the VC-list
# Further Implementation can be asked for at f_schm52@wwu.de
from .eNMR_Methods import _eNMR_Methods
import matplotlib.pyplot as plt
from .base import Measurement
from re import findall
import numpy as np
import pandas as pd
     

#class eNMR_Emma(eNMR_Measurement):
class eNMR_Emma(_eNMR_Methods):
    #'''
    #This is the subsubclass of Masurement() and subclass of eNMR_Methods specialised to process data obtained from the experimental Schönhoff set-up
    
    #path:
        #relative or absolute path to the measurements folder
    #measurement:
        #the to the experiment corresponding EXPNO
    #alias:
        #Here you can place an individual name relevant for plotting. If None, the path is taken instead.    
    #Uink:
        #voltage increment. Usually extracted from the title file if defined with e.g. "Uink = 10V"
        #If Uink cannot be found or is wrong it can be entered manually during the data import.
        #The voltage list is calculated from the voltage increment and the vc list when the incrementation loop is used in the pulse program
    #dependency:
        #'U': voltage dependent eNMR measurement
        #'G': fieldgradient dependent eNMR measurement

    #linebroadening:
        #setting a standard-value for the linebroadening.
    #'''
    def __init__(self, path, expno, Uink=None, dependency="U", alias=None, lineb=0.5, electrode_distance=2.2e-2):
        Measurement.__init__(self, path, expno, lineb=lineb, alias=alias)
        self.dependency = dependency.upper()
        
        self._x_axis = {"U": "U / [V]",
                   "G": "g in T/m",
                   "I": "I / mA",
                   'RI': 'RI / V'
                   }[self.dependency.upper()]
        
        #self._x_axis = {"G": "g in T/m",
                        #"U": "U / [V]"}[self.dependency.upper()]
        
        self.difflist = pd.read_csv(self.dateipfad+"/difflist",
                            names=["g in T/m"])*0.01
        
        self.vcList = pd.DataFrame()
        
        if self.dic['acqus']['PULPROG'][-3:] == 'var':
            polarity = 1
            print('this is a regular measurement! (non-_pol)')
        elif self.dic['acqus']['PULPROG'][-3:] == 'pol':
            polarity = -1
            print('this is a _pol-Measurement!')
        else:
            print("no var or pol PULPROG")
    
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
            
            self.vcList["U / [V]"] = [i*self.uInk*polarity for i in range(len(self.difflist))]
                
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
                
            #if Uink is not None:
                #self.uInk = Uink
            
            self.vcList["U / [V]"] = [self.uInk*polarity for i in range(len(self.difflist))]
            #self.vcList["U / [V]"] = [self.vcList["vc"][n]/2*self.uInk if self.vcList["vc"][n] % 2 == 0
                                    #else (self.vcList["vc"][n]+1)/2*self.uInk*-1
                                    #for n in range(len(self.data[:, 0]))]

        #if self.dependency.upper() == "U":
            #try:
                #self.vcList = pd.read_csv(self.dateipfad+"/vclist",
                                          #names=["vc"]).loc[:len(self.data[:, 0])-1]
                
            #except:
                #print("There is a Problem with the VC-list or you performed a gradient dependent measurement")
        #elif self.dependency.upper() == "G":
            #self.vcList = pd.DataFrame(np.ones((len(self.data[:, 0]), 1)),
                                       #columns=["vc"])
        #else:
            #print("The dependency is not properly selected, try again!")

        #self.difflist = pd.read_csv(self.dateipfad+"/difflist",
                                    #names=["g in T/m"])*0.01
        
            #if Uink is not None:
                #self.uInk = Uink
                
            #self.vcList["U / [V]"] = [self.vcList["vc"][n]/2*self.uInk if self.vcList["vc"][n] % 2 == 0
                                    #else (self.vcList["vc"][n]+1)/2*self.uInk*-1
                                    #for n in range(len(self.data[:, 0]))]
            
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
    
    def __add__(self, other):
        
        for obj in [self, other]:
            for k in obj.eNMRraw.columns:
                if k[:2] == 'ph':
                    obj.eNMRraw[k] -= obj.eNMRraw.loc[0, k]
                    print('%s normalized to 0V'%k)
                else:
                    pass
        self.eNMRraw = self.eNMRraw.append(other.eNMRraw)
        self.eNMRraw.sort_values('U / [V]', inplace=True)
        return self
