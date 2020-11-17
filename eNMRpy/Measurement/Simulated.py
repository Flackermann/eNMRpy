from . import eNMR_Methods

class Simulated(eNMR_Methods):
    '''
    sublass of eNMR_Methods and Measurement for the import of simulated data via the SpecSim class in the Phasefitting Module
    '''

    def __init__(self, sim, alias=None):
        self.g = sim.params.spec_par['g']
        self.d = sim.params.spec_par['d']
        self.Delta = sim.params.spec_par['Delta']
        self.delta = sim.params.spec_par['delta']
        self.dependency = 'U'
        self.eNMRraw = sim.eNMRraw
        self.data = sim.data
        self.data_orig = sim.data
        self.ppm = sim.ppm
        self.eNMRraw['g in T/m'] = self.g

        # the gamma_values in rad/Ts
        gamma_values = {'1H':26.7513e7,
                        '7Li': 10.3962e7,
                        '19F': 25.1662e7}
        self.gamma = gamma_values[sim.params.spec_par['NUC']]
        # Umrechnung von rad in °
        self.gamma = self.gamma/2/np.pi*360
        
        self.lin_res_dic = {}
        self.alias = alias
        self.path = 'simulated data/'
        self.expno = '0'
        self.uInk = self.eNMRraw.loc[1, 'U / [V]'] - self.eNMRraw.loc[0, 'U / [V]']
        self._x_axis = 'U / [V]'
        #self._x_axis = {"U": "U / [V]",
            #"G": "g in T/m"
            #}[self.dependency.upper()]
##################################
 
