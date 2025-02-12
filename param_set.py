"""
Utility functions and classes (including default parameters).
Author: Sara Mathieson, Rebecca Riley
Date: 9/27/22
"""

# python imports
import numpy as np
from scipy.stats import norm

class Parameter:
    """
    Holds information about evolutionary parameters to infer.
    Note: the value arg is NOT the starting value, just used as a default if
    that parameter is not inferred, or the truth when training data is simulated
    """

    def __init__(self, value, min, max, name):
        self.value = value
        self.min = min
        self.max = max
        self.name = name
        self.proposal_width = (self.max - self.min)/15 # heuristic

    def __str__(self):
        s = '\t'.join(["NAME", "VALUE", "MIN", "MAX"]) + '\n'
        s += '\t'.join([str(self.name), str(self.value), str(self.min),
            str(self.max)])
        return s

    def start(self):
        if isinstance(self.value, list):
            main_value = np.random.uniform(self.min, self.max)
            scnd_value = np.random.uniform(self.min, 1 - main_value)
            return [main_value, scnd_value]
        # random initialization
        return np.random.uniform(self.min, self.max)

    def start_range(self):
        if self.name == 'propJOINT':
            raise IOError
        start_min = np.random.uniform(self.min, self.max)
        start_max = np.random.uniform(self.min, self.max)
        if start_min <= start_max:
            return [start_min, start_max]
        return self.start_range()

    def fit_to_range(self, value):
        if isinstance(value, list):
            new_list = []
            new_list.append(max(min(value[0], self.max), self.min))
            new_list.append(max(min(value[1], 1 - new_list[0]), self.min))
            return new_list
        
        value = min(value, self.max)
        return max(value, self.min)

    def proposal(self, curr_value, multiplier):
        if multiplier <= 0: # last iter
            return curr_value

        if isinstance(curr_value, list):
            new_val = norm(curr_value[0], self.proposal_width * multiplier).rvs()
            new_comp = norm(curr_value[1], self.proposal_width * multiplier).rvs()
            new_value = [new_val, new_comp]
        # normal around current value (make sure we don't go outside bounds)
        else:
            new_value = norm(curr_value, self.proposal_width * multiplier).rvs()
        new_value = self.fit_to_range(new_value)
        # if the parameter hits the min or max it tends to get stuck
        if new_value == curr_value or new_value == self.min or new_value == self.max:
            return self.proposal(curr_value, multiplier) # recurse
        elif isinstance(new_value, list) and (new_value[0] == self.min or new_value[0] == self.max):
            return self.proposal(curr_value, multiplier)
        else:
            return new_value

    def proposal_range(self, curr_lst, multiplier):
        if self.name == 'propJOINT':
            raise IOError

        new_min = self.fit_to_range(norm(curr_lst[0], self.proposal_width *
            multiplier).rvs())
        new_max = self.fit_to_range(norm(curr_lst[1], self.proposal_width *
            multiplier).rvs())
        if new_min <= new_max:
            return [new_min, new_max]
        return self.proposal_range(curr_lst, multiplier) # try again

class ParamSet:

    def __init__(self):

        # default Ne and reco and mut
        self.Ne = Parameter(10000, 1000, 30000, "Ne")
        self.reco = Parameter(1.25e-8, 1e-9, 1e-7, "reco")
        self.mut = Parameter(1.25e-8, 1e-9, 1e-7, "mut")

        # IM
        self.N_anc = Parameter(15000, 1000, 25000, "N_anc")
        self.T_split = Parameter(2000, 500, 20000, "T_split")
        self.mig = Parameter(0.05, -0.2, 0.2, "mig")

        # IM and exp
        self.N1 = Parameter(9000, 1000, 30000, "N1")
        self.N2 = Parameter(5000, 1000, 30000, "N2")

        # exp
        self.growth = Parameter(0.005, 0.0, 0.05, "growth")

        # ooa2
        self.N3 = Parameter(12000, 1000, 30000, "N3")
        self.T1 = Parameter(2000, 1500, 5000, "T1")
        self.T2 = Parameter(350, 100, 1500, "T2")
        
        # ADMIX
        self.NeAFR = Parameter(14474, 4000, 30000, "NeAFR") 
        self.NeEUR = Parameter(34039, 10000, 50000, "NeEUR") 
        self.NeEAS = Parameter(45852, 10000, 60000, "NeEAS") 
        self.NeADMIX = Parameter(30000, 20000, 70000, "NeADMIX")
        self.growth_EUR = Parameter(0.0038, 0.0, 0.05, "growth_EUR") 
        self.growth_EAS = Parameter(0.0048, 0.0, 0.05,  "growth_EAS")
        self.growth_ADMIX = Parameter(0.05, 0.0, 0.2, "growth_EUR") 
        self.mAFEU = Parameter(0, -0.000025, 0.00025, "mAFEU")
        self.mAFEA = Parameter(0, -0.2, 0.2, "mAFEA")
        self.mEUEA = Parameter(0, -0.0000078, 0.00002, "mEUEA")
        self.Tadmix = Parameter(12, 5, 30, "Tadmix") 
        
        self.propJOINT = Parameter([0.2, 0.5], 0, 1, "propJOINT")
	# older msprime version 
        self.radmix= Parameter(.05, -0.05, 0.08, "radmix") # growth rate of admixed population


        # ooa3
        self.N_A = Parameter(7300, 1000, 30000, "N_A")
        self.N_B = Parameter(2100, 1000, 20000, "N_B")
        self.N_AF = Parameter(12300, 1000, 40000, "N_AF")
        self.N_EU0 = Parameter(1000, 100, 20000, "N_EU0")
        self.N_AS0 = Parameter(510, 100, 20000, "N_AS0")
        self.r_EU = Parameter(0.004, 0.0, 0.05, "r_EU")
        self.r_AS = Parameter(0.0055, 0.0, 0.05, "r_AS")
        self.T_AF = Parameter(8800, 8000, 15000, "T_AF")
        self.T_B = Parameter(5600, 2000, 8000, "T_B")
        self.T_EU_AS = Parameter(848, 100, 2000, "T_EU_AS")
        self.m_AF_B = Parameter(25e-5, 0.0, 0.01, "m_AF_B")
        self.m_AF_EU = Parameter(3e-5, 0.0,  0.01, "m_AF_EU")
        self.m_AF_AS = Parameter(1.9e-5, 0.0, 0.01, "m_AF_AS")
        self.m_EU_AS = Parameter(9.6e-5, 0.0, 0.01, "m_EU_AS")

    def update(self, names, values):
        """Based on generator proposal, update desired param values"""
        assert len(names) == len(values)

        for j in range(len(names)):
            param = names[j]

            # credit: Alex Pan (https://github.com/apanana/pg-gan)
            attr = getattr(self, param)
            if attr == None:
                sys.exit(param + " is not a recognized parameter.")
            else:
                attr.value = values[j]
