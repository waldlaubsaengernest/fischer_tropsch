from fischer_tropsch.catlib_mkm_parallel import MKM, Reaction, to_dict, to_list
from copy import deepcopy

class FT_new(MKM):
    def __init__(self, eads, temperature, gas_pressures, steady_state=True, sites=None, nguess=200, max_cost=1e-8, initial_guess=None, optcyc=2000, adsorption_steady_state=False,cov_dep_eads=False,cov_ads="CO",min_rates={("C4H10","C3H8","C3H6","C2H6","C2H4","CH4"):1e-8},parallel=True):
        super().__init__(eads, temperature, gas_pressures, steady_state, sites, nguess, max_cost, initial_guess, optcyc, adsorption_steady_state,cov_dep_eads,cov_ads,min_rates=min_rates,parallel=parallel)
    
    def initialize_reactions(self):
        self.reactions = [Reaction(a="CO",b="*",ab="CO*",reverse=True),
                     Reaction(a="H*",b="H*",ab="H2",reverse=False),
                     Reaction(a="C*",b="O*",ab="CO*",reverse=False),
                     Reaction(a="O*",b="H*",ab="OH*",reverse=True),
                     Reaction(a="OH*",b="H*",ab="H2O*",reverse=True),
                     Reaction(a="H2O",b="*",ab="H2O*",reverse=False),
                     Reaction(a="C*",b="H*",ab="CH*",reverse=True),
                     Reaction(a="CH*",b="H*",ab="CH4",reverse=True,nu={"H*":3}),
                     Reaction(a="C*",b="CH*",ab="CCH*",reverse=True),
                     Reaction(a="CCH*",b="H*",ab="CCH2*",reverse=True),
                     Reaction(a="CCH2*",b="H*",ab="CCH3*",reverse=True),
                     Reaction(a="CCH3*",b="C*",ab="CCCH3*",reverse=True),
                     Reaction(a="C3H8",b="*",c="CCCH3*",d="H*",ab=None,reverse=False, nu={"H*":5})
        ]
        self.read_reaction_energies()
    
    def read_reaction_energies(self):
        reactions = to_dict(self.reactions)
        eads = self.eads
        reactions["CO+*->CO*"].er         =   eads["CO"]       - eads["CO(g)"]
        reactions["CO+*->CO*"].ea         =   0
        reactions["H2->H*+H*"].er         =   eads["H"]*2      - eads["H2(g)"]
        reactions["H2->H*+H*"].ea         =   eads["H-H"]      - eads["H2(g)"]
        reactions["CO*->C*+O*"].er        =   eads["C"]        + eads["O"]      - eads["CO"]
        reactions["CO*->C*+O*"].ea        =   eads["C-O"]                       - eads["CO"]
        reactions["O*+H*->OH*"].er        =   eads["OH"]       - eads["O"]      - eads["H"]
        reactions["O*+H*->OH*"].ea        =   eads["O-H"]      - eads["O"]      - eads["H"]
        reactions["OH*+H*->H2O*"].er      =   eads["H2O"]      - eads["OH"]     - eads["H"]
        reactions["OH*+H*->H2O*"].ea      =   eads["OH-H"]     - eads["OH"]     - eads["H"]
        reactions["H2O*->H2O+*"].er       =   eads["H2O(g)"]   - eads["H2O"]
        reactions["H2O*->H2O+*"].ea       =   0
        reactions["C*+H*->CH*"].er        =   eads["CH"]       - eads["C"]      -  eads["H"]
        reactions["C*+H*->CH*"].ea        =   eads["C-H"]      - eads["C"]      -  eads["H"]
        reactions["CH*+3H*->CH4"].er      =   eads["CH4(g)"]   - eads["CH"]     -  3*eads["H"]
        reactions["CH*+3H*->CH4"].ea      =   eads["CH3-H_s2"] - eads["CH"]     -  3*eads["H"]
        reactions["C*+CH*->CCH*"].er      =   eads["CCH"]      - eads["C"]      -  eads["CH"]
        reactions["C*+CH*->CCH*"].ea      =   eads["C-CH"]     - eads["C"]      -  eads["CH"]
        reactions["CCH*+H*->CCH2*"].er    =   eads["CCH2"]     - eads["CCH"]    -  eads["H"]
        reactions["CCH*+H*->CCH2*"].ea    =   eads["CCH-H"]    - eads["CCH"]    -  eads["H"]
        reactions["CCH2*+H*->CCH3*"].er   =   eads["CCH3"]     - eads["CCH2"]   - eads["H"]
        reactions["CCH2*+H*->CCH3*"].ea   =   eads["CCH2-H"]   - eads["CCH2"]   - eads["H"]
        reactions["CCH3*+C*->CCCH3*"].er  =   eads["CCCH3"]    - eads["CCH3"]   - eads["C"]
        reactions["CCH3*+C*->CCCH3*"].ea  =   eads["CCH3-C"]   - eads["CCH3"]   - eads["C"]
        reactions["CCCH3*+5H*->C3H8+*"].er  =   eads["C3H8(g)"] - eads["CCCH3"] - 5*eads["H"]
        reactions["CCCH3*+5H*->C3H8+*"].ea  =   0
        self.reactions = to_list(reactions)

    def prepare_mkm(self):
        self.adsorbates            = ["CO*","H*","C*","CH*","O*","OH*","H2O*","CCH*","CCH2*","CCH3*","CCCH3*"]
        self.adsorption_steady_state = True
        self.in_steady_state       = deepcopy(self.adsorbates)
        self.include_initial_guess = ["CO*","H*","C*","O*","CH*"]
        self.ads_ignore            = []
        if self.sites is None:
            self.sites= {"h":["H*"],"t":["O*","CO*","OH*","H2O*"],"c":["C*"],"c2":["CH*","CCH*","CCH2*","CCH3*","CCCH3*"]}

    def get_final_rates_and_selectivities(self):
        reaction_rates = dict()
        reaction_rates["CO"]   = - self.reaction_rates["CO+*->CO*"]["net"]
        reaction_rates["H2"]   = -1/2 *self.reaction_rates["H2->H*+H*"]["net"]
        reaction_rates["C3H8"] = self.reaction_rates["CCCH3*+5H*->C3H8+*"]["net"]
        reaction_rates["CH4"]  = self.reaction_rates["CH*+3H*->CH4"]["net"]
        self.reaction_rates_nomult = deepcopy(reaction_rates)
        reaction_rates["C3H8"] *=3
        eps = 10**(-14)
        self.selectivity       = (max(eps,reaction_rates["C3H8"]))/(max(eps,reaction_rates["C3H8"]) + max(eps,reaction_rates["CH4"]))
        self.reaction_rates    = reaction_rates
