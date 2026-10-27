import numpy as np
class GasReference:
    def __init__(self,ref="(A)",reaction="CO_oxidation"):
        self.reaction=reaction
        self.ref = ref
        self.get_delta_e()
        self.set_reference_energies()
    
    def get_delta_e(self):
        self.delta_e = {}
        if self.reaction =="CO_oxidation":
            self.delta_e["2CO(g)+O2(g)->CO2(g)"] = -5.82
            self.delta_e["2CO(g)->2C_gr+O2(g)"]  = 2.21
            self.delta_e["CO(g)+3H2(g)->CH4(g)+H2O(g)"]  = -2.77
            self.delta_e["2H2(g)+O2(g)->2H2O(g)"]  = -5.43
        if self.reaction =="NO_reduction":
            self.delta_e["2H2(g)+O2(g)->2H2O(g)"]  = -5.43

    def set_reference_energies(self):
        if self.reaction =="CO_oxidation":
            if self.ref == "(A)": # Cref = CO; Oref = O2
                e_co  = 0
                e_o2  = 0
            elif self.ref == "(B)": # Cref = CO2; Oref = O2
                e_co  = -self.delta_e["2CO(g)+O2(g)->CO2(g)"]/2
                e_o2  = 0
            elif self.ref == "(C)": # Cref = C_gr; Oref = O2
                e_co  = -self.delta_e["2CO(g)->2C_gr+O2(g)"]/2
                e_o2  = 0
            elif self.ref == "(D)": # Cref = CH4; Oref = H2O --> E_CO = 2.77; E_O = 2.715; E_C = 0.062
                e_co  = -self.delta_e["CO(g)+3H2(g)->CH4(g)+H2O(g)"]
                e_o2  = -self.delta_e["2H2(g)+O2(g)->2H2O(g)"]
            self.eref = dict()
            self.eref["CO_g"]    = e_co
            self.eref["O2_g"]    = e_o2
        if self.reaction =="NO_reduction":
            if self.ref == "(A)": # Cref = CO; Oref = O2
                e_h2  = 0
                e_n2  = 0
                e_no  = 0
            self.eref["H2_g"]    = e_h2
            self.eref["N2_g"]    = e_n2
            self.eref["NO_g"]    = e_no  # eno = en + eo  ## eh2o = e_h2

class ThermoGasReference(GasReference):
    def __init__(self, ref="(A)", reaction="CO_oxidation",correct_H=True,correct_p=True):
        self.correct_H = correct_H
        self.correct_p = correct_p
        super().__init__(ref, reaction)

    def get_delta_e(self):
        super().get_delta_e()
        if self.reaction=="CO_oxidation":
            t = 600; pCO = 0.2; pO2=0.7; pCO2=0.1; pH2=0.000001; pCH4 = 0.000001; pH2O=0.000001
            self.delta_e["2CO(g)+O2(g)->CO2(g)"] += self.thermo_correct(t,p_product =[pCO2], p_educt=[pCO,pCO,pO2])
            self.delta_e["2CO(g)->2C_gr+O2(g)"]  += self.thermo_correct(t,p_product =[pO2], p_educt=[pCO,pCO])
            self.delta_e["CO(g)+3H2(g)->CH4(g)+H2O(g)"]  += self.thermo_correct(t,p_product=[pCH4,pH2O],p_educt=[pH2,pH2,pH2,pCO])
            self.delta_e["2H2(g)+O2(g)->2H2O(g)"]  += self.thermo_correct(t,p_product =[pH2O,pH2O], p_educt=[pH2,pH2,pO2])

    def thermo_correct(self,t,p_educt=[],p_product=[]):
        corr = 0
        for p in p_product:
            corr += self.correction(p,t)
        for p in p_educt:
            corr -= self.correction(p,t)
        return corr
    
    def correction(self,p,t):
        r = 8.314*10**(-3)/96.485
        corr = -t*0.002 
        if self.correct_H: corr += 7/2*r*t
        if self.correct_p: corr += r*t*np.log(p)
        return corr  
    
def set_ts_gas_pressures(reaction,gas_pressures):
     if reaction == "CO_oxidation":
          gas_species = ["CO","CO2","O2"]
          ts_gas_species = ["CO-*"]
          gas_pressures["CO-*"] = gas_pressures["CO"]
     if reaction == "NO_reduction":
          gas_species = ["H2","N2","H2O"]
          ts_gas_species = ["H2-*","N2-*"]
          gas_pressures["H2-*"] = gas_pressures["H2"]
          gas_pressures["N2-*"] = gas_pressures["N2"]
     return ts_gas_species, gas_species, gas_pressures

def set_eads(reaction,gas_ref,er,ea,thermo=False,correct_H=False,correct_p=False):
     eads = dict()
     gas_energies = GasReference(ref=gas_ref)
     if thermo:
         gas_energies = ThermoGasReference(ref=gas_ref,correct_H=correct_H,correct_p=correct_p)
     if reaction == "CO_oxidation":
          eads = dict()
          eads["CO"]    = gas_energies.eref["CO_g"]
          eads["O2"]    = gas_energies.eref["O2_g"]
          eads["CO2"]   = 0.5*gas_energies.delta_e["2CO(g)+O2(g)->CO2(g)"] + gas_energies.eref["CO_g"] + 0.5* gas_energies.eref["O2_g"]
          eads["CO*"]   = er["CO+*->CO*"] + gas_energies.eref["CO_g"]
          eads["O*"]    = 0.5*er["O2->O*+O*"] + 0.5*gas_energies.eref["O2_g"]
          eads["O-O*"]  = ea["O2->O*+O*"] + gas_energies.eref["O2_g"]
          eads["O-CO*"] = ea["O*+CO*->CO2"] + eads["CO*"] + eads["O*"]
     if reaction == "NO_reduction2":
          eads = dict()
          eads["NO"]    = gas_energies.eref["NO_g"]
          eads["H2"]    = gas_energies.eref["H2_g"]
          eads["N2"]    = gas_energies.eref["N2_g"]
          eads["H2O"]   = 0.5*gas_energies.delta_e["2H2(g)+O2(g)->2H2O(g)"]/2 + gas_energies.eref["H2"] + 0.5* gas_energies.eref["O2_g"]
          eads["CO*"]   = er["CO+*->CO*"] + gas_energies.eref["CO_g"]
          eads["O*"]    = 0.5*er["O2->O*+O*"] + 0.5*gas_energies.eref["O2_g"]
          eads["O-O*"]  = ea["O2->O*+O*"] + gas_energies.eref["O2_g"]
          eads["O-CO*"] = ea["O*+CO*->CO2"] + eads["CO*"] + eads["O*"]
     return eads

def set_path(reaction,eads):
     if reaction == "CO_oxidation":
          path = [[["CO_g","CO*","2CO*","2CO*+O2(g)","2CO*+2O*","CO*+O*+CO2_g","2CO2_g"]]]

          ts = dict()
          ts["2CO*+O2(g)"]   = "2CO*+O-O"
          ts["2CO*+2O*"]     = "CO*+O*+O-CO*"
          ts["CO*+O*+CO2_g"] = "O-CO*+CO2_g"

          e=dict()
          e["CO_g"]         = eads["CO"]
          e["CO*"]          = eads["CO*"]
          e["2CO*"]         = 2*eads["CO*"]
          e["2CO*+O2(g)"]   = 2*eads["CO*"] + eads["O2"]
          e["2CO*+O-O"]     = 2*eads["CO*"] + eads["O-O*"]
          e["2CO*+2O*"]     = 2*eads["CO*"] + 2*eads["O*"]
          e["CO*+O*+O-CO*"] = eads["CO*"] + eads["O*"] + eads["O-CO*"]
          e["CO*+O*+CO2_g"] = eads["CO*"] + eads["O*"] + eads["CO2"]
          e["O-CO*+CO2_g"]  = eads["O-CO*"] + eads["CO2"]
          e["2CO2_g"]       = 2*eads["CO2"]
     else: path = None; e=None; ts=None
     return path, e, ts

 # Sabotier Analysis functions for FT
def Rate_FT(k,ki,K,pCO,pH2,pH2O,pCH4,RDS=3):
    lambda_ch4   = 1/K[10-1] * pCH4
    lambda_h     = np.sqrt(K[2-1]*pH2)
    lambda_ch3   = 1/K[9-1] * lambda_ch4/lambda_h
    lambda_ch2   = 1/K[8-1] * lambda_ch3/lambda_h
    lambda_ch    = 1/K[7-1] * lambda_ch2/lambda_h
    lambda_c     = 1/K[6-1] * lambda_ch/lambda_h
    lambda_cch   = K[11-1] * lambda_c * lambda_ch
    lambda_c2h2  = K[12-1] * lambda_h * lambda_cch
    lambda_chch2 = K[13-1] * lambda_h * lambda_c2h2
    lambda_c2h4  = K[14-1] * lambda_h * lambda_chch2
    p_C2H4       = K[15-1] * lambda_c2h4
    if RDS==3:
        lambda_oh = 1/K[5-1] * pH2O/lambda_h
        lambda_o  = 1/K[4-1] * lambda_oh/lambda_h
        lambda_co = 10000 #K[1-1] * pCO
        #print(lambda_co)
        theta_star  = 1/(1+lambda_ch4+lambda_h+lambda_ch3+lambda_ch2+lambda_ch+lambda_c+lambda_cch+lambda_c2h2+lambda_chch2+lambda_c2h4+lambda_o+lambda_oh+lambda_co)
        theta_star2 = 1/(1+lambda_ch4+lambda_ch3+lambda_ch2+lambda_ch+lambda_c+lambda_cch+lambda_c2h2+lambda_chch2+lambda_c2h4+lambda_o+lambda_oh)
        lambda_co = 1/theta_star
        rate = k[3-1] * lambda_co*theta_star*theta_star2 - ki[3-1] * lambda_c * lambda_o * theta_star2**2
        print(rate)
        #print(k3 * lambda_co*theta_star**2,  k3i * lambda_c * lambda_o * theta_star**2)
        #print([lambda_*theta_star for lambda_ in [lambda_ch4,lambda_c2h4,lambda_h,lambda_ch3,lambda_ch2,lambda_ch,lambda_c,lambda_cch,lambda_c2h2,lambda_chch2,lambda_o,lambda_oh,lambda_co]])
    if RDS==4:
        lambda_oh = 1/K[5-1] * pH2O/lambda_h
        lambda_co = 1 #K[1-1] * pCO
        lambda_o  = K[3-1]*lambda_co/lambda_c
        theta_star  = 1/(1+lambda_ch4+lambda_h+lambda_ch3+lambda_ch2+lambda_ch+lambda_c+lambda_cch+lambda_c2h2+lambda_chch2+lambda_c2h4+lambda_o+lambda_oh+lambda_co)
        #print(lambda_o*theta_star,lambda_h*theta_star)
        theta_star2 = 1/(1+lambda_ch4+lambda_ch3+lambda_ch2+lambda_ch+lambda_c+lambda_cch+lambda_c2h2+lambda_chch2+lambda_c2h4+lambda_oh+lambda_h)
        #print(lambda_o*theta_star,lambda_h*theta_star2,theta_star2)
        rate = k[4-1] * lambda_o*lambda_h*theta_star*theta_star2 - ki[4-1] * lambda_oh*theta_star2**2
    return rate, p_C2H4
