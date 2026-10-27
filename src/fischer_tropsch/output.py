from ase.io import read
import os
import numpy as np
from ase.mep import NEBTools
from ase.thermochemistry import HarmonicThermo,IdealGasThermo
from ase.calculators.dftd3 import DFTD3
import pickle
from fischer_tropsch.catlib import SimpleMKM
from copy import deepcopy

class Output:
    def __init__(self,dname="",nebpath=[],project="Fischer-Trops",spin=True,defpath=None,nco=0,temperature=523,ef = None,nat=None, d3=False,
            gas_pressures = {'CO':6,'H2':12,'H2O':0.1,'CH4':0.0001, 'C2H4':0.0001 ,'C2H6':0.0001,'C2H2':0.0001,'CO2':0.0001,"C4H10":0.0001},
            smaller_model=False, larger_model=False, eslab = None, eslab0 = None, adv_thermo=False, not_ready=False,
            constructor=None):
        
        if constructor is not None:
            self.project       = project
            self.spin          = spin
            if self.project in ["Fischer-Trops","Steam-Reforming","CO2RR"]:
                self.adsorbate             = dname.split("_")[2]
                self.site_name             = dname.split("_")[1]
            for attr in ["dname","compound","project","facet","gas_pressures","temperature","adv_thermo","atoms","energy","eslab","freqs","thermo_ref_dname","transition_state","neb","neb_nimages","ts_energy","ts_freqs","thermo_ref_dname_ts","gas_names","gas_energies","gas_frequencies", "gas_geometries", "gas_symmetrie_numbers","gas_structures","ef_slab"]:
                if attr in constructor:
                    setattr(self,attr,constructor[attr])
                else:
                   setattr(self,attr,None)
            self.gibbs_gas()
            thermo_contr = 0
            if self.freqs is not None:
                thermo_contr = self.calculate_gibbs(self.freqs)
            self.gibbs_energy =  + thermo_contr + self.energy
            if self.adv_thermo:
                self.thermo_ref_multiplicator = constructor.get("thermo_ref_multiplicator", 1)
                if self.thermo_ref_dname is not None:
                    if isinstance(self.thermo_ref_dname, (list, tuple)):
                        ref_constructors = constructor["thermo_ref"]
                        self.thermo_ref = [Output(dname=dname, constructor=ref_constructor) for dname, ref_constructor in zip(self.thermo_ref_dname, ref_constructors)]
                        self.gibbs_energy += sum((ref.gibbs_energy-ref.energy)*mult for ref, mult in zip(self.thermo_ref, self.thermo_ref_multiplicator))
                    else:
                        self.thermo_ref_multiplicator = int(self.thermo_ref_multiplicator)
                        self.thermo_ref = Output(dname=self.thermo_ref_dname,constructor=constructor["thermo_ref"])
                        self.gibbs_energy += (self.thermo_ref.gibbs_energy-self.thermo_ref.energy)*self.thermo_ref_multiplicator
            self.eads, self.count, self.gads  = self.get_eads()
            self.ts_gads = None
            if self.transition_state is not None:
                ts_term_corr                  = self.calculate_gibbs(self.ts_freqs,ts=True)
                if ts_term_corr is None: self.ts_gibbs_energy = None
                else: self.ts_gibbs_energy          =  ts_term_corr + self.ts_energy
                if self.adv_thermo:
                    if self.thermo_ref_dname_ts is not None:
                        self.thermo_ref_ts_multiplicator = constructor.get("thermo_ref_ts_multiplicator", 1)
                        if isinstance(self.thermo_ref_dname_ts, (list, tuple)):
                            ref_constructors = constructor["thermo_ref_ts"]
                            self.thermo_ref_ts = [Output(dname=dname, constructor=ref_constructor) for dname, ref_constructor in zip(self.thermo_ref_dname_ts, ref_constructors)]
                            multiplicators = self.thermo_ref_ts_multiplicator if isinstance(self.thermo_ref_ts_multiplicator, (list, tuple)) else [self.thermo_ref_ts_multiplicator] * len(self.thermo_ref_ts)
                            if self.ts_gibbs_energy is not None:
                                self.ts_gibbs_energy += sum((ref.gibbs_energy-ref.energy)*mult for ref, mult in zip(self.thermo_ref_ts, multiplicators))
                        else:
                            self.thermo_ref_ts = Output(dname=self.thermo_ref_dname_ts,constructor=constructor["thermo_ref_ts"])
                            self.thermo_ref_ts_multiplicator = int(self.thermo_ref_ts_multiplicator)
                            if (self.thermo_ref_ts.freqs is not None or self.thermo_ref_ts.gibbs_energy is not None) and self.ts_gibbs_energy is not None:
                                self.ts_gibbs_energy += (self.thermo_ref_ts.gibbs_energy-self.thermo_ref_ts.energy)*self.thermo_ref_ts_multiplicator
                self.ts_eads, self.ts_count, self.ts_gads = self.get_eads(ts=True)
            
        else:
            self.fname                 = dname + "/vasp.out"
            if not_ready:
                self.fname             = dname + "/OSZICAR"
            self.dname                 = dname
            self.compound              = dname.split("_")[0]
            self.facet                 = dname.split("_")[1]
            self.thermo_dname          = "thermo/" + dname
            self.thermo_fname          = "thermo/" + dname + "/energy.dat"
            self.adv_thermo = adv_thermo
            self.temperature           = temperature
            self.gas_pressures         = gas_pressures
            self.d3                    = d3
            self.thermo_ref_dname      = None
            self.gibbs_energy          = None
            self.freqs                 = None
            if self.adv_thermo:
                self.thermo_ref_multiplicator = 1
                self.thermo_ref_dname = self.get_thermo_ref("details.txt")
                if self.thermo_ref_dname is not None:
                    if isinstance(self.thermo_ref_dname, (list, tuple)):
                        self.thermo_ref = [Output(ref_dname,spin=spin,temperature=temperature,gas_pressures=gas_pressures,smaller_model=smaller_model,larger_model=larger_model,eslab=eslab,adv_thermo=adv_thermo) for ref_dname in self.thermo_ref_dname]
                    else:
                        self.thermo_ref = Output(self.thermo_ref_dname,spin=spin,temperature=temperature,gas_pressures=gas_pressures,smaller_model=smaller_model,larger_model=larger_model,eslab=eslab,adv_thermo=adv_thermo)
                    self.thermo_dname = "thermo/new/" + dname
                    self.thermo_fname = "thermo/new/" + dname + "/energy.dat"
                    if self.thermo_move_atoms == []:
                        refs = self.thermo_ref if isinstance(self.thermo_ref, (list, tuple)) else [self.thermo_ref]
                        multiplicators = self.thermo_ref_multiplicator if isinstance(self.thermo_ref_multiplicator, (list, tuple)) else [self.thermo_ref_multiplicator] * len(refs)
                        if all(os.path.isfile(ref.thermo_fname) for ref in refs):
                            self.gibbs_energy = sum((ref.gibbs_energy-ref.energy)*mult for ref, mult in zip(refs, multiplicators))
            self.contcar               = dname + "/CONTCAR"
            self.project               = project
            self.atoms                 = None
            self.nco                   = nco
            self.eslab                 = eslab
            self.eslab0                = eslab0
            self.ef_slab               = ef
            self.nat                   = nat
            self.gas_names             = ["CO","H2","H2O","CH4","C2H6","C2H4","C2H2","CH3OH","CO2","C3H8","C4H10","C3H6"]
            for g in self.gas_names:
                if g not in self.gas_pressures:
                    self.gas_pressures[g] = 0.0001
            self.gas_energies          = {"CO":-14.44945908-0.18,"H2":-6.97818589+0.09,"H2O":-14.15581388-0.21,"CH4":-24.03087622,"C2H6":-40.26121986, "C2H4":-.31688127E+02, "C2H2":-.22580473E+02, "CH3OH": -.29940972E+02 ,  "CO2":-22.308091, "C3H6": -.48063156E+02, "C3H8": -56.57543631, "C4H10":-72.811436}
            self.gas_frequencies       = {"CO":[0,0,0,0,0,2101.5],"H2":[0,0,0,0,0,4311.0],"H2O":[0,0,0,0,0,0,1588.7,3691.1,3802.6],"CH4":[0, 0, 0, 0, 0, 0, 1288.8, 1289.9, 1291.1, 1509.3, 1509.5, 2961.5, 3075.1, 3082.6, 3083.4],"C2H6":[0, 0, 0, 0, 0, 0, 292.1,796.9,797.4,984.4,1177.4,1178.4,1358.8,1368.0,1452.6,1452.7,1454.9,1456.7,2965.9,2978.6,3020.4,3021.8,3045.8,3046.6], "C2H4":[0, 0, 0, 0, 0, 0, 594.4, 930.4, 941.3, 1030.9, 1119.1, 1266.5, 1432.3, 1620.0, 2914.8, 3073.7, 3116.3, 3172.0], "C2H2":[0,0,0,0,0,585.0,585.4,731.0,736.5,2001.6,3366.9,3455.7], "CH3OH":[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 306.3, 998.3, 1051.8, 1128.9, 1339.0, 1433.3, 1448.2, 1462.9, 2935.4, 2979.5, 3035.7, 3734.6] ,"CO2":[627.4,627.8,1304.6,2347.7], "C3H6": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 99.9, 406.0, 566.2, 893.4, 905.2, 948.7, 991.7, 1018.3, 1137.6, 1308.3, 1355.8, 1414.3, 1442.8, 1444.3, 1634.3, 2959.4, 2992.8, 3015.3, 3041.0, 3064.4, 3163.0], "C3H8":[31.7, 53.1, 206.9, 251.4, 342.9, 711.9, 831.9, 858.3, 895.9, 1032.4, 1100.2, 1164.5, 1260.6, 1290.6, 1320.0, 1370.3, 1415.5, 1444.5, 1451.5, 1454.1, 1467.2, 2927.3, 2939.5, 2955.3, 2960.0, 2985.4, 3026.5, 3027.5, 3063.5], "C4H10": [18.7,64.6,99.0,203.6,252.6,287.9,410.3,683.4,768.8,797.2,911.1,950.1,965.5,1037.7,1108.6,1153.3,1246.7,1263.1,1306.3,1327.1,1355.8,1368.8,1439.4,1441.6,1458.3,1460.4,1464.6,1526.4,2720.0,2930.8,2941.4,2948.6,2960.5,2981.1,2994.2,3008.8,3029.5,3353.7]}
            self.gas_structures        = {"CO":read("/home/cat/bengr/Fischer_Trops/gas_contcar/CO.CONTCAR"),"H2":read("/home/cat/bengr/Fischer_Trops/gas_contcar/H2.CONTCAR"),"H2O":read("/home/cat/bengr/Fischer_Trops/gas_contcar/H2O.CONTCAR"),
                                        "CH4":read("/home/cat/bengr/Fischer_Trops/gas_contcar/CH4.CONTCAR"),"C2H6":read("/home/cat/bengr/Fischer_Trops/gas_contcar/C2H6.CONTCAR"),"C2H4":read("/home/cat/bengr/Fischer_Trops/gas_contcar/C2H4.CONTCAR"),
                                        "C2H2":read("/home/cat/bengr/Fischer_Trops/gas_contcar/C2H2.CONTCAR"),"CH3OH":read("/home/cat/bengr/Fischer_Trops/gas_contcar/CH3OH.CONTCAR") ,"CO2":read("/home/cat/bengr/Fischer_Trops/gas_contcar/CO2.CONTCAR"), "C3H6":read("/home/cat/bengr/Fischer_Trops/gas_contcar/C3H6.CONTCAR") ,"C3H8":read("/home/cat/bengr/Fischer_Trops/gas_contcar/C3H8.CONTCAR"),"C4H10":read("/home/cat/bengr/Fischer_Trops/gas_contcar/C4H10.CONTCAR")}
            self.gas_free_energies     = dict()
            self.gas_geometries        = {"CO":"linear","H2":"linear","H2O":"nonlinear","CH4":"nonlinear", "C2H4":"nonlinear" , "C2H6":"nonlinear", "C2H2":"linear", "CH3OH":"nonlinear" ,"CO2":"linear","C3H6":"nonlinear","C3H8":"nonlinear","C4H10":"nonlinear"}
            self.gas_symmetrie_numbers = {"CO":1,"H2":2,"H2O":2,"CH4":4,"C2H4":2,"C2H6":3,"C2H2":1, "CH3OH":1 ,"CO2":2, "C3H6":1 , "C3H8":2, "C4H10":1}
            
            if self.d3:
                #self.gas_energies          = {"CO":-14.44945908,"H2":-6.97818589,"H2O":-14.15581388,"CH4":-24.03087622,"C2H6":-40.26121986, "C2H4":-.31688127E+02, "C2H2":-.22580473E+02}
                for g in self.gas_names:
                    self.gas_energies[g] += self.dftd3(self.gas_structures[g])

            if spin != None:
                self.spin = spin
            else:
                if self.compound == "Ni":
                    self.spin = False
                else:
                    self.spin = True
            self.smaller_model = smaller_model
            self.larger_model  = larger_model

            self.energy                = self.get_energy()
            if "slab" in self.dname:
                self.gibbs_energy = self.energy
                self.eads = self.energy - eslab
                self.gads = self.eads
                return
            if self.gibbs_energy is not None: self.gibbs_energy += self.energy
            self.gibbs_gas()
            if os.path.isfile(self.thermo_fname):
                self.freqs        = self.get_freqs()
                self.gibbs_energy = self.calculate_gibbs(self.freqs) + self.energy
                if self.adv_thermo and self.thermo_ref_dname is not None:
                    refs = self.thermo_ref if isinstance(self.thermo_ref, (list, tuple)) else [self.thermo_ref]
                    multiplicators = self.thermo_ref_multiplicator if isinstance(self.thermo_ref_multiplicator, (list, tuple)) else [self.thermo_ref_multiplicator] * len(refs)
                    if all(ref.freqs is not None or ref.gibbs_energy is not None for ref in refs):
                        self.new_gibbs_energy = self.calculate_gibbs(self.freqs)
                        self.gibbs_energy = self.new_gibbs_energy + sum((ref.gibbs_energy-ref.energy)*mult for ref, mult in zip(refs, multiplicators)) + self.energy

            if os.path.isfile(dname + "/final.traj"):
                self.atoms  = read(dname + "/final.traj")
                self.energy = self.atoms.get_potential_energy()
            else:
                self.atoms = read(dname + "/CONTCAR")

            if self.d3:
                d3_energy = self.dftd3(self.atoms)
                self.energy += d3_energy
                if self.gibbs_energy is not None:
                    self.gibbs_energy += d3_energy

            if self.project in ["Fischer-Trops","CO2RR","Steam-Reforming"]:
                self.adsorbate             = dname.split("_")[2]
                self.site_name             = dname.split("_")[1]
            elif self.project=="Alloys":
                self.adsorbate             = dname.split("_")[1]
                self.site_name             = "bcc111"
            self.eads, self.count, self.gads      = self.get_eads()

            self.neb                   = None; self.transition_state = None; self.ts_energy=None; self.ts_eads=None; self.ts_count=None;self.barrier=None;self.barrier_reverse=None; self.ts_thermo_fname = None; self.ts_gads =None
            self.nebpath               = None; self.ts_freqs = None
            if dname in nebpath:
                self.nebpath = nebpath[dname][0]
            if self.nebpath!=None and ((os.path.isfile(self.nebpath + "/neb.traj") and os.stat(self.nebpath + "/neb.traj").st_size != 0) or (os.path.isfile(self.nebpath + "/predicted.traj") and os.stat(self.nebpath + "/predicted.traj").st_size != 0)):
                nebtraj = "neb.traj"
                if nebpath[dname][1]=="auto":
                    self.neb_nimages = get_neb_nimages(self.nebpath+"/neb.traj")
                    self.ts_image_nr = self.get_ts_image_number(self.neb_nimages,self.nebpath+"/neb.traj")
                elif nebpath[dname][1]=="mlneb":
                    self.neb_nimages = get_neb_nimages(self.nebpath+"/predicted.traj")
                    if len(nebpath[dname])>2:
                        self.neb_nimages = nebpath[dname][2]
                    self.ts_image_nr = self.get_ts_image_number(self.neb_nimages,self.nebpath+"/predicted.traj")
                    nebtraj = "predicted.traj"
                else:
                    self.neb_nimages  = nebpath[dname][1]
                    self.ts_image_nr = nebpath[dname][2]
                if len(nebpath[dname])>3:
                    cycle = nebpath[dname][3]
                    cont=True
                    if nebpath[dname][1]=="mlneb":
                        cont=False
                        if len(nebpath[dname])>4:
                            cycle = nebpath[dname][4]
                            cont = True
                    if cont:
                         self.neb = []
                         for i in range(0,self.neb_nimages+2):
                             a = (cycle-1)*(self.neb_nimages+2) + i
                             self.neb.append(read(self.nebpath + '/' + nebtraj  + '@'+str(a)))
                    else:
                        self.neb = read(self.nebpath + "/" + nebtraj + "@-"+str(self.neb_nimages+2)+":")
                else:
                    self.neb = read(self.nebpath + "/" + nebtraj + "@-"+str(self.neb_nimages+2)+":")

                self.thermo_dname_ts = "thermo/" + dname + "_ts"
                self.ts_gibbs_energy = None; self.barrier_g = None; self.barrier_reverse_g = None; self.ts_freqs=None

                if self.adv_thermo:
                    self.thermo_ref_ts_multiplicator = 1
                    self.thermo_ref_dname_ts = self.get_thermo_ref("details.txt",ts=True)
                    if self.thermo_ref_dname_ts is not None:
                        if isinstance(self.thermo_ref_dname_ts, (list, tuple)):
                            self.thermo_ref_ts = [Output(ref_dname,spin=spin,temperature=temperature,gas_pressures=gas_pressures,smaller_model=smaller_model,larger_model=larger_model,eslab=eslab,adv_thermo=adv_thermo) for ref_dname in self.thermo_ref_dname_ts]
                        else:
                            self.thermo_ref_ts = Output(self.thermo_ref_dname_ts,spin=spin,temperature=temperature,gas_pressures=gas_pressures,smaller_model=smaller_model,larger_model=larger_model,eslab=eslab,adv_thermo=adv_thermo)
                        self.thermo_dname_ts = "thermo/new/" + dname + "_ts"
                        if self.thermo_move_atoms_ts == []:
                            refs = self.thermo_ref_ts if isinstance(self.thermo_ref_ts, (list, tuple)) else [self.thermo_ref_ts]
                            multiplicators = self.thermo_ref_ts_multiplicator if isinstance(self.thermo_ref_ts_multiplicator, (list, tuple)) else [self.thermo_ref_ts_multiplicator] * len(refs)
                            if all(os.path.isfile(ref.thermo_fname) for ref in refs):
                                self.ts_gibbs_energy = sum((ref.gibbs_energy-ref.energy)*mult for ref, mult in zip(refs, multiplicators))

                self.thermo_fname_ts        = self.thermo_dname_ts + "/energy.dat"
                self.transition_state       = self.neb[self.ts_image_nr]
                self.ts_energy              = self.transition_state.get_potential_energy()
                
                if self.d3:
                    d3_energy = self.dftd3(self.transition_state)
                    self.ts_energy += d3_energy

                if self.ts_gibbs_energy is not None: self.ts_gibbs_energy += self.ts_energy
                
                if os.path.isfile(self.thermo_fname_ts):
                    self.ts_freqs                 = self.get_freqs(ts=True)
                    self.ts_gibbs_energy          = self.calculate_gibbs(self.ts_freqs,ts=True) + self.ts_energy
                    if self.adv_thermo:
                        if self.thermo_ref_dname_ts is not None:
                            refs = self.thermo_ref_ts if isinstance(self.thermo_ref_ts, (list, tuple)) else [self.thermo_ref_ts]
                            multiplicators = self.thermo_ref_ts_multiplicator if isinstance(self.thermo_ref_ts_multiplicator, (list, tuple)) else [self.thermo_ref_ts_multiplicator] * len(refs)
                            if all(ref.freqs is not None or ref.gibbs_energy is not None for ref in refs):
                                self.ts_gibbs_energy += sum((ref.gibbs_energy-ref.energy)*mult for ref, mult in zip(refs, multiplicators))

                self.ts_eads, self.ts_count, self.ts_gads = self.get_eads(ts=True)
                self.barrier                = self.ts_energy - self.energy #self.neb[0].get_potential_energy()
                if os.path.isfile(self.thermo_fname_ts):
                    if self.gibbs_energy!=None: self.barrier_g              = self.ts_gibbs_energy - self.gibbs_energy #self.neb[0].get_potential_energy()
                    else: self.barrier_g = None
                if self.smaller_model:# and not self.spin and self.compound=="Ni":
                    self.neb[0] = self.atoms
            self.dbe = None
            if defpath!=None:
                self.defpath = defpath
                self.dbe = self.get_dbe()

    def get_energy(self,fname=None):
        if fname is None:
            fname = self.fname
        f=open(fname)
        l=f.readlines()
        e=0
        for line in l:
            spl = line.split()
            if len(spl) >= 5:
                if spl[3]=="E0=":
                    e = float(spl[4])
        return e

    def dftd3(self,atoms,functional="RPBE"):
        calc = DFTD3(damping='bj', xc=functional)
        atoms.calc = calc
        d3_energy = atoms.get_potential_energy()
        return d3_energy

    def get_dbe(self):
        e_adsorbate = self.defpath[1]
        fname = self.defpath[0] + "/" + self.dname + "_defect/vasp.out"
        if os.path.isfile(fname):
            e_def = self.get_energy(fname)
            e_dbe = e_def + e_adsorbate - self.energy
            return round(e_dbe,2)
        else: return None

    def get_g_minus_e(self,ts=False):
        if ts: f= open(self.ts_thermo_fname)
        else: f= open(self.thermo_fname)
        l = f.readlines()
        g = float(l[1].split()[-2])
        e = float(l[0].split()[-2])
        return g - e

    def get_eref(self,ts=False):
        if ts: f= open(self.ts_thermo_ref_fname)
        else: f= open(self.thermo_ref_fname)
        l = f.readlines()
        e = float(l[0].split()[-2])
        return e

    def get_thermo_ref(self,fname,ts=False):
        dname = self.dname
        if ts: dname+="_ts"
        f = open(fname)
        l = f.readlines()
        for line in l:
            spl = line.split()
            if dname.endswith("+"+spl[0]) or dname.endswith("_"+spl[0]):
                if spl[1]=="[]": numbers=[]
                else: numbers = list(spl[1].split(","))
                if not ts: self.thermo_move_atoms = [int(float(n)) for n in numbers]
                else: self.thermo_move_atoms_ts = [int(float(n)) for n in numbers]

                # >>> CHANGED: [] may define several references and multiplicators
                # New format:
                #   ADS [] [ref1,ref2] [mult1,mult2]
                # Old single-reference format remains unchanged.
                if spl[1] == "[]" and spl[2].startswith("[") and spl[2].endswith("]"):
                    refs = [ref.strip() for ref in spl[2][1:-1].split(",") if ref.strip()]
                    if len(spl) > 3 and spl[3].startswith("[") and spl[3].endswith("]"):
                        multiplicators = [int(float(mult.strip())) for mult in spl[3][1:-1].split(",") if mult.strip()]
                    else:
                        multiplicators = [1] * len(refs)
                    if len(refs) != len(multiplicators):
                        raise ValueError("Number of thermo references and multiplicators differs for " + dname)
                    if not ts:
                        self.thermo_ref_multiplicator = multiplicators
                    else:
                        self.thermo_ref_ts_multiplicator = multiplicators
                    return [self.compound + "_" + self.facet + "_" + ref for ref in refs]
                # <<< CHANGED

                if len(spl) > 3:
                    if not ts:
                        self.thermo_ref_multiplicator = int(float(spl[3]))
                    else:
                        self.thermo_ref_ts_multiplicator = int(float(spl[3]))
                return self.compound + "_" + self.facet + "_" + spl[2]
        return None

    def get_eads(self,ts=False):
        e_slab = self.eslab
        if self.project is None: self.project = "Fischer-Trops"

        if self.project in ["Fischer-Trops", "CO2RR","Steam-Reforming"]:
            e_h2o  = self.gas_energies["H2O"]; g_h2o  = self.gas_free_energies["H2O"]
            e_h2   = self.gas_energies["H2"];  g_h2  = self.gas_free_energies["H2"] 
            e_co   = self.gas_energies["CO"];  g_co  = self.gas_free_energies["CO"]  
            e_ch4   = self.gas_energies["CH4"];  g_ch4  = self.gas_free_energies["CH4"]  
            e_c2h6   = self.gas_energies["C2H6"];  g_c2h6  = self.gas_free_energies["C2H6"]  
        elif self.project == "Alloys":
            e_h2   = -6.97818726      ; g_h2  = -7.341075962009784
            e_n2   = -16.23720107     ; g_n2  = -17.4818241779517
            e_nh3  =  -19.46522364    ; g_nh3 = -19.9305653774276

        el_count = self.ads_count()
        count = dict()
        if self.project in ["Fischer-Trops","CO2RR"]:
            count["CO"]  = el_count["C"]
            count["H2O"] = el_count["O"] - el_count["C"]
            count["H2"]  = 1/2*(el_count["H"]) - count["H2O"]
        
        elif self.project == "Steam-Reforming":
            count["CH4"]  = el_count["C"]
            count["H2O"]  = el_count["O"]
            count["H2"]   = 1/2*(el_count["H"] - 4* count["CH4"] - 2*count["H2O"])
        
        elif self.project == "Alloys":
            count["N2"] = 0
            count["NH3"] = el_count["N"]
            count["H2"]  = (el_count["H"] - 3*el_count["N"])/2
            if count["NH3"] == 0 and count["H2"]==0:
                count["N2"]  = -1; count["H2"]=-3; count["NH3"]=2
        
        if self.project in ["Fischer-Trops","CO2RR"]:
            e_ads = self.energy - e_slab - count["CO"]*e_co - count["H2O"]*e_h2o - count["H2"]*e_h2
            if ts: e_ads = self.ts_energy - e_slab - count["CO"]*e_co - count["H2O"]*e_h2o - count["H2"]*e_h2

            if (ts and self.ts_gibbs_energy == None) or  (not ts and self.gibbs_energy == None):
                return round(e_ads,2), count, None
            
            if ts: g_ads = self.ts_gibbs_energy - e_slab - count["CO"]*g_co - count["H2O"]*g_h2o - count["H2"]*g_h2
            else: g_ads = self.gibbs_energy - e_slab - count["CO"]*g_co - count["H2O"]*g_h2o - count["H2"]*g_h2

        elif self.project =="Steam-Reforming":
            e_ads = self.energy - e_slab - count["CH4"]*e_ch4 - count["H2O"]*e_h2o - count["H2"]*e_h2
            if ts: e_ads = self.ts_energy - e_slab - count["CH4"]*e_ch4 - count["H2O"]*e_h2o - count["H2"]*e_h2

            if (ts and self.ts_gibbs_energy == None) or  (not ts and self.gibbs_energy == None):
                return round(e_ads,2), count, None
            
            if ts: g_ads = self.ts_gibbs_energy - e_slab - count["CH4"]*g_ch4 - count["H2O"]*g_h2o - count["H2"]*g_h2
            else: g_ads = self.gibbs_energy - e_slab - count["CH4"]*g_ch4 - count["H2O"]*g_h2o - count["H2"]*g_h2

        if self.project =="Alloys":
            e_ads = self.energy - e_slab - count["NH3"]*e_nh3 - count["H2"]*e_h2 - count["N2"]*e_n2
            if ts: e_ads = self.ts_energy - e_slab - count["NH3"]*e_nh3 - count["H2"]*e_h2 - count["N2"]*e_n2

            if (ts and self.ts_gibbs_energy == None) or  (not ts and self.gibbs_energy == None):
                if not self.dname.endswith("slab"):
                    return round(e_ads,2), count, None
                else:
                    g_ads = self.energy - e_slab - count["NH3"]*g_nh3 - count["H2"]*g_h2 - count["N2"]*g_n2
                    return round(e_ads,2), count, round(g_ads,2)

            g_ads = self.gibbs_energy - e_slab - count["NH3"]*g_nh3 - count["H2"]*g_h2 - count["N2"]*g_n2
            if ts:
                g_ads = self.ts_gibbs_energy - e_slab - count["NH3"]*g_nh3 - count["H2"]*g_h2 - count["N2"]*g_n2

        return round(e_ads,2), count, round(g_ads,2)

    def ads_count(self):
        i=0
        count = dict()
        count["C"] = 0
        count["H"] = 0
        count["O"] = 0
        count["N"] = 0

        if hasattr(self,"smaller_model") and self.smaller_model:
            if self.compound == "Ni":
                if not self.larger_model:
                    count["C"]-=3
                    count["O"]-=3
                    count["H"]-=4
                else:
                    count["C"]-=7
                    count["O"]-=7
            if self.compound in ["Co"]:
                if self.facet in ["hcp10m15","hcp10m11"]: count["C"]-=7;count["O"]-=7
                elif self.facet == "hcp11m21": count["C"]-=4;count["O"]-=4
            if self.compound in ["Ru"]:
                if self.facet == "hcp10m15": count["C"]-=7;count["O"]-=7
                elif self.facet == "hcp11m21": count["C"]-=3;count["O"]-=3

        if hasattr(self,"nco") and self.nco>0: count["C"]-= self.nco; count["O"]-=self.nco

        if self.adsorbate =="slab": return count
        for c in self.adsorbate:
            if c in ["1","2","3","4","5","6","7","8","9","0"]:
                do=True
                if len(self.adsorbate) <= i+2: do=False

                 # Multi-atomic adsorbates with multi counting (at the moment only CO (n_CO<100) and NH3 for (n_NH3<9) supported)

                if do and self.adsorbate[i-1] == "+" and self.adsorbate[i+1] == "C" and self.adsorbate[i+2] == "O":  # e.g. +2CO
                    count["C"] += int(c) - 1
                    count["O"] += int(c) - 1
                elif do and self.adsorbate[i-1] == "+" and self.adsorbate[i+2] == "C" and self.adsorbate[i+3] == "O": # e.g. +20CO
                    count["C"] += 10*int(c) 
                    count["O"] += 10*int(c)
                elif do and self.adsorbate[i-2] == "+" and self.adsorbate[i+1] == "C" and self.adsorbate[i+2] == "O": # e.g. +24CO (count also the 4)
                    count["C"] += int(c) - 1
                    count["O"] += int(c) - 1

                elif do and self.adsorbate[i-1] == "+" and self.adsorbate[i+1] == "N" and self.adsorbate[i+2] == "H" and self.adsorbate[i+2] == "3":
                    count["N"] += int(c) - 1
                    count["H"] += 3*(int(c) - 1)
                
                # One-atomic adsorbates with multi counting

                elif self.adsorbate[i-1] == "+":
                    for ads in ["C","H","O"]:
                        if self.adsorbate[i+1] == ads:   # e.g. +2H
                            count[ads] += int(c) - 1
                        elif i+2 < len(self.adsorbate) and self.adsorbate[i+2] == ads: # e.g. +20 H
                            count[ads] += 10*int(c)
                elif self.adsorbate[i-2] == "+":
                    found=False
                    for ads in ["C","H","O"]:
                        if self.adsorbate[i-1] in ["1","2","3","4","5","6","7","8","9","0"] and i+1 < len(self.adsorbate):
                            if self.adsorbate[i+1] == ads:
                                count[ads] += int(c) -1   # e.g. +24 H (count also the 4)
                                found=True
                    if not found:
                        count[self.adsorbate[i-1]]+=int(c)-1

                # Any other adsorbates without multi-counting

                else:
                    count[self.adsorbate[i-1]]+=int(c)-1  

            elif c=="+" or c=="-":
                i+=1;continue
            else:
                count[c]+=1
            i+=1
        return count

    def auto_detect_neb_dat(self,ase_skript,traj):
        f = open(ase_skript);l=f.readlines()
        n_images = 0
        for line in l:
            spl = line.split()
            if len(spl)<3: continue
            if spl[0] == "n_images" and spl[1] == "=":
                n_images = spl[2]
        images = read(traj+"@-"+str(int(n_images)+2)+":")
        e = []
        for i in range(1,len(images)-1):
            e.append([images[i].get_potential_energy(),i])
        e = sorted(e, key=lambda x: x[0], reverse=True)
        return int(n_images), int(e[0][1])

    def get_ts_image_number(self,n_images,traj):
        images = read(traj+"@-"+str(int(n_images)+2)+":")
        e = []
        for i in range(1,len(images)-1):
            e.append([images[i].get_potential_energy(),i])
        e = sorted(e, key=lambda x: x[0], reverse=True)
        return int(e[0][1])

    def nebplot(self):
        images = self.neb
        nebtools = NEBTools(images=images)
        nebtools.plot_bands(len(images))
        fname = self.nebpath.replace("/","_") + ".pdf"
        os.system("mv nebplots.pdf " + fname)

    def nebplot(self):
        images = self.neb
        nebtools = NEBTools(images=images)
        nebtools.plot_bands(len(images))
        fname = self.nebpath.replace("/","_") + ".pdf"
        os.system("mv nebplots.pdf " + fname)

    def read_freq_lines(self,filename):
        with open(filename, "r") as f:
            lines = f.readlines()
        zpe_indices = [
            i for i, line in enumerate(lines)
            if line.startswith("Zero-point energy:")
        ]

        if len(zpe_indices) <= 1:
            return lines
        else:
            last_idx = zpe_indices[-2]
            return lines[last_idx + 1:]

    def get_freqs(self,fname="",ts=False):
        if fname == "":
            fname = self.thermo_fname
            if ts:
                fname = self.thermo_fname_ts
        freqs = []
        fname = fname.replace("energy.dat","vib.dat")
        l = self.read_freq_lines(fname)
        for line in l:
            spl = line.split()
            if spl[-1] !="cm^-1" and len(spl) == 3:
                if "i" in spl[-1]: freqs.append("0")
                else: freqs.append(spl[-1])
        return freqs
    
    def gibbs_gas(self,update=False):
        for specie in self.gas_names:
            if  hasattr(self, "gas_free_energies") and specie in self.gas_free_energies and not update:
                continue

            atoms = self.gas_structures[specie]
            freqs = self.gas_frequencies[specie]
            mass = sum(atoms.get_masses())
            
            #------Handle frequencies, geometry------
            geom = self.gas_geometries[specie]
            if geom == 'nonlinear':
                nb_freqs = 3 * len(atoms) - 6

            elif geom == 'linear':
                nb_freqs = 3 * len(atoms) - 5

            elif geom == 'monoatomic':
                nb_freqs = 0
                
            else:
                raise Exception('Incorrect geometry for ' + specie)
                
            if len(freqs) < nb_freqs:
                raise Exception('Too few frequencies for ' + specie)
            
            freqs = freqs[-nb_freqs:] if nb_freqs > 0 else [] #Check
            
            #-------------Apply cutoff-------------
            for i,freq in enumerate(freqs):
                if isinstance(freq,complex) or freq < 50.0:
                    print('Warning! Low frequency for gas ' + specie + '. ' + str(freq) + \
                            ' replaced by cutoff frequency')
                    freqs[i] = 50.0
                else:
                    #Ascending values
                    break
            
            #-----------------------------------------
            invcm_in_eV = 1.23981e-4

            thermo = IdealGasThermo(vib_energies=[freq * invcm_in_eV for freq in freqs],
                                    potentialenergy=self.gas_energies[specie],
                                    geometry=geom,
                                    symmetrynumber=self.gas_symmetrie_numbers[specie],
                                    spin=0,atoms=atoms)

            pressure = 10**5 * self.gas_pressures[specie]
            if not hasattr(self,"gas_free_energies"): self.gas_free_energies = dict()
            self.gas_free_energies[specie] = thermo.get_gibbs_energy(self.temperature,pressure,verbose=False)

    def calculate_gibbs(self,freqs,ts=False):
        if freqs is None: return None
        freqs = [float(f) for f in freqs]
        if ts:
            if len(freqs) % 3 == 0 and (isinstance(freqs[0],complex) or freqs[0]<50.0):
                # Three freqs per atom, one two much for a transition state
                freqs = freqs[1:]
        else:
            if not len(freqs) % 3 == 0:
                raise Exception('Number of frequencies is not a multiple of 3')
                
        for i,freq in enumerate(freqs):
            if isinstance(freq,complex) or freq < 50.0:
                freqs[i] = 50.0
            else:
                # As values are ascending, we're done
                break

        invcm_in_eV = 1.23981e-4
        thermo = HarmonicThermo(vib_energies=[freq * invcm_in_eV for freq in freqs],
                                potentialenergy=0)
        if ts:
            if self.freqs is not None and len(self.freqs) == len(freqs):
                if max([(freqs[i]- float(self.freqs[i])) for i in range(len(freqs))]) == 0:
                    return self.gibbs_energy - self.energy
            return thermo.get_helmholtz_energy(self.temperature,verbose=False)
        else:
            return thermo.get_helmholtz_energy(self.temperature,verbose=False)

import numpy as np
from ase.io import iread

def get_neb_nimages(predicted="predicted.traj", tol=1e-6):
    images = list(iread(predicted, index=":"))
    if len(images) < 2:
        raise RuntimeError("Not enough frames.")

    first_pos = images[0].positions
    first_nums = images[0].numbers

    hits = []
    for i, atoms in enumerate(images[1:], start=1):
        if len(atoms) != len(images[0]):
            continue
        if not np.array_equal(atoms.numbers, first_nums):
            continue
        if np.max(np.abs(atoms.positions - first_pos)) < tol:
            hits.append(i)

    if not hits:
        n_img = len(images)
    else:
        n_img = hits[0]

    if len(images) % n_img != 0:
        raise RuntimeError(
            f"Detected {n_img} images/path, but total frames {len(images)} "
            "is not divisible by that."
        )

    return n_img -2


def catmap_input(outputs=None):
    if outputs==None:
        with open("structures.pkl","rb") as fin:
            outputs = pickle.load(fin)
    f=open("energy_inputs.txt","w")
    f.write("\t".join(["surface_name","site_name","species_name","formation_energy","frequencies","reference"])+"\n")
    for out in outputs:
        if out.freqs==None: continue
        freqs = out.freqs
        freq_string = "[" + ",".join([str(f) for f in freqs]) + "]"
        f.write("\t".join([out.compound,out.site_name,out.adsorbate,str(out.eads),freq_string,"This work"]))
        if out!= outputs[-1]: f.write("\n")
        if out.nebpath!=None:
            freqs = out.ts_freqs
            if freqs is None: continue
            freq_string = "[" + ",".join([str(f) for f in freqs]) + "]"
            f.write("\t".join([out.compound,out.site_name,out.adsorbate+"_ts",str(out.ts_eads),freq_string,"This work"]))
            f.write("\n")
    return outputs

import ast

def add_reference(obj, constructors, ts=False):
    name = "thermo_ref_dname"
    ref  = "thermo_ref"
    if ts: name+="_ts"; ref+="_ts"
    if name in constructors[obj]:
        ref_target = constructors[obj][name]
        constructors[obj][ref] = add_reference(ref_target, constructors)
        
    return constructors[obj]

def read_outputs(path,gas_path):
    files = os.listdir(path)
    outputs = []
    constructors = {}
    for f in files:
        if f.endswith(".txt"):
            dname = f.replace(".txt","")
            constructors[dname] = dict()
            file = open(os.path.join(path,f))
            l = file.readlines()
            for line in l:
                line = line.replace("\n","")
                if not " " in line: continue
                attr, value = line.split(" ",1)
                try:
                    constructors[dname][attr] = ast.literal_eval(value.strip())
                except (ValueError, SyntaxError):
                    constructors[dname][attr] = value.strip()

            constructors[dname]["atoms"] = read(os.path.join(path,dname + ".traj"))
            if os.path.isfile(os.path.join(path,dname) + "_ts.traj"):
                constructors[dname]["transition_state"] = read(os.path.join(path,dname + "_ts.traj"))
            if os.path.isfile(os.path.join(path,dname) +"_neb.traj"):
                constructors[dname]["neb"] = read(os.path.join(path,dname) + "_neb.traj"+"@-"+str(constructors[dname]["neb_nimages"])+":")
    for c in constructors:
        constructors[c] = add_reference(c,constructors)
        if "transition_state" in constructors[c]:
            constructors[c] = add_reference(c,constructors,ts=True)
        constructors[c]["gas_structures"] = dict()
        for g in constructors[c]["gas_names"]:
            constructors[c]["gas_structures"][g] = read(os.path.join(gas_path,g+".CONTCAR"))

    outputs = []
    for c in constructors:
        outputs.append(Output(dname=c,constructor=constructors[c]))
    return outputs

import re
import copy

def change_comps(comp,comp_new,comp_old):
    i=0
    for c in comp_new:
        comp[c] = copy.deepcopy(comp[comp_old[i]])
        comp[c].nickname = c
        i+=1
    return comp

TERM_RE = re.compile(
    r"""
    (\d*)                                   # optionaler Koeffizient
    (
        [A-Za-z][A-Za-z0-9_]*(?:\(g\))?     # erster Block (optional mit (g))
        (?:-[A-Za-z][A-Za-z0-9_]*(?:\(g\))?)*  # optionale -Blöcke (TS), ebenfalls mit (g)
    )
    """,
    re.VERBOSE
)

def free_energy_of_combo(combo, comps, free_energies = None,gn=None):
    """
    combo: z.B. 'A2B-CD4+2CO2+A-B'
    comps: dict mit comps[species].free_energy
    """
    total = 50.0

    for term in combo.split("+"):
        if "(g)" in term:
            if gn is None: continue
        match = TERM_RE.fullmatch(term)
        if not match:
            raise ValueError(f"Ungültiger Term: '{term}'")

        n_str, species = match.groups()
        n = int(n_str) if n_str else 1

        if species in comps:
            if free_energies is not None:
                free_energies[species] = comps[species].free_energy
            if comps[species].free_energy is not None:
                if total==50: total=0
                total += n * comps[species].free_energy
        elif gn is not None:
            if total==50: total=0
            if "(g)" in species:
                total += n*gn[species.replace("(g)","")]

    if free_energies is not None:
        free_energies[combo] = total
        return free_energies
    
    return total

class Compound:
    def __init__(self,adsorbate,ts,comp_ref,nickname,outputs,gas_e,gas_g):
        self.adsorbate = adsorbate
        self.metal = outputs[0].compound
        self.spin  = outputs[0].spin
        self.facet = outputs[0].facet
        self.comp_ref  = comp_ref
        self.nickname = nickname
        self.gas_e = gas_e
        self.gas_g = gas_g
        self.ts = ts
        self.atoms     = None
        self.atoms_ref = None
        self.neb       = None
         
        self.init_energy      = self.get_g(outputs,self.adsorbate,ts=self.ts,type="e")
        self.energy = None
        self.free_energy = None

        self.init_free_energy   = self.get_g(outputs,self.adsorbate,ts=self.ts,type="g")
        self.init_correction    = self.get_g(outputs,self.adsorbate,ts=self.ts,type="g-e")
        if self.comp_ref is not None and self.comp_ref != "None":
            self.ref_energy       = self.get_g(outputs,self.comp_ref,ts=False,type="e")
            self.energy = self.init_energy - self.ref_energy
            self.ref_free_energy  = self.get_g(outputs,self.comp_ref,ts=False,type="g")
            self.ref_correction   = self.get_g(outputs,self.comp_ref,ts=False,type="g-e")
            if self.init_free_energy is None:
                self.free_energy = None
                self.correction = None
            else:
                if self.ref_free_energy is None:
                    self.free_energy is None
                else:
                    self.free_energy   = self.init_free_energy - self.ref_free_energy
                    self.correction    = self.init_correction - self.ref_correction
            self.gf = None
            try:
                if self.ef_slab is not None:
                    self.gf = self.get_gf()
            except AttributeError:
                return
        else:
            self.energy = self.init_energy - outputs[0].eslab
            #print(self.energy,self.adsorbate,self.init_energy)
            if self.init_free_energy is None:
                self.free_energy = None
                self.correction = None
            else:
                self.free_energy = self.init_free_energy -  outputs[0].eslab
                self.correction_energy = self.init_correction
            if hasattr(outputs[0],"ef_slab") and outputs[0].ef_slab is not None:
                self.gf = self.get_gf()

    def get_g(self,outputs,ads,ts=False,type="g",add=0):
        for out in outputs:
            if out.dname.endswith("_"+ads):
                count = out.ads_count()
                if hasattr(out,"ef_slab"):
                    self.ef_slab = out.ef_slab
                else: self.ef_slab = None
                if ads == self.adsorbate and self.atoms is None and not ts:
                    self.atoms     = out.atoms
                if ads == self.adsorbate and self.atoms is None and ts:
                    self.neb     = out.neb
                    self.atoms   = out.transition_state
                if ads == self.comp_ref and self.atoms is None:
                    self.atoms_ref = out.atoms

                if type == "e" or type =="g-e":
                    if not ts:
                        e = out.energy - count["O"]*self.gas_e["O"] - count["H"]*self.gas_e["H"] - count["C"]*self.gas_e["C"]  + add
                    else:
                        e =  out.ts_energy - count["O"]*self.gas_e["O"] - count["H"]*self.gas_e["H"] - count["C"]*self.gas_e["C"]  + add
                    if type =="e": return e
                if type == "g" or type =="g-e":
                    if not ts:
                        if out.gibbs_energy is None: return None
                        g = out.gibbs_energy - count["O"]*self.gas_g["O"] - count["H"]*self.gas_g["H"] - count["C"]*self.gas_g["C"]  + add
                    else:
                        if out.ts_gibbs_energy is None: return None
                        g = out.ts_gibbs_energy - count["O"]*self.gas_g["O"] - count["H"]*self.gas_g["H"] - count["C"]*self.gas_g["C"]  + add
                    if type =="g": return g
                return g-e + add
        #if ads != "None":
        #    print(ads,self.metal,self.facet,self.spin,self.adsorbate)
        return 0
    
    def get_gf(self):
        if self.free_energy is None: return self.ef_slab + self.energy
        if not hasattr(self,"ef_slab"): return 0
        gf = self.ef_slab + self.free_energy
        return gf
    
def read_comps(t,p,op,path):
    ge,gfe = get_gas_data(t,p,op)
    en,gn,ge,gfe = get_relative_gas_data(ge,gfe)
    metal = os.path.basename(os.path.normpath(path))
    f = open(os.path.join(path,metal+".compdata"))
    lines = f.readlines()
    comps = dict()
    for l in lines:
        spl = l.split()
        if len(spl)!=4: continue
        nickname  = spl[0]
        adsorbate = spl[1]
        ts        = (spl[2]=="True")
        comp_ref  = spl[3]
        if comp_ref=="None": comp_ref = None
        comps[nickname]= Compound(adsorbate,ts,comp_ref,nickname,op,ge,gfe)
    return comps

def read_data_structure_class(path):
    metal = os.path.basename(os.path.normpath(path))
    f = open(os.path.join(path,metal+".dscdata"))
    lines = f.readlines()
    constructor = dict()
    for l in lines:
        spl = l.split(" ",1)
        if len(spl)!=2: continue
        try:
            constructor[spl[0]] = ast.literal_eval(spl[1].strip())
        except (ValueError, SyntaxError):
            constructor[spl[0]] = spl[1].strip()
    return constructor

import inspect
import importlib.util
import sys

class DataStructure():
    def __init__(self, name=None, comp=None, out=None, compute_barrier_function=None, color=None,
                 from_constructor=True,path=".",paths=None,ts=None,gas_path=None,compute_barrier_source=None):
        if not from_constructor:
            self.name                     = name
            self.comp                     = comp
            self.op                       = out
            self.paths                    = paths
            self.ts                       = ts
            self.compute_barrier_function = compute_barrier_function
            self.compute_barrier_source   = compute_barrier_source
            self.color                    = color
        else:
            self.name = os.path.basename(os.path.normpath(path))
            gp = gas_path
            if gas_path is None: gp = path.replace(self.name,"gas_contcar")
            self.op   = read_outputs(path,gp)
            t = self.op[0].temperature
            p = self.op[0].gas_pressures
            self.comp = read_comps(t,p,self.op,path)
            constructor = read_data_structure_class(path)
            for attr in ["color","paths","ts"]:
                if attr in constructor:
                    setattr(self,attr,constructor[attr])
            spec = importlib.util.spec_from_file_location("compute_barrier_function",os.path.join(path,"compute_barrier_function.py"))
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            self.compute_barrier_function = module.compute_barrier_function
    
    def unpickle(self,dir="."):
        metal = self.name
        if not os.path.isdir(os.path.join(dir,"structures",metal)):
            os.mkdir(os.path.join(dir,"structures",metal))
        f = open(os.path.join(dir,"structures",metal,metal+".dscdata"),"w")
        for attr in ["color","paths","ts"]:
            if hasattr(self,attr):
                f.write(attr + " " + str(getattr(self,attr))+"\n")
        f.close()
        barrier_func=self.compute_barrier_source
        barrier_func=barrier_func.replace(self.compute_barrier_function.__name__,"compute_barrier_function",1)
        with open(os.path.join(dir,"structures",metal,"compute_barrier_function.py"),"w") as f:
            f.write(barrier_func)
        if self.paths is not None:
            mode = "w"
            for path in self.paths:
                self.get_structures(self.paths[path],self.ts[path],self.comp,metal,self.op,mode=mode,dir=dir)
                mode = "a"
    
    def flat(self,array):
        anew = []
        for e1 in array:
            for e2 in e1:
                for e3 in e2:
                    #for e4 in e3:
                    anew.append(e3)
        return anew

    def get_tref(self,dname,op):
        for o in op:
            if o.dname == dname:
                return o.thermo_ref_dname

    def get_structures(self, path, ts, comps, metal, op, mode="w",dir="."):
        new_cnames = self.flat(path) + [ts[v] for v in ts]
        comps_dat = dict()
        op_dict = {o.dname: o for o in op}

        def get_all_refs(start_dname):
            names = set()

            def recurse(dname):
                if not dname or dname == "None" or dname in names:
                    return
                names.add(dname)
                if dname not in op_dict:
                    return
                o = op_dict[dname]
                # >>> CHANGED: recurse through scalar or list-valued thermo references
                for atr in ["thermo_ref_dname", "thermo_ref_dname_ts"]:
                    if hasattr(o, atr) and getattr(o, atr) is not None:
                        refs = getattr(o, atr)
                        if not isinstance(refs, (list, tuple)):
                            refs = [refs]
                        for ref in refs:
                            recurse(str(ref))
                # <<< CHANGED

            recurse(start_dname)
            return names

        for combo in new_cnames:
            for term in combo.split("+"):
                match = TERM_RE.fullmatch(term)
                if not match:
                    raise ValueError(f"Ungültiger Term: '{term}'")

                n_str, c = match.groups()
                n = int(n_str) if n_str else 1

                if c in comps:
                    comps_dat[comps[c].nickname] = " ".join([comps[c].nickname, comps[c].adsorbate, str(comps[c].ts), str(comps[c].comp_ref)]) + "\n"

                    dname = comps[c].metal + "_" + comps[c].facet + "_" + comps[c].adsorbate
                    ref_dname = comps[c].metal + "_" + comps[c].facet + "_" + str(comps[c].comp_ref)

                    needed_dnames = set()
                    needed_dnames.update(get_all_refs(dname))
                    needed_dnames.update(get_all_refs(ref_dname))

                    metal2 = comps[c].metal; facet2 = comps[c].facet
                    
                    import re
                    base = f"{metal2}_{facet2}_CO"
                    pattern = re.compile(rf"^{re.escape(base)}(?:\+(?:\d+)?CO)?$")
                    needed_dnames.update(d for d in op_dict if pattern.fullmatch(d))

                    for odname in needed_dnames:
                        if odname not in op_dict:
                            continue

                        o = op_dict[odname]

                        f = open(os.path.join(dir,"structures", metal, o.dname + ".txt"), "w")
                        for atr in ["dname", "compound", "project", "facet", "gas_pressures", "temperature", "adv_thermo", "energy", "eslab", "freqs", "thermo_ref_dname", "ts_energy", "ts_freqs", "thermo_ref_dname_ts", "neb_nimages", "thermo_ref_multiplicator", "thermo_ref_ts_multiplicator", "gas_names", "gas_energies", "gas_frequencies", "gas_geometries", "gas_symmetrie_numbers", "ef_slab"]:
                            if hasattr(o, atr) and getattr(o, atr) is not None:
                                f.write(atr + " " + str(getattr(o, atr)) + "\n")
                        f.close()

                        if o.dname == dname:
                            o.atoms.write(os.path.join(dir,"structures", metal, dname + ".traj"))

                            if comps[c].ts:
                                o.transition_state.write(os.path.join(dir,"structures", metal, dname + "_ts.traj"))
                                traj = Trajectory(os.path.join(dir,"structures", metal, dname + "_neb.traj"), "w")
                                for image in o.neb:
                                    traj.write(image)
                                traj.close()

                        else:
                            o.atoms.write(os.path.join(dir,"structures", metal, o.dname + ".traj"))

        fcomp = open(os.path.join(dir,"structures", metal, metal + ".compdata"), mode=mode)
        for c in comps_dat:
            fcomp.write(comps_dat[c])
        fcomp.close()


def print_comps(comps):
    for comp in comps:
        c = comp
        eads ="    "; gads = "    "
        if c.energy is not None:
            eads = round(c.energy,2)
        if c.free_energy is not None:
            gads = round(c.free_energy,2)
        name = c.nickname + " "*(10-len(c.nickname))
        dname = "_".join([c.metal,c.facet,c.adsorbate])
        print("\t".join([name,str(eads),str(gads),dname,c.comp_ref]))


def update_comps_pt (compounds,outputs,t,p,cref="H2O"):
    ge,gfe = get_gas_data(t,p,outputs)
    if outputs[0].project=="Steam-Reforming": cref="CH4_reform"
    elif outputs[0].project=="CO2RR": cref="CO2RR"
    en,gn,ge,gfe = get_relative_gas_data(ge,gfe,cref=cref)
    op2 = [outputs]; comps2 = [compounds]
    for op in op2:
        op = update_outputs_gibbs_energy(op,t,p)
    count=0
    for c2 in comps2:
        for adsorbate in c2:
            c = c2[adsorbate]
            comps2[count][adsorbate] = Compound(c.adsorbate,c.ts,c.comp_ref,c.nickname,op2[count],ge,gfe)
        count+=1
    compounds = comps2[0]
    outputs   = op2[0]
    return compounds,outputs

def flat(array):
    anew = []
    for e1 in array:
        for e2 in e1:
            for e3 in e2:
                #for e4 in e3:
                anew.append(e3)
    return anew

def get_comp_elist(path,ts,compounds,outputs,t,p,cref="H2O",update_comp=True):
    if outputs[0].project=="Steam-Reforming":
        cref="CH4_reform"
    elif outputs[0].project=="CO2RR":
        cref="CO2RR"
    if update_comp:
        compounds, outputs = update_comps_pt(compounds,outputs,t=t,p=p)
    ge,gfe = get_gas_data(t,p,outputs)
    en,gn,ge,gfe = get_relative_gas_data(ge,gfe,cref=cref)
    new_cnames = flat(path) + [ts[v] for v in ts]
    count=0
    pfe = dict()
    for combo in new_cnames:
        pfe = free_energy_of_combo(combo,compounds,pfe,gn)
    for c in gn:
        pfe[c+"(g)"] = gn[c]
    count+=1
    return pfe
    

def get_gas_data(temperature,gas_pressures,outputs):
    gas_energies    = outputs[0].gas_energies
    for gs in outputs[0].gas_structures:
        if gs=="C2H2": continue
        outputs[0].gas_pressures[gs] = gas_pressures[gs]
    outputs[0].temperature = temperature
    outputs[0].gibbs_gas(update=True)
    gas_free_energies = dict()
    for gs in outputs[0].gas_structures:
        gas_free_energies[gs] = outputs[0].gas_free_energies[gs]
    return gas_energies, gas_free_energies

def comp_erel(gas_free_energies,cref):
    g_new=dict()
    if cref == "CH4_reform":
        g_new["CH4"]   = 0
        g_new["H2O"]   = 0
        g_new["H2"]    = 0
        g_new["C2H4"]  = gas_free_energies["C2H4"]  - (gas_free_energies["CH4"]   + gas_free_energies["C"])
        g_new["C2H6"]  = gas_free_energies["C2H6"]  - (gas_free_energies["CH4"]   + gas_free_energies["C"] + gas_free_energies["H2"])
        g_new["C3H6"]  = gas_free_energies["C3H6"]  - (gas_free_energies["CH4"]   + 2*gas_free_energies["C"] + gas_free_energies["H2"])
        g_new["C3H8"]  = gas_free_energies["C3H8"]  - (2*gas_free_energies["CH4"] + gas_free_energies["C"])
        g_new["C4H10"] = gas_free_energies["C4H10"] - (2*gas_free_energies["CH4"] + 2*gas_free_energies["C"] + gas_free_energies["H2"])
        g_new["CH3OH"] = gas_free_energies["CH3OH"] - (gas_free_energies["H2O"]   + gas_free_energies["C"] + gas_free_energies["H2"])
        g_new["CO2"]   = gas_free_energies["CO2"]   - (gas_free_energies["C"]     + 2*gas_free_energies["O"])
        g_new["CO"]    = gas_free_energies["CO"]    - (gas_free_energies["C"]     + gas_free_energies["O"])

    if cref == "CO2RR":
        g_new["H2"]     = 0
        g_new["CO2"]    = 0
        g_new["C2H4"]   = gas_free_energies["C2H4"]    -  2*gas_free_energies["H2"]   - 2*gas_free_energies["C"]
        g_new["C2H6"]   = gas_free_energies["C2H6"]    -  3*gas_free_energies["H2"]   - 2*gas_free_energies["C"]
        g_new["CH4"]    = gas_free_energies["CH4"]     -  2*gas_free_energies["H2"]   - gas_free_energies["C"]
        g_new["H2O"]    = gas_free_energies["H2O"]     -    gas_free_energies["H2"]   - gas_free_energies["O"]
        g_new["CO"]     = gas_free_energies["CO2"]     -    gas_free_energies["O"]
        g_new["C3H8"]   = gas_free_energies["C3H8"]    -  4*gas_free_energies["H2"]   - 3*gas_free_energies["C"]
        g_new["C3H6"]   = gas_free_energies["C3H6"]    -  3*gas_free_energies["H2"]   - 3*gas_free_energies["C"]
        g_new["CH3OH"]  = gas_free_energies["CH3OH"]   -  2*gas_free_energies["H2"]   - gas_free_energies["C"]  - gas_free_energies["O"]
        g_new["C4H10"]  = gas_free_energies["C4H10"]   -  5*gas_free_energies["H2"]   - 4*gas_free_energies["C"]

    if cref in ["CH4","H2O","0.5co"]:
        g_new["H2"]   = 0
        g_new["CO"]   = 0
        g_new["C2H4"] = gas_free_energies["C2H4"] - 2*gas_free_energies["H2"] - 2*gas_free_energies["C"]
        g_new["C2H6"] = gas_free_energies["C2H6"] - 3*gas_free_energies["H2"] - 2*gas_free_energies["C"]
        g_new["CH4"]  = gas_free_energies["CH4"] -  2*gas_free_energies["H2"] - gas_free_energies["C"]
        g_new["H2O"]  = gas_free_energies["H2O"] -    gas_free_energies["H2"] - gas_free_energies["O"]
        if "CO2" in gas_free_energies:
            g_new["CO2"]  = gas_free_energies["CO2"] -    gas_free_energies["CO"] - gas_free_energies["O"]
        if "C3H8" in gas_free_energies:
            g_new["C3H8"]  = gas_free_energies["C3H8"] -    4*gas_free_energies["H2"] - 3*gas_free_energies["C"]
        if "C3H6" in gas_free_energies:
            g_new["C3H6"]  = gas_free_energies["C3H6"] -    3*gas_free_energies["H2"] - 3*gas_free_energies["C"]
        if "CH3OH" in gas_free_energies:
            g_new["CH3OH"]  = gas_free_energies["CH3OH"] -    2*gas_free_energies["H2"] - gas_free_energies["CO"]
        if "C4H10" in gas_free_energies:
            g_new["C4H10"]  = gas_free_energies["C4H10"] -    5*gas_free_energies["H2"] - 4*gas_free_energies["C"]
    return g_new

def get_relative_gas_data(gas_energies,gas_free_energies,cref="H2O"):
    if cref=="CH4_reform":
        for g in [gas_free_energies,gas_energies]:
            g["C"] = g["CH4"] - 2* g["H2"]
            g["O"] = g["H2O"] - g["H2"]
            g["H"] = g["H2"]/2

    if cref=="CO2RR":
        for g in [gas_free_energies,gas_energies]:
            g["C"] = g["CO2"] + 2*g["H2"] - 2*g["H2O"]
            g["O"] = g["H2O"] - g["H2"]
            g["H"] = g["H2"]/2

    if cref in ["H2O","CH4","0.5co"]:
        for g in [gas_free_energies,gas_energies]:
            g["H"]   = g["H2"]/2
            g["C"]   = g["CO"] + g["H2"] - g["H2O"]
            g["O"]   = g["H2O"] - g["H2"]
            if cref == "0.5co":
                g["C"] = g["CO"]/2
                g["O"] = g["CO"]/2
            if cref == "ch4":
                g["C"] = g["CH4"] - 2* g["H2"]
                g["O"] = g["CO"] + 2*g["H2"] - g["CH4"]

    e_new = comp_erel(gas_energies,cref)
    g_new = comp_erel(gas_free_energies,cref)

    return e_new,g_new,gas_energies,gas_free_energies

def collect_freqs(o):
    freqs = []
    if o.freqs is not None: freqs = list(o.freqs)

    if hasattr(o, "thermo_ref") and o.thermo_ref is not None:
        for i in range(0,o.thermo_ref_multiplicator):
            freqs += collect_freqs(o.thermo_ref)

    return freqs

def collect_gibbs(o,temperature):
    o.temperature = temperature
    g_total = 0
    if o.freqs is not None:
        g_total = o.calculate_gibbs(o.freqs) + o.energy
    else:
        g_total = o.energy

    if hasattr(o, "thermo_ref") and o.thermo_ref is not None:

        # >>> CHANGED: multiple thermo_ref support
        if isinstance(o.thermo_ref, list):
            for ref, mult in zip(o.thermo_ref, o.thermo_ref_multiplicator):
                g_ref = collect_gibbs(ref,temperature)
                g_total += (g_ref-ref.energy)*mult

        # <<< CHANGED
        else:
            g_ref = collect_gibbs(o.thermo_ref,temperature)
            g_total += (g_ref-o.thermo_ref.energy)*o.thermo_ref_multiplicator

    return g_total

def update_outputs_gibbs_energy(outputs,temperature,pressures): 
    for o in outputs:
        o.temperature = temperature
        o.gas_pressures = pressures
        #gibbs0 = o.gibbs_energy
        o.gibbs_energy = collect_gibbs(o,temperature)
        o.gibbs_gas(update=True)
        o.eads, o.count, o.gads  = o.get_eads() 
        
        #print(o.gibbs_energy,gibbs0,o.dname)
        if o.transition_state is not None:
            #gibbs0 = o.ts_gibbs_energy
            if hasattr(o, "thermo_ref_ts") and o.thermo_ref_ts is not None:
                if o.ts_freqs is not None:
                    gibbs = o.calculate_gibbs(o.ts_freqs,ts=True) + o.ts_energy

                    # >>> CHANGED: sum one or multiple TS thermo references
                    refs = o.thermo_ref_ts if isinstance(o.thermo_ref_ts, (list, tuple)) else [o.thermo_ref_ts]
                    multiplicators = o.thermo_ref_ts_multiplicator if isinstance(o.thermo_ref_ts_multiplicator, (list, tuple)) else [o.thermo_ref_ts_multiplicator] * len(refs)
                    add_gibbs = sum(mult * (collect_gibbs(ref,temperature) - ref.energy) for ref, mult in zip(refs, multiplicators))
                    # <<< CHANGED

                    o.ts_gibbs_energy = gibbs + add_gibbs
                    o.ts_eads, o.ts_count, o.ts_gads = o.get_eads(ts=True)
            
                #print(o.ts_gibbs_energy,gibbs0,o.dname)
    return outputs

from ase.io.trajectory import Trajectory

def get_structures(outputs,format="traj"):
    
    if not os.path.isdir(format): os.mkdir(format)

    for out in outputs:
        struct = out.atoms
        ts = out.transition_state
        name   = out.dname
        struct.write(format + "/" + name + "." + format)
        if ts is not None:
            ts.write(format + "/" + name + "_ts." + format)
            neb=out.neb
            traj = Trajectory(format + "/" + name + "_neb.traj","w")
            for image in neb:
                traj.write(image)
            traj.close()

def make_simple_mkm_plots(temp,pres,metals,comp,op,path,ts,barrier_function,tstart=450,tstop=1050,pstart=1,pstop=20,
                          colors=None,x_min=-0.3,x_max=0.7,y_min=0,y_max=0.6,names=None,plot_only=False,tplot=True,
                          pplot=True,plot_fname="plot_simpleMKM.pdf",model="simple",model_parameters=None,parallel=False,
                          modes=None, all_barriers=False, barplot_groups=None, barplot_name=None):

    if len(metals)>1:
        for name,obj in [("path",path),("ts",ts),("comp",comp),("op",op)]:
            if not isinstance(obj,dict):
                raise ValueError(f"{name} must be a dict when more than one metal is used")

    if model != "simple":
        path = {metals[0]:path}
        ts   = {metals[0]:ts}
        comp = {metals[0]:comp}
        op   = {metals[0]:op}

    if model_parameters is None: model_parameters = {}

    if parallel:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()
    else:
        comm = None
        rank = 0
        size = 1

    eads = dict()

    if not parallel or rank==0:
        for metal in metals:
            if model == "simple": pres2 = pres[metal]
            else: pres2 = pres
            eads[metal] = get_comp_elist(path[metal],ts[metal],comp[metal],op[metal],temp,pres2)
            if not plot_only:
                if tplot:
                    f = open("rates_"+metal+".txt","w")
                    g = open("selectivity_"+metal+".txt","w")
                    f.write("Temperature TOF/(1/s)\n")
                    g.write("Temperature Selectivity\n")
                    f.close()
                    g.close()
                if pplot:
                    f2 = open("rates_"+metal+"_p.txt","w")
                    g2 = open("selectivity_"+metal+"_p.txt","w")
                    f2.write("p_tot TOF/(1/s)\n")
                    g2.write("p_tot Selectivity\n")
                    f2.close()
                    g2.close()

    if parallel: eads = comm.bcast(eads if rank==0 else None,root=0)

    if model=="simple":
        smkm = SimpleMKM(eads,temp,pres,barrier_function,modes,all_barriers=all_barriers,barplot_groups=barplot_groups,barplot_name=barplot_name)

        if rank==0:
            smkm.plot(x_min=x_min,x_max=x_max,y_min=y_min,y_max=y_max,names=names,fname=plot_fname)

    else:
        smkm = dict()

        for metal in metals:
            params = deepcopy(model_parameters)
            params["temperature"] = temp
            params["gas_pressures"] = pres
            smkm[metal] = model(eads=eads[metal],**params)

    if plot_only:
        if parallel: comm.Barrier()
        return


    # ------------------------------------------------------------------
    # Prepare all temperature-dependent eads in parallel.
    # ------------------------------------------------------------------

    eads_T = {}

    if tplot:
        eads_T_local = {}

        for i,t in enumerate(range(tstart,tstop)):
            if not parallel or i % size == rank:
                eads_T_local[t] = {}

                for metal in metals:
                    if model == "simple": pres2 = pres[metal]
                    else: pres2 = pres
                    comp_t = deepcopy(comp[metal])
                    op_t = deepcopy(op[metal])
                    comp_t,op_t = update_comps_pt(comp_t,op_t,t,pres2)
                    eads_T_local[t][metal] = get_comp_elist(path[metal],ts[metal],comp_t,op_t,t,pres2)
                    print("done",metal,t)

        if parallel:
            eads_T_parts = comm.allgather(eads_T_local)

            for part in eads_T_parts:
                eads_T.update(part)
        else:
            eads_T = eads_T_local


    # ------------------------------------------------------------------
    # Temperature loop: only MKM + output remains here.
    # ------------------------------------------------------------------

    if tplot:
        for t in range(tstart,tstop):
            eads = eads_T[t]

            if model=="simple":
                smkm = SimpleMKM(eads,t,pres,barrier_function)

            else:
                for metal in metals:
                    params = deepcopy(model_parameters)
                    params["initial_guess"] = smkm[metal].theta
                    params["temperature"] = t
                    params["gas_pressures"] = pres
                    smkm[metal] = model(eads=eads[metal],**params)

            if rank==0:
                for metal in metals:
                    f = open("rates_"+metal+".txt","a")
                    g = open("selectivity_"+metal+".txt","a")

                    if model=="simple":
                        f.write(" ".join([str(t),str(smkm.reaction_rate[metal][0]),str(smkm.reaction_rate[metal][1])])+"\n")
                        g.write(" ".join([str(t),str(smkm.selectivity[metal]/100)])+"\n")
                    else:
                        f.write(" ".join([str(t),str(smkm[metal].reaction_rates["CH4"]),str(smkm[metal].reaction_rates["C3H8"])])+"\n")
                        g.write(" ".join([str(t),str(smkm[metal].selectivity)])+"\n")

                    f.close()
                    g.close()

            if parallel: comm.Barrier()


    # ------------------------------------------------------------------
    # Prepare pressure-dependent eads before the pressure loop.
    # Pressure range is normally small, so rank 0 does this once.
    # ------------------------------------------------------------------

    eads_p = {}
    pres_p = {}

    if pplot:
        if not parallel or rank==0:
            for p in range(pstart,pstop):
                pCO  = p/3
                p_H2 = 2*p/3
                pres_i = deepcopy(pres)

                if model=="simple":
                    for metal in metals:
                        pres_i[metal]["H2"] = p_H2
                        if "CO2" in metal: pres_i[metal]["CO2"] = pCO
                        else: pres_i[metal]["CO"]  = pCO
                else:
                    pres_i["H2"] = p_H2
                    if "CO2" in metals[0]:
                        pres_i["CO2"] = pCO
                    else:
                        pres_i["CO"] = pCO

                pres_p[p] = pres_i
                eads_p[p] = {}

                for metal in metals:
                    if model == "simple": pres2_i = pres_i[metal]
                    else: pres2_i = pres_i
                    comp_p = deepcopy(comp[metal])
                    op_p = deepcopy(op[metal])
                    comp_p,op_p = update_comps_pt(comp_p,op_p,temp,pres2_i)
                    eads_p[p][metal] = get_comp_elist(path[metal],ts[metal],comp_p,op_p,temp,pres2_i)

        if parallel:
            eads_p = comm.bcast(eads_p if rank==0 else None,root=0)
            pres_p = comm.bcast(pres_p if rank==0 else None,root=0)


    # ------------------------------------------------------------------
    # Pressure loop.
    # ------------------------------------------------------------------

    if pplot:
        for p in range(pstart,pstop):
            eads = eads_p[p]
            pres_i = pres_p[p]

            if model=="simple":
                smkm = SimpleMKM(eads,temp,pres_i,barrier_function)

            else:
                for metal in metals:
                    params = deepcopy(model_parameters)
                    params["initial_guess"] = smkm[metal].theta
                    params["temperature"] = temp
                    params["gas_pressures"] = pres_i
                    smkm[metal] = model(eads=eads[metal],**params)

            if rank==0:
                for metal in metals:
                    f2 = open("rates_"+metal+"_p.txt","a")
                    g2 = open("selectivity_"+metal+"_p.txt","a")

                    if model=="simple":
                        f2.write(" ".join([str(p),str(smkm.reaction_rate[metal][0]),str(smkm.reaction_rate[metal][1])])+"\n")
                        g2.write(" ".join([str(p),str(smkm.selectivity[metal]/100)])+"\n")
                    else:
                        f2.write(" ".join([str(p),str(smkm[metal].reaction_rates["CH4"]),str(smkm[metal].reaction_rates["C3H8"])])+"\n")
                        g2.write(" ".join([str(p),str(smkm[metal].selectivity)])+"\n")

                    f2.close()
                    g2.close()

            if parallel: comm.Barrier()

    if rank==0:
        plot_rs(metals,colors,names=names,tplot=tplot,pplot=pplot)

def plot_rs(metals,colors=None,names=None,tplot=True,pplot=True,rate_ymin=-5,rate_ymax=None):
    if colors is None:
        colors2 = {"Ni":"green","Ru":"orange","Co":"violet","Co_CO2":"blue"}
        for metal in metals:
            if metal not in colors2: colors2[metal] = "black"
        colors = [colors2[c] for c in metals]
    else:
        colors = [colors[c] for c in metals]
    files = ["rates_" + metal + ".txt" for metal in metals]
    if tplot: plot_rates_and_selectivity(files,colors,labels=metals,names=names,rate_ymin=rate_ymin, rate_ymax=rate_ymax)
    if pplot: plot_rates_and_selectivity(files,colors,labels=metals,variable="p",names=names,rate_ymin=rate_ymin, rate_ymax=rate_ymax)

import matplotlib.pyplot as plt
import numpy as np

def plot_rates_and_selectivity(
    files,
    colors,
    labels=None,
    Tmin=500,
    Tmax=1000,
    rate_ymin=-5,
    rate_ymax=None,
    variable="T",
    pmin = 1,
    pmax = 20,
    names = None
):
    if labels is None:
        labels = [f"File {i+1}" for i in range(len(files))]
    if names is None:
        names = dict()
        for label in labels:
            names[label] = label

    if variable == "p":
        Tmin = pmin; Tmax = pmax
    # ======================
    # RATES-PLOT
    # ======================
    fig_rates, ax_rates = plt.subplots()

    for fname, color, label in zip(files, colors, labels):
        actual_fname = fname
        if variable == "p":
            actual_fname = actual_fname.replace(".txt","_p.txt")
        data = np.loadtxt(actual_fname, skiprows=1)

        T = data[:, 0]
        rate1 = data[:, 1]
        rate2 = data[:, 2]

        mask_T = np.ones_like(T, dtype=bool)
        if Tmin is not None:
            mask_T &= T >= Tmin
        if Tmax is not None:
            mask_T &= T <= Tmax

        T = T[mask_T]
        rate1 = rate1[mask_T]
        rate2 = rate2[mask_T]

        mask1 = rate1 > 0.0
        mask2 = rate2 > 0.0

        ax_rates.plot(
            T[mask1],
            np.log10(rate1[mask1]),
            "-",
            color = color,
            label = names[label] + " – CH$_4$",
        )
        ax_rates.plot(
            T[mask2],
            np.log10(rate2[mask2]),
            "--",
            color=color,
            label=names[label] +" – C$_{2+}$",
        )

    ax_rates.set_xlabel("Temperature/ K")
    if variable == "p": ax_rates.set_xlabel("p$_\\text{tot}$/ bar")
    ax_rates.set_ylabel("log$_{10}$(rate / s$^{-1}$)")
    if rate_ymin is not None: ax_rates.set_ylim(rate_ymin, rate_ymax)
    ax_rates.legend()
    fig_rates.tight_layout()
    if variable=="p":
        fig_rates.savefig("plot_rates_p.pdf")
    else:
        fig_rates.savefig("plot_rates.pdf")
    plt.close(fig_rates)

    # ======================
    # SELECTIVITY-PLOT
    # ======================
    fig_sel, ax_sel = plt.subplots()

    for fname, color, label in zip(files, colors, labels):
        sel_file = fname.replace("rates", "selectivity")
        actual_fname = sel_file
        if variable == "p":
            actual_fname = actual_fname.replace(".txt","_p.txt")
        data = np.loadtxt(actual_fname, skiprows=1)

        T = data[:, 0]
        S = data[:, 1]

        mask_T = np.ones_like(T, dtype=bool)
        if Tmin is not None:
            mask_T &= T >= Tmin
        if Tmax is not None:
            mask_T &= T <= Tmax

        T = T[mask_T]
        S = S[mask_T]

        mask = S >= 0.0
        if not np.any(mask):
            continue

        ax_sel.plot(
            T[mask],
            S[mask],
            color=color,
            label=names[label],
        )

    ax_sel.set_xlabel("Temperature/ K")
    if variable=="p": ax_sel.set_xlabel("p$_\\text{tot}$/ bar")
    ax_sel.set_ylabel("C$_{2+}$ Selectivity")
    ax_sel.set_ylim(0.0, 1.0)
    ax_sel.legend()
    fig_sel.tight_layout()
    if variable=="p":
        fig_sel.savefig("plot_selectivity_p.pdf")
    else:
        fig_sel.savefig("plot_selectivity.pdf")
    plt.close(fig_sel)


