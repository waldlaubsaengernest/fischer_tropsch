import matplotlib.pyplot as plt
import numpy as np
import re
import ase.db
from decimal import Decimal
from fischer_tropsch.catlib_mkm import reaction_names,mkm_class
from fischer_tropsch.catlib_addons    import set_path, set_eads, set_ts_gas_pressures,GasReference,ThermoGasReference
import dill as pickle

REACTION_NAMES = reaction_names()

class Surface:
    def __init__(self,formula,facet,reaction="CO_oxidation",gas_ref="(A)",reactions=None,reverse=None,energies_from="default",eads=None,path=None,ts=None,ts2=[],temperature=None,gas_pressures=None,steady_state=False,sites=None,nguess=200,max_cost=1e-8,initial_guess=None,optcyc=2000):
        self.formula = formula
        self.facet = facet
        self.name = self.formula + self.facet
        self.reaction = reaction
        self.reactions = reactions
        self.reverse = reverse
        self.gas_ref = gas_ref
        self.steady_state = steady_state
        if energies_from=="default":
            if reaction=="FT": energies_from = "DFT"
            else: energies_from = "catapp"
        self.energies_from = energies_from
        self.eads = eads
        self.path = path
        self.path_energies = eads
        self.ts  = ts
        self.ts2 = ts2
        self.constants = {"R" : 8.314*10**(-3)/96.485, "h":4.1356677 *10**(-15)}
        self.temperature = temperature
        self.gas_pressures = gas_pressures
        if self.energies_from is not None:
            self.get_reaction_energies()
            self.get_energies()
            self.set_reaction_path()
            self.sites = sites
            self.nguess = nguess
            self.max_cost = max_cost
            self.optcyc = optcyc
            self.initial_guess = initial_guess

    def get_reaction_energies(self):
        if self.reaction in REACTION_NAMES:
            mkm_model = mkm_class(self.reaction)
            mkm = mkm_model(self.eads,None,None)
            self.reactions = mkm.reactions
            self.reverse = [r.reverse for r in self.reactions]
            if self.reactions is not None:
                for reaction in self.reactions:
                    reaction.er, reaction.ea = self.get_er_ea(reaction.er,reaction.ea,False)
        self.set_er_ea()

    def set_er_ea(self):
        self.er = dict()
        self.ea = dict()
        if self.energies_from == "catapp":
            con = ase.db.connect('catapp.db')
            for reaction in self.reactions:
                for row in con.select(a=reaction.a,b=reaction.b,ab=reaction.ab,formula=self.formula,facet=self.facet):
                    reaction.er,reaction.ea = self.get_er_ea(row.er,row.ea,reaction.reverse)
                if reaction.er is None: reaction.er = 0
                if reaction.ea is None: reaction.ea = 0           
        if self.reactions is not None:
            for reaction in self.reactions:
                self.er[reaction.name] = reaction.er
                self.ea[reaction.name] = reaction.ea

    def get_er_ea(self,er,ea,reverse):
        if er is None or ea is None: return 0,0
        if reverse:
            er_new = -1*er
            ea_new = max([0,ea - er,-er])
        else:
            er_new = er
            ea_new = max([0,ea,er])
        return er_new,ea_new

    def get_energies(self):
        if self.reaction in REACTION_NAMES:
            if self.path is None:
                self.eads = set_eads(self.reaction,self.gas_ref,self.er,self.ea)

    def set_reaction_path(self):
        if self.reaction in REACTION_NAMES:
            if self.path is None:
                self.path, self.path_energies, self.ts = set_path(self.reaction,self.eads)
    
    def label_replace(self,label):
        label = formula_to_latex(label) #label.replace("H4","H$_4$").replace("H2","H$_2$").replace("H3","H$_3$").replace("C2","C$_2$")
        #if not "(g)" in label:
        #        label = "*" + label.replace("+","+*").replace("*3H","3*H")
        return label 

class FreeEnergyDiagram:
    def __init__(self,surfaces,titel="Free Energy Diagram",output_name="plot.pdf",xlabel_rotation=45,fontsize=12,reaction=None,colors=None,ylim_top=None,ylim_bottom=None,max_x=None):
        self.titel = titel
        self.xlabel_rotation = xlabel_rotation
        if fontsize is None: fontsize=12
        self.fontsize = fontsize
        self.ylim_top    = ylim_top
        self.ylim_bottom = ylim_bottom
        self.reaction = reaction
        self.colors = colors
        self.max_x = [None] * len(surfaces)
        if max_x is not None:
            self.max_x = max_x
        self.output_name = output_name
        plt.clf()
        self.path = self.initialize_path(surfaces)

    def create_path(self,path,ts=[],e=None,ts2=[]):
        path_plt = []
        count=0
        count2=0
        for p in path:
            for p2 in p:
                count3=0
                for compound in p2:
                    p_plt = dict()
                    p_plt["label"] = compound#.replace("O2","O$_2$").replace("O2","O$_2$")
                    p_plt["x"] = count
                    p_plt["g"] = e[compound]
                    p_plt["ts_g"] = None
                    if compound in ts:
                        p_plt["ts_g"] = e[ts[compound]]
                    if compound in ts2:
                        if ts2[compound] is None: p_plt["ts2_g"] = None
                        else: p_plt["ts2_g"] = e[ts2[compound]]
                    lnk = []
                    if compound == p2[-1] and p2==p[-1]:
                        if len(p)>1:
                            if len(path)> count2+1:
                                lnk.append(count+1)
                                if len(path[count2+1])>1:
                                    lnk.append(count+1 + len(path[count2+1][0]))
                        elif len(path) > count2+1:
                            lnk.append(count+1)
                            if len(path[count2+1])>1:
                                    lnk.append(count+1 + len(path[count2+1][0]))
                    elif len(p2) > count3+1:
                        lnk.append(count+1)
                    p_plt["links"] = lnk
                    path_plt.append(p_plt)
                    count+=1
                    count3+=1
            count2+=1
        for i in range(len(path_plt)):
            y = []
            for l in path_plt[i]["links"]:
                y.append(path_plt[l]["g"])
            path_plt[i]["y"] = y

        return path_plt

    def interpolate(self,y1,y2,ymax,x1,x2):
        xmax = (x1+x2)/2
        converged = False
        cycle=0
        while not converged:
            A = [[x1**2,x1,1],[x2**2,x2,1],[xmax**2,xmax,1]]
            B = [y1,y2,ymax]
            x = np.linalg.solve(A,B)
            #print(cycle,x,xmax)
            if abs(2*x[0]*xmax + x[1])<0.0001:
                converged=True
            xmax = -x[1]/(2*x[0])
            if xmax<x1: xmax=x1+0.01
            if xmax>x2: xmax=x2-0.01
            cycle+=1
            if cycle>20: print(xmax,ymax);return x,xmax
        return x,xmax

    def plot_path(self,path,color,ls,label,width=0.5,lw=1.5,
                TS_quad_plotpoints=40):
        #   
        X_l = path["g"]
        X_r = path["y"]
        x0  = 0

        plt.plot((path["x"]-0.5*width,path["x"]+0.5*width),(X_l,X_l),lw=2.,c=color,ls=ls,label=label)
        
        i = path["x"]
        rxn_label = path["label"]
        coef = 1.0

        x_A = x0+i+0.5*width #connection point left side of rxn and barrier
        x_r = []
        x_C = []

        if path["links"]==[]:
            return
        
        for shift in path["links"]:
            x_r.append(x0 + shift) #center for right side of rxn
            x_C.append(x0 + shift-0.5*width) #connection point barrier and right side of rxn
        
        ts_g = [path["ts_g"]]
        if "ts2_g" in path:
            ts_g.append(path["ts2_g"])

        j=0
        for tsg in ts_g:
            if tsg is None: Xa=0
            else:  Xa = tsg - path["g"]
            for y in X_r:
                if y is None: continue
                if (y-path["g"]) > Xa: Xa = 0
            if Xa<=0:
                linestyle = "--"
                if j==1: linestyle = ":"
                plt.plot((x_A,x_C[j]),(X_l,X_r[j]),ls=linestyle,lw=lw,c=color)
                plt.plot((x_r[j]-0.5*width,x_r[j]+0.5*width),(X_r[j],X_r[j]),lw=lw*1.4,c=color,ls="--")
            else:
                x,xmax = self.interpolate(y1=X_l,y2=X_r[j],ymax=tsg,x1=x_A,x2=x_C[j])
                xs = np.linspace(x_A,x_C[j],TS_quad_plotpoints,endpoint=True)
                ys = x[0] * xs**2 + x[1] * xs + x[2]
                linestyle = "--"
                if j==1: linestyle = ":"
                plt.plot(xs,ys,lw=lw,c=color,ls=linestyle)
                plt.plot((x_r[j]-0.5*width,x_r[j]+0.5*width),(X_r[j],X_r[j]),lw=lw*1.4,c=color,ls="--")
            j+=1

    def initialize_path(self,surfaces):
        count=0  
        for surf in surfaces:
            if self.colors is None:
                self.colors= ["blue","orange","green","red","purple","brown"]
            if self.reaction is None: self.reaction = surf.reaction
            ts2 = []
            if hasattr(surf,"ts2"):
                ts2 = surf.ts2
            path = self.create_path(surf.path,surf.ts,surf.path_energies,ts2)

            formula = surf.formula
            if "_" in surf.formula:
                formula = "_".join([formula_to_latex(f) for f in formula.split("_")]).replace("*","")
            if "-" in surf.formula:
                formula = "-".join([formula_to_latex(f) for f in formula.split("-")]).replace("*","")
            
            i=0

            for p in path:
                if i==0: label = formula + " " + format_surfname(surf.facet.replace("m","-"))
                else:
                    label =""
                if self.max_x[count] is not None:
                    if p["x"] > self.max_x[count]: continue
                self.plot_path(p,self.colors[count],ls="-",label=label,width=0.4,lw=1.3)
                i+=1
            count+=1

            legend = True
            
            if "Free Energy" in self.titel or "FED" in self.titel:
                ylabel = "Free Energy (eV)"
            else:
                ylabel = "Energy (eV)"
            plt.ylabel(ylabel,fontsize=self.fontsize*12/10)
            plt.title(self.titel,fontsize=self.fontsize*16/10)

            if legend:
                plt.legend(fontsize=10)
            
            if self.reaction == "CO_oxidation":
                labels = [p["label"] for p in path]
            else:
                labels = [surfaces[0].label_replace(p["label"]) for p in path]

            ticks = [i for i in range(len(labels))]

            plt.xticks(ticks,labels,rotation=self.xlabel_rotation,fontsize=self.fontsize)
            plt.yticks(fontsize=self.fontsize)
            plt.tight_layout()
            if self.ylim_bottom is not None and self.ylim_top is not None:
                plt.ylim(self.ylim_bottom,self.ylim_top)
            plt.savefig(self.output_name)


TERM_RE = re.compile(
    r"""
    (\d*)                                   
    (
        [A-Za-z][A-Za-z0-9_]*(?:\(g\))?     
        (?:-[A-Za-z][A-Za-z0-9_]*(?:\(g\))?)*  
    )
    """,
    re.VERBOSE
)

from ase.io import Trajectory
import os

class SimpleMKM():
    def __init__(self,eads,temperature, gas_pressures, barrier_function,modes=None,uncertainty=None,**kwargs):
        self.eads             = eads
        self.temperature      = temperature
        self.gas_pressures    = gas_pressures
        self.barrier_function = barrier_function
        self.reaction = dict()
        self.x_dict   = dict()
        self.y_dict   = dict()
        self.barriers = dict()
        self.theta_CO = dict()
        self.theta_H  = dict()
        self.reaction_rate  = dict()
        self.selectivity    = dict()
        self.modes    = modes

        self.uncertainty = uncertainty
        self.barrier_uncertainty = dict()
        self.reaction_rate_uncertainty = dict()
        self.selectivity_uncertainty = dict()

        r = 8.314*10**(-3)/96.485
        for metal in eads:
            self.reaction[metal]  = "FT"
            if "CO2" in metal: self.reaction[metal]  = "CO2RR"
            mode = None
            if self.modes is not None: mode = modes[metal]
            self.barriers[metal], self.theta_CO[metal], self.theta_H[metal] = self.get_barriers_and_theta(eads[metal],barrier_function[metal],self.reaction[metal], mode=mode, gas_pressures=self.gas_pressures[metal])
            print(metal,self.theta_CO[metal],self.theta_H[metal],self.barriers[metal])
            self.x_dict[metal] = (self.barriers[metal]["CO->C"] - self.barriers[metal]["C->CH4"])  + r * self.temperature *np.log(self.theta_H[metal]/self.theta_CO[metal])
            self.y_dict[metal] = (self.barriers[metal]["CO->C"] - self.barriers[metal]["C+C->CCH3"]) + r * self.temperature *np.log(1/self.theta_CO[metal])
            self.reaction_rate[metal] = self.rate(self.barriers[metal]["CO->C"],self.barriers[metal]["C->CH4"],self.barriers[metal]["C+C->CCH3"],self.temperature,theta_CO=self.theta_CO[metal],theta_H=self.theta_H[metal])
            self.selectivity[metal]   = self.z(self.x_dict[metal],self.y_dict[metal],self.temperature)
            if self.uncertainty is not None:
                self.barrier_uncertainty[metal] = {key: self.uncertainty for key in self.barriers[metal]}
                self.reaction_rate_uncertainty[metal] = self.propagate_rate_uncertainty(self.barriers[metal], self.theta_CO[metal], self.theta_H[metal])
                self.selectivity_uncertainty[metal] = self.propagate_selectivity_uncertainty(self.x_dict[metal], self.y_dict[metal])
        
    def plot(self,x_min=-0.5,x_max=0.75,y_min=-0.35,y_max=0.55,fname="plot_simpleMKM.pdf",names=None):
        
        if self.uncertainty is not None:
            fname = fname[:-4] + "_uncertainty_"+str(self.uncertainty)+".pdf" if fname.endswith(".pdf") else fname + "_uncertainty.pdf"

        x_grid = np.linspace(x_min, x_max, 100)
        y_grid = np.linspace(y_min, y_max, 100)
        X, Y = np.meshgrid(x_grid, y_grid)
        Z = self.z(X, Y)

        plt.figure(figsize=(8, 6))
        
        contour = plt.contourf(X, Y, Z, levels=20, cmap='viridis', alpha=0.8)
        plt.colorbar(contour, label='C$_{2+}$ Selectivity/%')

        keys = self.x_dict.keys()
        x_points = [self.x_dict[k] for k in keys]
        y_points = [self.y_dict[k] for k in keys]

        plt.scatter(x_points, y_points, color='red', edgecolors='white', label='Dict Punkte')
        
        # Layout-Details
        plt.xlabel('d$_1$(E$_{a_1}$-E$_{a_2}$) (eV)')
        plt.ylabel('d$_2$(E$_{a_1}$-E$_{a_3}$) (eV)')
        plt.title('Selectivity Plot')
            
        label_positions = []

        for key in self.x_dict.keys():
            px, py = self.x_dict[key], self.y_dict[key]

            # d1 = Ea1-Ea2, d2 = Ea1-Ea3. For independent barrier errors sigma:
            # Cov(d1,d2) = sigma^2 * [[2,1],[1,2]], giving axes sqrt(3)*sigma and sigma at 45 degrees.
            if self.uncertainty is not None:
                from matplotlib.patches import Ellipse
                ellipse = Ellipse((px,py), width=2*np.sqrt(3)*self.uncertainty, height=2*self.uncertainty, angle=45, fill=False, edgecolor='red', linewidth=1.5, alpha=0.8, zorder=4)
                plt.gca().add_patch(ellipse)

            plt.scatter(px, py, color='red', edgecolors='white', zorder=5)

            prt = key if names is None else names[key]

            # Standard: rechts oben
            dx, dy = 5, 5

            # Wenn ein Punkt in der Nähe schon beschriftet wurde -> links
            if any(abs(px-x) < 0.05 and abs(py-y) < 0.05 for x, y in label_positions):
                dx = -5
                ha = 'right'
            else:
                ha = 'left'

            plt.annotate(prt, (px, py),
                        textcoords="offset points",
                        xytext=(dx, dy),
                        ha=ha,
                        fontsize=9,
                        fontweight='bold',
                        bbox=dict(boxstyle='round,pad=0.3', fc='white', alpha=0.6))

            label_positions.append((px, py))

        plt.grid(True, linestyle='--', alpha=0.5)
        plt.savefig(fname)

    def get_barriers_and_theta(self,eads,barrier_function,reaction="FT",mode=None,gas_pressures=None):
        if mode is None:
            barriers = barrier_function(eads)
        else:
            barriers = barrier_function(eads,mode)
        r = 8.314*10**(-3)/96.485
        if reaction == "CO2RR":
            eads["CO"] = 0
            p = gas_pressures["CO2"]
        else:
            p = gas_pressures["CO"]
        k_co = np.exp(-1*eads["CO"]/(r*self.temperature))
        k_h2 = 2*np.exp(-1*eads["H"]/(r*self.temperature))
        theta_CO = k_co*p/(1+k_co*p)
        theta_H  = np.sqrt(k_h2*gas_pressures["H2"])/(1+np.sqrt(k_h2*gas_pressures["H2"]))

        return barriers, theta_CO, theta_H
    
    def propagate_selectivity_uncertainty(self,x,y,npoints=10000):
        phi = np.linspace(0,2*np.pi,npoints,endpoint=False)
        major = np.sqrt(3)*self.uncertainty
        minor = self.uncertainty
        dx = (major*np.cos(phi)-minor*np.sin(phi))/np.sqrt(2)
        dy = (major*np.cos(phi)+minor*np.sin(phi))/np.sqrt(2)
        nominal = float(self.z(x,y,self.temperature))
        values = np.asarray(self.z(x+dx,y+dy,self.temperature),dtype=float)
        lower = nominal - np.min(values)
        upper = np.max(values) - nominal
        return [float(lower),float(upper)]
    
    def propagate_rate_uncertainty(self,barriers,theta_CO=1,theta_H=1):
        import itertools
        nominal = np.asarray(self.rate(barriers["CO->C"],barriers["C->CH4"],barriers["C+C->CCH3"],self.temperature,theta_CO=theta_CO,theta_H=theta_H,verbose=False),dtype=float)
        values = []
        for signs in itertools.product([-1,1],repeat=3):
            ea1 = barriers["CO->C"] + signs[0]*self.uncertainty
            ea2 = barriers["C->CH4"] + signs[1]*self.uncertainty
            ea3 = barriers["C+C->CCH3"] + signs[2]*self.uncertainty
            values.append(self.rate(ea1,ea2,ea3,self.temperature,theta_CO=theta_CO,theta_H=theta_H,verbose=False))
        values = np.asarray(values,dtype=float)
        lower = nominal - np.min(values,axis=0)
        upper = np.max(values,axis=0) - nominal
        return [[float(lower[i]),float(upper[i])] for i in range(len(nominal))]

    def z(self,x,y,t=523):
        r = 8.314*10**(-3)/96.485
        ye = np.exp(y/(r*t)) #* 1/theta_CO
        xe = np.exp(x/(r*t)) #* theta_H/theta_CO
        fraction = 1/ (  1/2 * (1+1/xe) * (np.sqrt(1 + 8* ye/(1+xe)**2)-1) )
        return 100*(1/(1+fraction))
    
    def rate(self,ea1,ea2,ea3,t=523,theta_CO=1,theta_H=1,verbose=True):
        r = 8.314*10**(-3)/96.485
        h = 4.1356677 *10**(-15)
        k1 = r*t/h * np.exp(-ea1/(r*t)) * theta_CO
        k2 = r*t/h * np.exp(-ea2/(r*t)) * theta_H
        k3 = r*t/h * np.exp(-ea3/(r*t))
        theta = (k1 + k2 - np.sqrt((k1+k2)**2 + 8*k3*k1))/(-4*k3)
        r1 = k1*(1-theta)
        r2 = k2*(theta)
        r3 = 2*k3*(theta**2)
        if t==523 and verbose: print(ea1,ea2,theta,r1,r2,r3,100*(r3)/(r3+r2))
        return [r2,r3]
    
def space_string(s):
    s2 = str(round(s,2))
    if not s2.startswith("-"): s2 = " "+s2
    return s2 + " "*(5-len(s2))

def print_surface_table(surfaces):
    print("    ".join(["[eV]  ","Er_1","Ea_1","Er_2","Ea_2","Er_3","Ea_3"]))
    for surf in surfaces:
        reactions = surf.reactions
        print("   ".join([surf.name,"   ".join([space_string(reac.er) + "   " + space_string(reac.ea)  for reac in reactions])]))


def print_thermo_surface_table(surfaces):
    print("    ".join(["[eV]   ","Er_1","Ea_1","Er_2","Ea_2","Er_3","Ea_3","Rates"]) + "            Log10")
    for surf in surfaces:
        reactions = surf.reactions
        print("   ".join([surf.name,"   ".join([space_string(reac.er) + "   " + space_string(reac.ea)  for reac in reactions])]))# + "   " , ",".join([str('%.2E' % Decimal(rate)) for rate in surf.reaction_rates]) + "   " , ",".join([str(round(np.log10(rate),2)) for rate in surf.reaction_rates]) )


class ThermoSurface(Surface):
    def __init__(self, formula, facet, t, gas_pressures,reaction="CO_oxidation",gas_ref="(A)",reactions=None,correct_H = True, correct_p=True):
        self.temperature = t
        self.gas_pressures = gas_pressures
        self.correct_H = correct_H
        self.correct_p = correct_p
        super().__init__(formula, facet, reaction, gas_ref, reactions,temperature=t,gas_pressures=gas_pressures)
        self.get_rate_constants()

    def get_reaction_energies(self):
        if self.reaction in REACTION_NAMES:
            mkm_model = mkm_class(self.reaction)
            mkm = mkm_model(self.eads,None,None)
            self.reactions = mkm.reactions
            self.reverse = [r.reverse for r in self.reactions]
        self.set_er_ea()

    def set_er_ea(self):
        super().set_er_ea()
        for reaction in self.reactions:
            self.thermo_correct(reaction)
    
    def get_energies(self):
        if self.reaction in REACTION_NAMES:
            self.eads = set_eads(self.reaction, self.gas_ref ,self.er, self.ea, thermo=True, correct_H=self.correct_H, correct_p=self.correct_p)

    def thermo_correct(self,reaction):
        r = self.constants["R"]
        if not reaction.reverse:
            educts = [reaction.ab]
            products = [reaction.a,reaction.b]
        else:
            educts = [reaction.a,reaction.b]
            products = [reaction.ab]
        
        if self.reaction in REACTION_NAMES:
            ts_gas_species, gas_species, self.gas_pressures = set_ts_gas_pressures(self.reaction,self.gas_pressures)

        transition_state = reaction.a + "-" + reaction.b
        t = self.temperature
        if transition_state in ts_gas_species:
            reaction.ea += self.correction(self.gas_pressures[transition_state],t)
        for p in products:
            if p in gas_species:
                reaction.er += self.correction(self.gas_pressures[p],t)
        for p in educts:
            if p in gas_species: 
                reaction.er -= self.correction(self.gas_pressures[p],t)
                reaction.ea -= self.correction(self.gas_pressures[p],t)
        reaction.ea = max(0,reaction.ea,reaction.er)
        self.er[reaction.name] = reaction.er
        self.ea[reaction.name] = reaction.ea
    
    def correction(self,p,t):
        r = 8.314*10**(-3)/96.485
        corr = -t*0.002 
        if self.correct_H: corr += 7/2*r*t
        if self.correct_p: corr += r*t*np.log(p)
        return corr


TERM_RE = re.compile(r'(\d*)(.*)')

def formula_to_latex(formula: str, pd=False) -> str:
    if "$" in formula:
        return formula

    terms = []

    for term in formula.split("+"):
        match = TERM_RE.fullmatch(term)
        if not match:
            terms.append(term)
            continue

        coeff, species = match.groups()

        # Gasphase → kein Stern
        is_gas = "(g)" in species

        # Stern schon vorhanden?
        has_star = species.startswith("*")

        # Stern ergänzen, falls nötig
        if not is_gas and not has_star and not pd:
            species = "*" + species

        # Stöchiometrie korrekt platzieren
        if coeff:
            species = f"{coeff}{species}"

        terms.append(species)

    formula = "+".join(terms)

    # nur echte Indizes tiefstellen (keine führende Stöchiometrie)
    formula = re.sub(
        r'(?<=[A-Za-z\)])(\d+(?:\.\d+)?)',
        r'$_{\1}$',
        formula
    )

    return formula

import re

def format_tabelle(input_str, latex=False, split=None):
    griechisch = sorted([
        "alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta", 
        "theta", "iota", "kappa", "lambda", "mu", "nu", "xi", 
        "pi", "rho", "sigma", "tau", "phi", "chi", "psi", "omega"
    ], key=len, reverse=True)

    def transform_cell(s, is_header=False):
        if not latex: return s
        
        # 1. Chemie-Logik für Header (behandelt auch Pfeile wie O->H2O)
        if is_header:
            if "->" in s:
                parts = s.split("->")
                # Verarbeitet beide Seiten: Falls es wie Chemie aussieht -> \ce{}, sonst normal
                new_parts = []
                for p in parts:
                    if re.search(r'[A-Z]', p): # Enthält mindestens einen Großbuchstaben
                        new_parts.append(rf"\ce{{{p}}}")
                    else:
                        new_parts.append(p)
                return r" $\rightarrow$ ".join(new_parts)
            elif re.search(r'([A-Z][a-z]?\d+|\d+[A-Z]|[A-Z]{2,})', s):
                return rf"\ce{{{s}}}"

        # 2. Standard Ersetzungen für alle Zeilen
        s = s.replace("->", r"$\rightarrow$")
        
        for b in griechisch:
            s = re.sub(rf'\b{b}\b', rf'$\\{b}$', s)
        
        if "_" in s and "$" not in s:
            s = re.sub(r'([a-zA-Z0-9\\]+)_([a-zA-Z0-9]+)', r'$\1_{\2}$', s)
        elif "_" in s:
            s = re.sub(r'_([a-zA-Z0-9]+)', r'_{\1}', s)
            
        return s

    # Daten einlesen
    zeilen_raw = [zeile.split() for zeile in input_str.strip().split('\n')]
    if not zeilen_raw: return ""
    
    daten = []
    for i, zeile in enumerate(zeilen_raw):
        daten.append([transform_cell(w, is_header=(i == 0)) for w in zeile])

    def render_block(block_daten):
        if not block_daten: return ""
        num_cols = max(len(r) for r in block_daten)
        
        # Korrigierte Breitenberechnung
        breiten = [0] * num_cols
        for r in block_daten:
            for i, cell in enumerate(r):
                if i < num_cols:
                    breiten[i] = max(breiten[i], len(cell))
        
        res = []
        if latex:
            layout = "c|" + " c" * (num_cols - 1)
            res.append(f"\\begin{{tabular}}{{{layout}}}")
        
        for idx, row in enumerate(block_daten):
            row_full = row + [""] * (num_cols - len(row))
            separator = " & " if latex else "  "
            line = separator.join(row_full[i].ljust(breiten[i]) for i in range(num_cols))
            
            if latex:
                res.append(line + (r" \\ \hline" if idx == 0 else r" \\"))
            else:
                res.append(line)
                
        if latex: res.append(r"\end{tabular}")
        return "\n".join(res)

    max_cols = max(len(z) for z in daten)
    if split is not None and split < max_cols:
        t1 = [r[:split] for r in daten]
        t2 = [[r[0]] + r[split:] for r in daten]
        return render_block(t1) + "\n\n" + render_block(t2)
    
    return render_block(daten)

def format_zahl(n):
    if n == 0:
        return "0.00"
    if abs(n) < 0.01:
        return f"{n:.1e}"
    else:
        return f"{n:.2f}"

def format_surfname(text):
    # 1. Trenne Buchstaben von der Zahlenfolge
    match = re.match(r"([a-zA-Z]+)(.*)", text)
    
    if not match:
        return text
    
    letters = match.group(1).upper()
    indices_raw = match.group(2)
    
    # 2. Ersetze negative Zahlen wie '-1' durch '\overline{1}'
    # \d+ findet eine oder mehrere Ziffern, die direkt auf ein '-' folgen
    formatted_indices = re.sub(r"-(\d)", r"\\overline{\1}", indices_raw)
    
    # 3. Rückgabe im LaTeX-Mathemodus (umschlossen von $)
    return f"{letters} (${formatted_indices}$)"

def read_pickle(pkl):
    with open(pkl,"rb") as fin:
        eads = pickle.load(fin)
    return eads

def write_pickle(comps,pkl):
    with open(pkl,"wb") as fout:
        pickle.dump(comps,fout)
    return

def coverage_plots(op, t=523, ntot=None, outliers=None, groups=None, filename="plot_coverage.pdf",colors=None):
    import re
    import os
    import matplotlib.pyplot as plt
    import numpy as np

    adsorbate = "CO"
    ntot = {} if ntot is None else ntot
    outliers = {} if outliers is None else outliers
    plot_groups = [None] if groups is None else groups
    colors = {} if colors is None else colors

    for group in plot_groups:
        labels = list(op) if group is None else [label for label in op if group in label]

        fig, ax = plt.subplots()

        for label in labels:
            outputs = op[label]

            if not outputs:
                continue

            metal = outputs[0].compound
            facet = outputs[0].facet
            nsites = ntot.get(label, ntot.get(metal, 8))
            excluded = outliers.get(label, outliers.get(metal, []))

            base = f"{metal}_{facet}_{adsorbate}"
            pattern = re.compile(rf"^{re.escape(base)}(?:\+(?:(\d+)?{adsorbate}))?$")

            data = []

            for out in outputs:
                if pattern.fullmatch(out.dname) is None:
                    continue

                suffix = out.dname[len(f"{metal}_{facet}_"):]
                if suffix in excluded:
                    continue

                if adsorbate not in out.count or out.gads is None:
                    continue

                count = out.count[adsorbate]

                if count/nsites > 1.2:
                    continue

                sconf = 0
                #if count in [0, nsites]:
                #    sconf = 0
                #else:
                #    theta = count / nsites
                #    sconf = -nsites * 8.314e-3 / 96.485 * (theta * np.log(theta) + (1-theta) * np.log(1-theta))

                data.append((count, out.gads - t * sconf))

            data.sort(key=lambda x: x[0])

            xg = []
            yg = []

            for i in range(len(data) - 1):
                dn = data[i+1][0] - data[i][0]

                if dn <= 0:
                    continue

                xg.append(data[i+1][0] / nsites)
                yg.append((data[i+1][1] - data[i][1]) / dn)

            if xg:
                color = colors.get(label, colors.get(metal, None))
                ax.plot(xg, yg, "x--", markersize=10, label=label, color=color)

        ax.set_ylabel(r"$\Delta G_{ads}$ (eV)", fontsize=15)
        ax.set_xlabel(r"$\theta_{CO}$", fontsize=15)
        ax.axhline(y=0, color="black", linestyle="--", linewidth=1)
        ax.legend()

        plt.tight_layout()

        if group is None:
            outfile = filename
        else:
            root, ext = os.path.splitext(filename)
            outfile = f"{root}_{group}{ext}"

        plt.savefig(outfile)
        plt.close()
