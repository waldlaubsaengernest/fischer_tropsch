from fischer_tropsch.output import DataStructure,get_comp_elist
from fischer_tropsch.catlib import SimpleMKM
import os

comps           = dict()
op              = dict()
barrier_func    = dict()
path            = dict()
ts              = dict()
colors          = dict()
names           = dict()
eads            = dict()
pressures       = dict()

metals = ["Ni","Co","Ru","Co_CO2","Ru_CO2"]
temp = 523

for metal in metals:
    data = DataStructure(path=os.path.join("..","structures",metal))
    comps[metal]        = data.comp
    op[metal]           = data.op
    barrier_func[metal] = data.compute_barrier_function
    path[metal]         = data.paths["vsimple"]
    ts[metal]           = data.ts["vsimple"]
    colors[metal]       = data.color
    names[metal]        = metal.replace("CO2","CO$_2$")
    pressures[metal]    = {'CO':6,'H2':12,'H2O':0.001,'CH4':0.001,'C2H2':0.001,'C2H4':0.001 ,'C2H6':0.001,'C3H8':0.001,'C4H10':0.001, "C3H6":0.001, "CO2":0.001, "CH3OH":0.001}
    if "CO2" in metal: pressures[metal] = {'CO':0.001,'H2':12,'H2O':0.001,'CH4':0.001,'C2H2':0.001,'C2H4':0.001 ,'C2H6':0.001,'C3H8':0.001,'C4H10':0.001, "C3H6":0.001, "CO2":6, "CH3OH":0.001}

    eads[metal]         = get_comp_elist(path[metal],ts[metal],comps[metal],op[metal],temp,pressures[metal])

smkm = SimpleMKM(eads,temp,gas_pressures=pressures,barrier_function=barrier_func,modes=None,uncertainty=0.1)
x_max=0.55;x_min=-0.4;y_min=-0.3;y_max=0.55
smkm.plot(x_min=x_min,x_max=x_max,y_min=y_min,y_max=y_max,names=names,fname="plot_simpleMKM.pdf")
    
