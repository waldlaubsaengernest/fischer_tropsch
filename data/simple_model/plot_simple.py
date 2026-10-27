from fischer_tropsch.output import make_simple_mkm_plots, DataStructure,get_comp_elist
import os

comps           = dict()
op              = dict()
barrier_func    = dict()
path            = dict()
ts              = dict()
colors          = dict()
names           = dict()
pressures       = dict()

metals = ["Ni","Co","Ru","Co_CO2","Ru_CO2"]

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

make_simple_mkm_plots(temp=523,pres=pressures,metals=metals,comp=comps,op=op,path=path,ts=ts,barrier_function=barrier_func,colors=colors,names=names,x_max=0.45,x_min=-0.25,y_min=-0.2,y_max=0.45,plot_only=False,tplot=True,pplot=True)
