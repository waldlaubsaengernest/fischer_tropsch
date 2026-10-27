from fischer_tropsch.catlib import FreeEnergyDiagram, Surface
from fischer_tropsch.output import get_comp_elist, DataStructure
import os
from copy import deepcopy

pres0={'CO':6,'H2':12,'H2O':0.001,'CH4':0.001,'C2H2':0.001,'C2H4':0.001 ,'C2H6':0.001,'CO2':0.001,'C3H8':0.001,'C4H10':0.001, "C3H6":0.001, "CO2":0.001, "CH3OH":0.001}
metals = ["Ni_211","Co","Ru","Co_CO2","Ru_CO2"]

for metal in metals:
    if not os.path.isdir(os.path.join("..","structures",metal)): continue
    data = DataStructure(path=os.path.join("..","structures",metal))
    path = data.paths["vsimple"]
    ts    = data.ts["vsimple"]
    color = data.color
    pres = deepcopy(pres0)
    if "CO2" in metal: pres["CO"]=0.00001; pres["CO2"] = 6
    eads = get_comp_elist(path,ts,data.comp,data.op,523,pres)
    reaction = "FT"
    ylim_bottom = -4.5
    if "CO2" in metal:
        reaction ="CO2RR"
        ylim_bottom = -4.5
    surf = Surface(metal.replace("_211",""),data.op[0].facet,"FT_simple",eads=eads,energies_from=None,path=path,ts=ts)
    gases = ["CO","H2"]
    if "CO2" in metal: gases= ["CO2","H2"]
    titel   = "FED at T=" + str(523) + " K, " + ", ".join(["p$_{"+g.replace("2","_2")+"}$=" +str(round(pres[g],2)) for g in gases]) + " bar" 

    FreeEnergyDiagram([surf],titel,"plot_"+metal.replace("_211","")+"_simple.pdf",colors=[color],ylim_bottom=ylim_bottom,ylim_top=1.5,reaction=reaction,xlabel_rotation=90)
