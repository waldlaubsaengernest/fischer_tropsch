from fischer_tropsch.catlib import FreeEnergyDiagram, Surface
from fischer_tropsch.output import get_comp_elist, DataStructure
import os

pres={'CO':6,'H2':12,'H2O':0.001,'CH4':0.001,'C2H2':0.001,'C2H4':0.001 ,'C2H6':0.001,'CO2':0.001,'C3H8':0.001,'C4H10':0.001, "C3H6":0.001, "CH3OH":0.001, "CO2":6}
metals = ["Ni_211","Co","Ru","Co_CO2","Ru_CO2"]

for metal in metals:
    if not os.path.isdir(os.path.join("..","structures",metal)): continue
    data = DataStructure(path=os.path.join("..","structures",metal))
    path  = data.paths["full"]
    ts    = data.ts["full"]
    color = data.color
    eads = get_comp_elist(path,ts,data.comp,data.op,523,pres)
    surf = Surface(metal.replace("_211",""),data.op[0].facet,"FT",eads=eads,energies_from=None,path=path,ts=ts)
    gases = ["CO","H2"]
    if "CO2" in metal: gases= ["CO2","H2"]
    titel   = "FED at T=" + str(523) + " K, " + ", ".join(["p$_{"+g.replace("2","_2")+"}$=" +str(round(pres[g],2)) for g in gases]) + " bar" 

    if not "CO2" in metal:
        ylim_bottom = -5.5; ylim_top = 1.5;fontsize=12
    else:
        ylim_bottom = -5.5; ylim_top = 1.6;fontsize=10
        if "hc" in metal:
            ylim_bottom = -1.1; ylim_top = 1.9;fontsize=12

    FreeEnergyDiagram([surf],titel,"plot_"+metal.replace("_211","")+"_full.pdf",colors=[color],ylim_bottom=ylim_bottom,ylim_top=ylim_top,xlabel_rotation=90,fontsize=fontsize)

