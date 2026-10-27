from fischer_tropsch.output import make_simple_mkm_plots, DataStructure, get_comp_elist
from fischer_tropsch.catlib import format_tabelle
from fischer_tropsch.ft_co import FT_new
import os
import dill
from copy import deepcopy

from mpi4py import MPI

comps           = dict()
op              = dict()
barrier_func    = dict()
path            = dict()
ts              = dict()
colors          = dict()
names           = dict()

rank = MPI.COMM_WORLD.Get_rank()

if rank == 0:
    from threadpoolctl import threadpool_info
    print(threadpool_info(), flush=True)

metals = ["Co"]#,"Co","Ru","Co_CO2","Ru_CO2"]
pres1   = {'CO':6,'H2':12,'H2O':0.00001,'CH4':0.00001,'C2H2':0.00001,'C2H4':0.00001 ,'C2H6':0.00001,'C3H8':0.00001,'C4H10':0.00001, "C3H6":0.00001, "CO2":0.00001, "CH3OH":0.00001}
table = ""

for metal in metals:

    data                = DataStructure(path=os.path.join("../../structures",metal))
    comps[metal]        = data.comp
    op[metal]           = data.op
    path                = data.full_path
    ts                  = data.full_ts

    pres = deepcopy(pres1)

    eads   = get_comp_elist(path,ts,data.comp,data.op,523,pres,update_comp=False)
    model = FT_new
    mkm    = model(eads,temperature=523,gas_pressures=pres,nguess=10000,max_cost=2.5e-06,min_rates={("CH3_s2","CCCH3"):1e-8})
    table +=mkm.print_head
    table += metal+"    "+mkm.print_line

    if rank==0:
        print(format_tabelle(table,split=len(mkm.theta),latex=False))

if rank==0:
    print(format_tabelle(table,split=len(mkm.theta),latex=False))
    
#make_simple_mkm_plots(temp=523,pres=pres
#,metals=metals,comp=comps,op=op,path=path,ts=ts,barrier_function=barrier_func,colors=colors,names=names,x_max=0.85,y_min=0.0,y_max=0.65,plot_only=True,tplot=True,pplot=False)
