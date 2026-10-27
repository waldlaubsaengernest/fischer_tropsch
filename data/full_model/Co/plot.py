from fischer_tropsch.output import DataStructure, plot_rs
import os

metals = ["Co"]#,"Co","Ru","Co_CO2","Ru_CO2"]
pres1   = {'CO':6,'H2':12,'H2O':0.00001,'CH4':0.00001,'C2H2':0.00001,'C2H4':0.00001 ,'C2H6':0.00001,'C3H8':0.00001,'C4H10':0.00001, "C3H6":0.00001, "CO2":0.00001, "CH3OH":0.00001}
table = ""

colors = dict()
for metal in metals:
    data                = DataStructure(path=os.path.join("../../structures",metal))
    colors[metal] = data.color

plot_rs(metals,colors=colors,names=None,tplot=True,pplot=True)
