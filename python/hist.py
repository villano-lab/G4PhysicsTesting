#Import Libraries
try:
    import uproot
    defaultmethod = 'root'
except(ModuleNotFoundError):
    print('Dependency `uproot` for `hist` not found. Only text-based files will be available.')
    defaultmethod = 'csv'
try:
    import pandas
except(ModuleNotFoundError):
    if defaultmethod == 'csv': #If uproot also didn't exist.
        raise ModuleNotFoundError("Local module `hist.py` found neither `pandas` nor `uproot`. Please install one of these two dependencies to continue.")
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../')
#print(sys.path)
import nc_kinematics as nck
import lindhard as lin
import R68_yield as R68y

#############
#Definitions#
#############

#Define the main function we'll be re-using
def f(model,material='Si',k=0.178,q=0.00075):
    if model=='Lindhard':
        if material in ['Si', 'si', 'silicon', 'Silicon']:
            return lin.getLindhardSi_k(k)
        elif material in ['Ge', 'ge', 'germanium', 'Germanium']:
            return lin.getLindhardGe_k(k)
    elif model=='Sorenson':
        return lambda Er: R68y.ySor(Er,k,0.00075)
    else:
        print("Unrecognized. Available models: Lindhard, Sorenson")

#Get total deposit for the given cascade:
def Eitot(i,l,en,en_dep,c_id,model,material='Si'): #Get yield total from a given cid and k, and choose an instance of that cascade
    #Determine values
    energy = en[c_id==i][l]
    delEnergy = en_dep[c_id==i][l]
    #Create array and iterate
    returnval = 0
    j = 0
    while j < len(energy):
        returnval += f(model,material)(energy[j])*energy[j]
        returnval -= f(model,material)(energy[j] - delEnergy[j])*(energy[j] - delEnergy[j])
        j += 1
    return returnval

#Determine what lands inside
threshold = 50

onestep = lambda method,material: f(method,material='Si')(nck.D1s)*nck.D1s

#Fitted resolution function
def res(E,scalefactor=1):
    res0 = 10.55
    Fano = 0.6111
    eps  = 3.8
    return np.sqrt(res0**2 + Fano*eps*E)*scalefactor

def histogramable(file,binsize=8,binmin=0,binmax=425,labels=[],model='Lindhard',material='Si',method=defaultmethod,resolution='normal',val=None,scalefactor=1,seed=None): #Build a histogram from the file
    if seed != None:
        np.random.seed(seed)
    if method=='root':
        x = uproot.open(file)
        cas = x['cascade']
        en = cas.array('E')
        en_dep = cas.array('delE')
        c_id = cas.array('cid')
        totalpoints = cas.numentries
    elif method=='csv':
        print('Placeholder function. Sorry!')
    else:
        raise ValueError("Bad value for `method`. Available methods: `root`, `csv`.")
    if resolution=='normal':
        resolution = lambda E: res(E,scalefactor)
    elif resolution =='none':
        resolution = lambda E: 0 #Define as function of E to keep code down to one line @ 81
    elif resolution == 'offset':
        if val is None:
            offset = float(input('Choose your added value for the resolution'))
        else:
            offset = val
        resolution = lambda E: res(E,scalefactor) + offset
    elif resolution == 'constant':
        if val is None:
            resolution = lambda E: float(input('Choose the resolution'))
        else:
            resolution = lambda E: val
    else:
        print("Unrecognized. Available resolutions: normal, none, offset, constant")
    maxcidlength = 0
    for i in range(max(c_id)):
        if len(c_id[c_id==i]) > maxcidlength:
            maxcidlength = len(c_id[c_id==i])
    plottable = np.zeros([max(c_id)+1,maxcidlength]) 
    for i in range(max(c_id)+1):
        for j in range(len(c_id[c_id==i])):                   #Up to the number of events,
            plottable[i,j] = Eitot(i,j,en,en_dep,c_id,model,material)  #Calculate the energy deposit.
        for j in range(len(c_id[c_id==i]),len(c_id[c_id==0])):#For any leftover points,
            plottable[i,j] = np.nan                           #The value is not a number
    a_list = []
    for i in range(max(c_id)):
        a_list.append(plottable[i,:])
        labels.append('ID: '+str(i))
    for i,vec in np.ndenumerate(a_list):
        x,y = i #Restructure enumeration
        a_list[x][y] += np.random.normal(scale=resolution(a_list[x][y]))
    return a_list,labels
"""
    bins = np.arange(binmin,binmax,binsize)
    n_sim,nx_sim = np.histogram(a_list,bins)
    xc = (nx_sim[:-1] + nx_sim[1:]) / 2 
    width = (nx_sim[1] - nx_sim[0])
    n_rel = n_sim/(len(a_list)*width)
    print(len(a_list[0]))
    return xc,n_rel,n_sim
"""