#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on

@author: a.mikkonen@iki.fi
"""
import time
import scipy as sp
from scipy import optimize
from matplotlib import pyplot as plt
import os
import json
#import pyswarm
import shutil
import glob
from natsort import natsorted
import scipy.interpolate as interpolate
from PyPDF2 import PdfFileMerger


s2min = 60
s2h   = 60**2
s2d   = 24*60**2



#######
# TODO
##########3

# init different
# materials


def get_data(label):
    
    material_root_path = "measurementData/materialProperties/"
    #A4_Isover_RKL-EJ/wc.txt
    material_paths = glob.glob(material_root_path+"*")
    
#    print(label, material_paths)
    
#    for material_path in material_paths[-4:]:
#    for material_path in material_paths[:1]:
    
    
    hankalat = [
            "measurementData/materialProperties/C5_Luja_A",
            "measurementData/materialProperties/A4_Isover_RKL-EJ",
            "measurementData/materialProperties/A7_Vital-levy",
            "measurementData/materialProperties/A3_Isover_RKL",
            "measurementData/materialProperties/D5_Pellavaeriste_T3",
            "measurementData/materialProperties/A12_Tuulensuojaluja",
            "measurementData/materialProperties/D3_Vital",
            ]
    
    for material_path in material_paths:
        if os.path.split(material_path)[-1].split("_")[0] == label:
            break
         
#    print(material_path)
    
#    name = os.path.split(material_path)[-1]
    
    RHr = []
    wr = []
    
    with open(os.path.join(material_path, "wc.txt"),"r") as ifile:
        lines = ifile.readlines()
    for line in lines[2:]:
#            print(line)
        parts = [round(float(part),4) for part in line.split(",")]
        RHr.append(parts[0])
        wr.append(parts[1])
        
#    if len(wr) == 1:
#        continue
#    if sp.isclose(wr[-1],0):
#        continue
    
    
    if wr[-1] >= 2*wr[-2]:
        w = wr[:-1]
        RH = RHr[:-1]
    else:
        w = wr[:]
        RH = RHr[:]
   
    
    

    # 3 
    
    
    if material_path in hankalat:
        w = [(w[0]+w[1])/2, w[2]] + w
        RH = [0.22, 0.51] + RH
        xpoints = [min(RH)]*4 + [0.61] + [max(RH)]*4
    else:
        w = [w[0]] + w
        RH = [0.1] + RH
        xpoints = [min(RH)]*4 + [0.8] + [max(RH)]*4
    
    
#        w = [(w[0]+w[1])/2, w[2]] + w
#        RH = [0.22, 0.51] + RH
#        xpoints = [min(RH)]*4 + [0.61] + [max(RH)]*4
    
    
    def tomin(x):
        tcki = [xpoints,
            list(x) + [0]*5,
            3]
        
        
        return (w - interpolate.splev(RH, tcki, der=0))#*penalty


    x0 = [1]*5
#    print(tomin(x))
    res = optimize.least_squares(tomin, x0)
    
    
    
    tck = [xpoints,
            list(res.x) + [0]*5,
            3]
    
#        RHsp = sp.linspace(0,RH[-1],1000)
    RHsp = sp.linspace(min(RH),max(RH),1000)
    
    wsp = interpolate.splev(RHsp, tck, der=0)
    xisp = interpolate.splev(RHsp, tck, der=1)
    
    # Filter
    wsp[wsp<0] = 0
    
    
    
    
    if sp.isclose(RHsp[-1],1):
        RHextra = []    
        wextra = []
        xiextra = []
    else:
        RHextra = sp.linspace(RHsp[-1]+(RHsp[-1]-RHsp[-2]), 2, 200)
#            RHextra = sp.linspace(0.8, 1, 100000)
        
        #
        xi0 = xisp[-1]#(wsp[-1]-wsp[-2]) / (RHsp[-1]-RHsp[-2])

        x0 = RHsp[-1]
        x1 = RHr[-1]
        y0 = wsp[-1]
        y1 = wr[-1]

        a = (x0**2*x1*xi0 + x0**2*y1 - x0*x1**2*xi0 - 2*x0*x1*y0 + x1**2*y0)/(x0**2 - 2*x0*x1 + x1**2)
        b = (-x0**2*xi0 + 2*x0*y0 - 2*x0*y1 + x1**2*xi0)/(x0**2 - 2*x0*x1 + x1**2)
        c = (-xi0*(x0 - x1) + y0 - y1)/(x0**2 - 2*x0*(x0 - x1) - x1**2)
        
        # funcs
        wextra = a + b*RHextra + c*RHextra**2
        xiextra = b + 2*c*RHextra
            
            
    return sp.array(list(RHsp)+list(RHextra)), sp.array(list(wsp)+list(wextra)), sp.array(list(xisp)+list(xiextra))
            


#
def fitw(material, figurek):
    
#    source_path = os.path.join("measurementData","2a.json")
#    with open(source_path, "r") as ifile:
#        idata = json.load(ifile)
        
#    for k, material in enumerate(idata["materials"][:]):
        
#    label = material["label"]
#    print(label)
    
    w = sp.array(material["w_sorption"])
#    if type(material["w_sorption"]) == list:# and label == "D2":
    RH = sp.array(material["w_sorption_RH"])
    
    wmax = w[-1]
    #    146/(1+(-8e-8*R_water*T_ref*rho_w*sp.log(phi))**1.6)**0.375
    def fit_func(phi, a,b,c):
        return wmax / (1+(a*sp.log(phi))**b)**c 

    
    def to_min(parms):
        return (((w[:-1]-fit_func(RH[:-1], *parms))**2)**0.5).sum()
    try:
        popt, pcov = optimize.curve_fit(fit_func, 
                            RH[:], 
                            w[:],
                            [-2,
                             0.7,
                             2],
#                                [-67744.35379546453,
#                                 0.38825626197169383,
#                                 1.006448450321765]
                            )
    except RuntimeError:
        lb = [-1e4, 0,1]
        ub = [0, 1, 3]
        popt, fopt = pyswarm.pso(to_min, lb, ub,maxiter=1000 )
    
#        print(fopt)
    
    RHfit = sp.linspace(0,1,10000)
    wfit=fit_func(RHfit, *popt)
        
    RHfit = sp.append(RHfit, 1000)
    wfit = sp.append(wfit, (w[-1] - w[-2])/(RH[-1] - RH[-2])* (RHfit[-1] - RHfit[-2]) )
    
#    plt.figure(figurek)
#    plt.plot(RH[:], w[:], "d")
#    plt.plot(RHfit, wfit, "-")
#    plt.xlim(0,1.1)
#    plt.ylim(0,w.max()*1.3)
    
    
    return RHfit, wfit    
    
def writeline(ofile, key, val):
    ofile.write(key+" " + str(val) + ";\n")
    
    
def filter_boundary_RH(data):
     # Filter test
    if not data["wallNro"] in [8, 11 ]: 
        RH_cavityInsideSurface = sp.array(data["RH_cavityInsideSurface"])
        RH_insideSurface = sp.array(data["RH_insideSurface"])
        diffLim = 0.1
        
        
        for k in range(len(data["time_boundary"])):
            if data["time_boundary"][k] > 5*24*60**2:
                firstTimek = k
                break
        for k in range(firstTimek, len(RH_cavityInsideSurface)-1):
            if RH_cavityInsideSurface[k+1] > (1+diffLim) * RH_cavityInsideSurface[k]:
                RH_cavityInsideSurface[k+1] = RH_cavityInsideSurface[k]
            elif RH_cavityInsideSurface[k+1] < (1-diffLim) * RH_cavityInsideSurface[k]:
                RH_cavityInsideSurface[k+1] = RH_cavityInsideSurface[k]
                
        for k in range(firstTimek, len(RH_insideSurface)-1):
            if RH_insideSurface[k+1] > (1+diffLim) * RH_insideSurface[k]:
                RH_insideSurface[k+1] = RH_insideSurface[k]
            elif RH_insideSurface[k+1] < (1-diffLim) * RH_insideSurface[k]:
                RH_insideSurface[k+1] = RH_insideSurface[k]
    else:
        RH_cavityInsideSurface = sp.array(data["RH_cavityInsideSurface"])
        RH_insideSurface = sp.array(data["RH_insideSurface"])
        
    return RH_insideSurface, RH_cavityInsideSurface
    
    

def gen_case(source, target_root, template):
    
#    print(source)
    with open(source, "r") as ifile:
        data = json.load(ifile)
    
    target = os.path.join(target_root, data["id"])
    print(target)
    
    ###########################
    # template
    ###########################
    if os.path.exists(target):
        shutil.rmtree(target)
    shutil.copytree(template, target)
    
#    for L in data["Ls"]:
#        L = round(L,3)
#        
#        
#        print(L)
    
    ###########################
    # parameters
    ###########################
    with open(os.path.join(target, "constant", "parameters"), "w") as ofile:
        # blockMesh
        L = round(sum(data["Ls"]),3)
        dx = 1e-3 
        n = L/dx # 10 ok, 1 causes errors
        if n % 1:
            raise ValueError()
        writeline(ofile, "L", L)
        writeline(ofile, "n", int(n))
        
        # controlDict
        writeline(ofile, "endTime", data["time_boundary"][-1])
#        writeline(ofile, "deltaT", 60*s2min)
#        writeline(ofile, "deltaT", 1)
        writeline(ofile, "deltaT", 60*s2min) # 10-60min ok
#        writeline(ofile, "deltaT", s2min)
        writeline(ofile, "writeInterval", s2d)
        # functionObjects
        writeline(ofile, "x_probeS", data["x_probeS"])
        writeline(ofile, "x_probeU", data["x_probeU"])
        
        # setFields and init values
        writeline(ofile, "T_initial", (data["T_probeS"][0] + data["T_probeU"][0])/2)
        writeline(ofile, "phi_initial ", (data["RH_probeS"][0] + data["RH_probeU"][0])/2)
        writeline(ofile, "phi_initial_frame", data["storage_RH_frame_and_innerBoard"])
        
        # topoSetDict
        if not data["id"] in ["29b","30b", "31b", "32b"]:
            writeline(ofile, "innerBoardLabel", data["structure"][0])
            writeline(ofile, "x_innerBoard",  round(data["Ls"][0],3))
            writeline(ofile, "moistureBarrierLabel", data["structure"][1])
            writeline(ofile, "x_moistureBarrier",  round(data["Ls"][0]+data["Ls"][1],3))
            writeline(ofile, "insulationLabel", data["structure"][2])
            writeline(ofile, "x_insulation",  round(data["Ls"][0]+data["Ls"][1]+data["Ls"][2],3))
            writeline(ofile, "outerBoardLabel", data["structure"][3])
            # ref
            writeline(ofile, "x_refStart", round(data["Ls"][0]-1e-3,3))
            writeline(ofile, "x_refEnd", round(data["Ls"][0]+data["Ls"][1]+1e-3,3))
        else:
            raise ValueError()
            
            
    ##############################################
    # transportProperties
    ##############################################
    with open(os.path.join(target, "constant", "transportProperties"), "r") as ifile:
        transportProperties_lines = ifile.readlines()
    with open(os.path.join(target, "constant", "transportProperties"), "w+") as ofile:
        for line in transportProperties_lines:
            ofile.write(line)
        for material_name in data["structure"]:
            ofile.write("#include \"materials/"+ material_name + "\"\n")
        ofile.write("\n\n// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //")    
           
            
            
    ##########################################
    # Boundary
    ##########################################
    
    
    time_boundary = data["time_boundary"]
    
    T_insideSurface = data["T_insideSurface"]
    T_cavityInsideSurface = data["T_cavityInsideSurface"]
    
    
    #T_init = (T_insideSurface[0] + T_cavityInsideSurface[0])/2
    
#    phi_insideSurface = data["RH_insideSurface"]
#    phi_cavityInsideSurface = data["RH_cavityInsideSurface"]
    
    
    phi_insideSurface, phi_cavityInsideSurface = filter_boundary_RH(data)
    
    
    
    #phi_init = (phi_insideSurface[0] + phi_cavityInsideSurface[0])/2
    
    
#    print(time_boundary[-1])
    # Initial
#    with open(os.path.join(target, "constant","initialConditions"), "w") as ofile:
#        ofile.write("T_initial %f;\n"  % T_init)
#        ofile.write("phi_initial %f;\n"  % phi_init)
    
    # T_insideSurface
    assert len(time_boundary) == len(T_insideSurface)
    with open(os.path.join(target, "constant","boundaryData", "T_insideSurface"), "w") as ofile:
        ofile.write("(\n")
        for k in range(len(time_boundary)-1):
            if time_boundary[k+1] > time_boundary[k]:
                ofile.write("    (%f %f)\n"  % (time_boundary[k], T_insideSurface[k]))    
        ofile.write(");")
            
    # T_cavityInsideSurface
    assert len(time_boundary) == len(T_insideSurface)
    with open(os.path.join(target, "constant","boundaryData", "T_cavityInsideSurface"), "w") as ofile:
        ofile.write("(\n")
        for k in range(len(time_boundary)-1):
            if time_boundary[k+1] > time_boundary[k]:
                ofile.write("    (%f %f)\n"  % (time_boundary[k], T_cavityInsideSurface[k]))    
        ofile.write(");")
            
    # phi_insideSurface
    assert len(time_boundary) == len(phi_insideSurface)
    with open(os.path.join(target, "constant","boundaryData", "phi_insideSurface"), "w") as ofile:
        ofile.write("(\n")
        for k in range(len(time_boundary)-1):
            if time_boundary[k+1] > time_boundary[k]:
                ofile.write("    (%f %f)\n"  % (time_boundary[k], phi_insideSurface[k]))    
        ofile.write(");")
            
    # phi_cavityInsideSurface
    assert len(time_boundary) == len(phi_insideSurface)
    with open(os.path.join(target, "constant","boundaryData", "phi_cavityInsideSurface"), "w") as ofile:
        ofile.write("(\n")
        for k in range(len(time_boundary)-1):
            if time_boundary[k+1] > time_boundary[k]:
                ofile.write("    (%f %f)\n"  % (time_boundary[k], phi_cavityInsideSurface[k]))    
        ofile.write(");")
                
            
    ###################################################################
    # Materials 
    ###################################################################
    figurek = 0
    for material in data["materials"]:
        label = material["label"]
#        print(label)
        with open(os.path.join(target, "constant","materials", label), "w") as ofile:
            ofile.write("%s\n" % label)
            ofile.write("{\n")
            ofile.write("    rho %e;\n" % material["rho"])
            ofile.write("    cp %e;\n" % material["c"])
            
            # k
            if type(material["kt"]) == list:
                ofile.write("    k table ( //phi     k\n")
                for k, phi in zip(material["kt"], material["kt_RH"]):
                    ofile.write("             (%e %e )\n" % (phi, k))
                ofile.write("            );\n")
            else:
                ofile.write("    k %e;\n" % material["kt"])
                
#            # w
#            if type(material["w_sorption"]) == list:
#                ofile.write("    w table ( //phi     w\n")
#                for w, phi in zip(material["w_sorption"], material["w_sorption_RH"]):
#                    ofile.write("             (%e %e )\n" % (phi, w))
#                ofile.write("            );\n")
#            else:
#                ofile.write("    w %e;\n" % material["w_sorption"])
            

            # FIT
            # w
#            if type(material["w_sorption"]) == list:
#                
#                phiwfit, wfit = fitw(material, figurek)
#                figurek += 1
#                ofile.write("    w table ( //phi     w\n")
#                for w, phi in zip(wfit, phiwfit):
#                    ofile.write("             (%e %e )\n" % (phi, w))
#                ofile.write("            );\n")
#            else:
#                ofile.write("    w %e;\n" % material["w_sorption"])
                
            
            # spline regression
            if type(material["w_sorption"]) == list:

                phis, ws, xis = get_data(label)
                
                ofile.write("    w table ( //phi     w\n")
                for w, phi in zip(ws, phis):
                    ofile.write("             (%e %e )\n" % (phi, w))
                ofile.write("            );\n")
            else:
                ofile.write("    w %e;\n" % material["w_sorption"])


                
            # Dw
            if type(material["Dw"]) == list:
                ofile.write("    Dw table ( //phi     Dw\n")
                for Dw, phi in zip(material["Dw"], material["Dw_RH"]):
                    ofile.write("              (%e %e )\n" % (phi, Dw))
                ofile.write("             );\n")
            else:
                ofile.write("    Dw %e;\n" % material["Dw"])

            # delta_p
            if type(material["delta_p"]) == list:
                ofile.write("    delta_p table ( //phi     delta_p\n")
                for delta_p, phi in zip(material["delta_p"], material["delta_p_RH"]):
                    ofile.write("             (%e %e )\n" % (phi, delta_p))
                ofile.write("             );\n")
            else:
                ofile.write("    delta_p %e;\n" % material["delta_p"])
            
            ##############################3
            # xi
            
#            if type(material["w_sorption"]) == list:
#                w = material["w_sorption"]
#                RH = material["w_sorption_RH"]
#                ofile.write("    xi table ( //phi     xi\n")
#                for k in range(1, len(w)):
#                    dw = w[k] - w[k-1]
#                    dRH = RH[k] - RH[k-1]
#                    xi = dw/dRH
#                    ofile.write("              (%e %e )\n" % ((RH[k] + RH[k-1])/2, xi))
#                ofile.write("             );\n")
#            else:
#                ofile.write("    xi %e;\n" % 0)
            
#            if type(material["w_sorption"]) == list:
#                ofile.write("    xi table ( //phi     xi\n")
#                for k in range(1, len(wfit)):
#                    dw = wfit[k] - wfit[k-1]
#                    dRH = phiwfit[k] - phiwfit[k-1]
#                    xi = dw/dRH
#                    ofile.write("              (%e %e )\n" % ((phiwfit[k] + phiwfit[k-1])/2, xi))
#                ofile.write("             );\n")
#            else:
#                ofile.write("    xi %e;\n" % 0)
#                
#            ofile.write("  TEST")
            

            # spline regressioon
            if type(material["w_sorption"]) == list:
                ofile.write("    xi table ( //phi     xi\n")
                for xi, phi in zip(xis, phis):
                    ofile.write("              (%e %e )\n" % (phi, xi))
                ofile.write("             );\n")
            else:
                ofile.write("    xi %e;\n" % 0)





            
            ###############################3
            
            
            # End
            ofile.write("}\n")
        
            
        

def gen_all(source_dir, target_root, template):
    sources = natsorted(glob.glob(os.path.join(source_dir,"*.json")))
    
    for source in sources:
        #if not source.split("/")[-1].split(".")[0] in ["29b","30b", "31b", "32b"]:
        if source.split("/")[-1].split(".")[0] in ["1a", "1b", "2a", "2b","3a", "3b", "4a", "4b"]:
            print(source)
            gen_case(source, target_root, template)
    


    
if __name__ == "__main__":
    start = time.time()
    print("START")
    
    source_dir = "measurementData/byCase/json"
    target_root = "cases"
    template="1DtransientTemplate"
    
    
    gen_all(source_dir,target_root,template)

    print("END %.4f s" % (time.time()-start))


