#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 19:20:11 2019

Steady state 1D validation. 

@author: Antti Mikkonen, a.mikkonen@iki.fi
"""

import scipy as sp
from scipy import optimize
import os
from matplotlib import pyplot as plt
from scipy import interpolate
import glob
import time
from natsort import natsorted
import json
import shutil
from PyPDF2 import PdfFileMerger

C2K = 273.15

def merge_pdfs(paths, out):
    merger = PdfFileMerger()
    for pdf in paths:
        merger.append(pdf)
    merger.write(out)

def read_postline():
    source_dir = os.path.join("postProcessing", "sampleLine")
    time_paths = [os.path.split(part)[-1] for part in glob.glob(source_dir+os.path.sep+"*")]
    times = sp.array([float(part) for part in time_paths])
    index = times.argsort()
    
    xs = []
    Ts = []
    phis = []
    rhos = []
    cps = []
    ks = []
    ws = []
    Dws = []
    delta_ps = [] 
    xis= []
    for time_path in time_paths:
        path = os.path.join(source_dir, time_path)


        with open(os.path.join(path, "out_T_phi_rho_cp_k_w_Dw_delta_p_xi.xy"), "r") as ifile:
            lines = ifile.readlines()
        
            x_cfd = []
            T_cfd = []
            phi_cfd = []
            rho_cfd = []
            cp_cfd = []
            k_cfd = []
            w_cfd = []
            Dw_cfd = []
            delta_p_cfd = []
            xi_cfd = []
            for line in lines:
                parts = [float(part) for part in line.split()]
                x_cfd.append(parts[0])
                T_cfd.append(parts[1])
                phi_cfd.append(parts[2])
                rho_cfd.append(parts[3])
                cp_cfd.append(parts[4])
                k_cfd.append(parts[5])
                w_cfd.append(parts[6])
                Dw_cfd.append(parts[7])
                delta_p_cfd.append(parts[8])
                xi_cfd.append(parts[9])
            xs.append(x_cfd)    
            Ts.append(T_cfd)
            phis.append(phi_cfd)
            rhos.append(rho_cfd)
            cps.append(cp_cfd)
            ks.append(cp_cfd)
            ws.append(w_cfd)
            Dws.append(Dw_cfd)
            delta_ps.append(delta_p_cfd)
            xis.append(xi_cfd)
    xs = sp.array(xs)[index]
    Ts = sp.array(Ts)[index]
    phis = sp.array(phis)[index]
    rhos = sp.array(rhos)[index]
    cps = sp.array(cps)[index]
    ks = sp.array(ks)[index]
    ws = sp.array(ws)[index]
    Dws = sp.array(Dws)[index]
    delta_ps = sp.array(delta_ps)[index]
    Dws = sp.array(Dws)[index]
    
    return times, xs, Ts, phis, rhos, cps, ks, ws, Dws, delta_ps

def plot_lines(times, xs, Ts, phis, rhos, cps, ks, ws, Dws, delta_ps):
    pass

#    
#    
#    fig, ax = plt.subplots(figsize=(7,3))
#    for k in range(len(times[:])):
#        x = xs[k]
#        T = Ts[k]
#        ax.plot(x*1e3,T, "-", label=times[k]/(60*60*24))
#        
#    # Prettyfy
#    ax.grid(True)
#    ax.set_xlim(0, 217)
#    ax.set_xlabel("x (mm)")
#    ax.set_ylabel("T (C)")
#    fig.tight_layout()
#    fig.savefig("profiles.pdf")
#    ax.legend(frameon=False)


def read_T_probes(case_path):
    with open(os.path.join(case_path, "postProcessing", "probes", "0", "T"), "r") as ifile:
            lines = ifile.readlines()
            
    times = []
    T_inner = []
    T_outer = []
    
    for line in lines[4:]:
        parts = [float(part) for part in line.split()]
        times.append(parts[0])
        T_inner.append(parts[1])
        T_outer.append(parts[2])
        
    return sp.array(times), sp.array(T_inner), sp.array(T_outer)

def read_phi_probes(case_path):
    with open(os.path.join(case_path, "postProcessing", "probes", "0", "phi"), "r") as ifile:
            lines = ifile.readlines()
            
    times = []
    phi_inner = []
    phi_outer = []
    
    for line in lines[4:]:
        parts = [float(part) for part in line.split()]
        times.append(parts[0])
        phi_inner.append(parts[1])
        phi_outer.append(parts[2])
        
    return sp.array(times), sp.array(phi_inner), sp.array(phi_outer)
        
def read_probes_T_comsol(comsol_path):
    with open(os.path.join(comsol_path, "comsolOut", "T_S.txt"), "r") as ifile:
            lines_s = ifile.readlines()
    with open(os.path.join(comsol_path, "comsolOut", "T_U.txt"), "r") as ifile:
            lines_u = ifile.readlines()
            
    times = []
    T_inner = []
    T_outer = []
    
    for line in lines_s[8:]:
        parts = [float(part) for part in line.split()]
#        print(parts)
#        times.append(parts[0])
#        T_inner.append(parts[1])
        T_outer.append(parts[1])

    for line in lines_u[8:]:
        parts = [float(part) for part in line.split()]
#        print(parts)
        times.append(parts[0])
        T_inner.append(parts[1])
#        T_outer.append(parts[1])
        
    return sp.array(times)*60*60, sp.array(T_inner), sp.array(T_outer)

def read_probes_phi_comsol(comsol_path):
    with open(os.path.join(comsol_path, "comsolOut", "RH_S.txt"), "r") as ifile:
            lines_s = ifile.readlines()
    with open(os.path.join(comsol_path, "comsolOut", "RH_U.txt"), "r") as ifile:
            lines_u = ifile.readlines()
            
    times = []
    phi_inner = []
    phi_outer = []
    
    for line in lines_s[8:]:
        parts = [float(part) for part in line.split()]
#        print(parts)
#        times.append(parts[0])
#        T_inner.append(parts[1])
        phi_outer.append(parts[1])

    for line in lines_u[8:]:
        parts = [float(part) for part in line.split()]
#        print(parts)
        times.append(parts[0])
        phi_inner.append(parts[1])
#        T_outer.append(parts[1])
        
    return sp.array(times)*60*60, sp.array(phi_inner), sp.array(phi_outer)


#def main():
##    
#    times_T_probe, T_inner_probe, T_outer_probe = read_T_probes()
#    times_phi_probe, phi_inner_probe, phi_outer_probe = read_phi_probes()
#    
#    times_probe_T_comsol, T_inner_probe_comsol, T_outer_probe_comsol = read_probes_T_comsol()
#    times_probe_phi_comsol, phi_inner_probe_comsol, phi_outer_probe_comsol = read_probes_phi_comsol()
#    
#    fig, axes = plt.subplots(4,1, figsize=(7,10))
#    
#    ax = axes[0]
#    ax.plot(times_probe_T_comsol/60/60/24, T_inner_probe_comsol-C2K, "k-", label="COMSOL")
#    ax.plot(times_T_probe/60/60/24, T_inner_probe-C2K, "--", label="OpenFOAM")
#    
#    
#    ax = axes[1]
#    ax.plot(times_probe_T_comsol/60/60/24, T_outer_probe_comsol-C2K, "k-", label="COMSOL")
#    ax.plot(times_T_probe/60/60/24, T_outer_probe-C2K, "--", label="OpenFOAM")
#    
#    ax = axes[2]
#    ax.plot(times_probe_phi_comsol/60/60/24, phi_inner_probe_comsol*100, "k-", label="COMSOL")
#    ax.plot(times_phi_probe/60/60/24, phi_inner_probe*100, "--", label="OpenFOAM")
#    
#    
#    ax = axes[3]
#    ax.plot(times_probe_phi_comsol/60/60/24, phi_outer_probe_comsol*100, "k-", label="COMSOL")
#    ax.plot(times_phi_probe/60/60/24, phi_outer_probe*100, "--", label="OpenFOAM")
#    
#    
#    
#    
#    # Prettyfy
#    for ax in axes.ravel()[:2]:
#        ax.set_xlim(0, (times_T_probe/60/60/24).max())
#        ax.set_ylabel("T (C)")
#        
#    for ax in axes.ravel()[2:]:
#        ax.set_xlim(0, (times_phi_probe/60/60/24).max())
#        ax.set_ylabel("RH (-)")    
#    
#    for ax in axes.ravel():
#        ax.grid(True)
#        ax.set_xlabel("time (d)")
#        ax.legend(frameon=False)
#    
#    fig.tight_layout()
#    fig.savefig("probes.pdf")
    
#    times, xs, Ts, phis, rhos, cps, ks, ws, Dws, delta_ps = read_postline()
#    plot_lines(times, xs, Ts, phis, rhos, cps, ks, ws, Dws, delta_ps)


#    
#    
#    fig, ax = plt.subplots(figsize=(7,3))
#    for k in range(len(times[:])):
#        x = xs[k]
#        T = Ts[k]
#        ax.plot(x*1e3,T, "-", label=times[k]/(60*60*24))
#        
#    # Prettyfy
#    ax.grid(True)
#    ax.set_xlim(0, 217)
#    ax.set_xlabel("x (mm)")
#    ax.set_ylabel("T (C)")
#    fig.tight_layout()
#    fig.savefig("profiles.pdf")
#    ax.legend(frameon=False)

    
def readMinMax(case_path):

    with open(os.path.join(case_path, "log", "solver"), "r") as ifile:
        lines = ifile.readlines()
    
    times = []
    T_min = []
    T_max = []
    phi_min = []
    phi_max = []
    for line in lines:
        # Time
        if line[:7] == "Time = ":
            times.append(float(line.split()[-1])) 
        # T_min
        elif line[:13] == "    min(T) = ":
            T_min.append(float(line.split()[2])) 
        # T_max
        elif line[:13] == "    max(T) = ":
            T_max.append(float(line.split()[2])) 
        # phi_min
        elif line[:15] == "    min(phi) = ":
            phi_min.append(float(line.split()[2])) 
        # phi_max
        elif line[:15] == "    max(phi) = ":
            phi_max.append(float(line.split()[2])) 
    
#    print(times)
#    print(T_min)
#    print(T_max)
#    print(phi_min)
#    print(phi_max)


    return sp.array(times), sp.array(T_min), sp.array(T_max), sp.array(phi_min), sp.array(phi_max)

    
    
    

def onecase(source, case_paths, comsol_path, post_path):
    with open(source, "r") as ifile:
        data = json.load(ifile)
        
    times_T_probes, T_inner_probes, T_outer_probes = [], [], []    
    times_phi_probes, phi_inner_probes, phi_outer_probes = [], [], []    
        
    for case_path in case_paths:
        times_T_probe, T_inner_probe, T_outer_probe = read_T_probes(case_path)
        times_phi_probe, phi_inner_probe, phi_outer_probe = read_phi_probes(case_path)
        #    times_min_max, T_min, T_max, phi_min, phi_max = readMinMax(case_path)
        times_T_probes.append(times_T_probe)
        T_inner_probes.append(T_inner_probe)
        T_outer_probes.append(T_outer_probe)    
        times_phi_probes.append(times_phi_probe)
        phi_inner_probes.append(phi_inner_probe)
        phi_outer_probes.append(phi_outer_probe)
    
    if comsol_root:
        times_probe_T_comsol, T_inner_probe_comsol, T_outer_probe_comsol = read_probes_T_comsol(comsol_path)
        times_probe_phi_comsol, phi_inner_probe_comsol, phi_outer_probe_comsol = read_probes_phi_comsol(comsol_path)
    
    
    fig, axes = plt.subplots(2)
    
    # MEASUREMENT
    ax = axes[0]
    ax.plot(sp.array(data["time_probe"])/60/60/24, sp.array(data["T_probeS"])-C2K, "k-", label="probe")
    ax.plot(sp.array(data["time_probe"])/60/60/24, sp.array(data["T_probeU"])-C2K, "k-", label="probe")
    
    
    
    ax = axes[1]
    ax.plot(sp.array(data["time_probe"])/60/60/24, sp.array(data["RH_probeS"])*100, "k-")#, label="probe")
    ax.plot(sp.array(data["time_probe"])/60/60/24, sp.array(data["RH_probeU"])*100, "k-")#, label="probe")

    # boundary
#    ax = axes[0]
#    ax.plot(sp.array(data["time_boundary"])/60/60/24, sp.array(data["T_insideSurface"])-C2K, "g-", label="boundary")
#    ax.plot(sp.array(data["time_boundary"])/60/60/24, sp.array(data["T_cavityInsideSurface"])-C2K, "g-", label="boundary")
#    ax = axes[1]
#    ax.plot(sp.array(data["time_boundary"])/60/60/24, sp.array(data["RH_insideSurface"])*100, "g-", label="boundary")
#    ax.plot(sp.array(data["time_boundary"])/60/60/24, sp.array(data["RH_cavityInsideSurface"])*100, "g-", label="boundary")
#
#    ax.plot(sp.array(data["time_boundary"])/60/60/24, RH_cavityInsideSurface*100, "c-", label="boundary")
#    ax.plot(sp.array(data["time_boundary"])/60/60/24, RH_insideSurface*100, "c-", label="boundary")




    # hamfoam
#        ax.fill_between(times_min_max/60/60/24, phi_min*100, phi_max*100,alpha=0.2, color="grey", linewidth=0 )#, label="min/max" )


    # COMSOL
    if comsol_root:
        ax = axes[0]
        ax.plot(times_probe_T_comsol/60/60/24, T_inner_probe_comsol-C2K, "b:", label="COMSOL")
        ax.plot(times_probe_T_comsol/60/60/24, T_outer_probe_comsol-C2K, "b:", label="COMSOL")
        ax = axes[1]
        ax.plot(times_probe_phi_comsol/60/60/24, phi_inner_probe_comsol*100, "b:")#, label="COMSOL")
        ax.plot(times_probe_phi_comsol/60/60/24, phi_outer_probe_comsol*100, "b:")#, label="COMSOL")

    #    
#    # hamFOAM
    hamcolors = ["r", "g", "c"]
    ax = axes[0]
    for times_T_probe, T_inner_probe, T_outer_probe, color in zip(times_T_probes, T_inner_probes, T_outer_probes, hamcolors):
        ax.plot(times_T_probe/60/60/24, T_inner_probe-C2K, color+"--", label="hamFOAM")
        ax.plot(times_T_probe/60/60/24, T_outer_probe-C2K, color+"--", label="hamFOAM")
    
    ax = axes[1]
    for times_phi_probe, phi_inner_probe, phi_outer_probe, color in zip(times_phi_probes, phi_inner_probes, phi_outer_probes, hamcolors):
        ax.plot(times_phi_probe/60/60/24, phi_inner_probe*100, color+"--")#, label="hamFOAM")
        ax.plot(times_phi_probe/60/60/24, phi_outer_probe*100, color+"--")#, label="hamFOAM")

    for ax in axes.ravel():
        ax.grid(True)
        ax.set_xlim(0, (sp.array(data["time_probe"])/60/60/24)[-1])
    axes[0].legend(frameon=False, loc='lower left')

    axes[1].set_xlabel("time (d)")
    axes[0].set_ylabel("T ($^\circ$C)")
    axes[1].set_ylabel("RH (%)")
    axes[0].set_title("Wall %s, %s, %s"    % (data["wallNro"],  data["id"], data["structure"]))
    
    # Text
    axes[0].annotate("Data by Tampere University, Building Physics",(0.37,0.03),xycoords='axes fraction')
    axes[1].annotate("alpha testing\na.mikkonen@iki.fi 2019",(0.01,0.055),xycoords='axes fraction')

    fig.tight_layout()
    figname = data["id"]+"_probes.pdf"
    figpath = os.path.join(case_path, figname)
    fig.savefig(figpath)
    #shutil.copyfile(figpath, os.path.join(post_path, figname))

    return figpath
    
        
def post_all(source_dir, calc_roots, comsol_root, post_path):
    sources = natsorted(glob.glob(os.path.join(source_dir,"*.json")))
    
    outfigpaths = []
    for source in sources[:]:
        
#        if not "21b" in source:
#            continue
        
        
        with open(source, "r") as ifile:
            data = json.load(ifile)
    
        
        calc_paths = [os.path.join(calc_root, data["id"]) for calc_root in calc_roots]
        
        if comsol_root:
            comsol_path = os.path.join(comsol_root, data["id"])
        else:
            comsol_path = None
        
        print(calc_paths)
        try:
            outfigpath = onecase(source, calc_paths, comsol_path, post_path)
            outfigpaths.append(outfigpath)
        except FileNotFoundError:
            print("/t FileNotFoundError")
    
    
    
    merge_pdfs(outfigpaths, os.path.join(post_path, "all.pdf"))
    
    
    
           
if __name__ == "__main__":
    start = time.time()
    print("START")
    

    source_dir = "measurementData/byCase/json"
    comsol_root = None 
    
    calc_roots = [
                    "cases",
                  ]
    post_path = "cases"
    
    post_all(source_dir, calc_roots, comsol_root, post_path)

    print("END %.4f s" % (time.time()-start))

    
    
    