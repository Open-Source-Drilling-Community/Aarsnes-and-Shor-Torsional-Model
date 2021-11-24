##
# Torsional Drillstring Model
#   Model Update
#
# @Authors: Ulf Jakob Aarsnes, Roman Shor, Jonathan McIntyre
# 
# Copyright 2021 Open Drilling Project
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software 
# and associated documentation files (the "Software"), to deal in the Software without restriction, 
# including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, 
# and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, 
# subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all copies or substantial 
# portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT 
# NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
#

import math
import numpy as np

def torsionalWaveStep(p, x0, p_dt, u_tau_Motor):

    # Define convenience parameters
    dxp = p.Lp / p.Pp
    dtp = dxp / p.c_t * 0.99    # As per the CFL cond.
    dxc = p.Lc / p.Pc
    dtc = dxc / p.c_t * 0.99;        # As per the CFL cond.

    n = math.ceil(p_dt / min(dtp, dtc))
    dt = p_dt / n
    
    # Total number of cells
    Pt = p.Pp + p.Pc

    # Parse input
    fp  = x0[0:p.Pp]               # Pipe strain
    fc  = x0[p.Pp:Pt]              # Collar strain
    op  = x0[Pt:Pt+p.Pp]           # Pipe angular velocity
    oc  = x0[Pt+p.Pp:Pt+p.Pp+p.Pc] # Collar angular velocity
    Otd = x0[2*Pt]                 # Top drive angular velocity
    
    # Compute Riemann invariants
    Palpha = op + p.c_t * fp # Pipe
    Pbeta  = op - p.c_t * fp
    Calpha = oc + p.c_t * fc # Collar
    Cbeta  = oc - p.c_t * fc
    ##
    t = 0
    for i in range(n):
        
        op = 1 / 2 * (Pbeta + Palpha)
        oc = 1 / 2 * (Cbeta + Calpha)
        
        ## Interfaces
        
        # Topside BC
        Palpha0 = -Pbeta[0] + 2 * Otd
        
        # Bit rock interaction
        tb = 0
        CbetaPp1 = Calpha[p.Pc - 1] - 2 * tb * p.c_t / (p.Jc * p.Gc)
        
        # Pipe - collar interface
        Zbar = p.Jc * p.Gc / (p.Jp * p.Gp)
        PbetaPp1 = 1 / (1 + Zbar) * (Palpha[p.Pp - 1] * (1 - Zbar) + 2 * Zbar * Cbeta[0])
        Calpha0  = 1 / (1 + Zbar) * (2 * Palpha[p.Pp - 1] - Cbeta[0] * (1 - Zbar))
        
        ## Augment Riemann invariants with interface values
        Palpha_pad = np.insert(Palpha, 0, Palpha0)
        Pbeta_pad  = np.append(Pbeta, PbetaPp1)
        Calpha_pad = np.insert(Calpha, 0, Calpha0)
        Cbeta_pad  = np.append(Cbeta, CbetaPp1)
        
        ## Source terms
        Pvf = dt * p.kt * (Pbeta + Palpha)   # Pipe viscous friction
        Cvf = dt * p.kt * (Cbeta + Calpha)   # Pipe viscous friction

        # Collar Coulomb inclusion as Torque per meter source terms
        Pci = (Palpha / dt - 1 / dxp * p.c_t * np.diff(Palpha_pad) + Pbeta / dt + 1 / dxp * p.c_t * np.diff(Pbeta_pad) - 2 * Pvf) / 2
        Pci = Pci * p.Jp * p.rho
        PThreshold = p.P_thresProfile
        PipeFc = (abs(op) < p.o_thres) * np.maximum(np.minimum(Pci, PThreshold), -PThreshold) + (abs(op) >= p.o_thres) * np.sign(op) * PThreshold * p.f_tRat
        
        # Collar Coulomb inclusion as Torque per meter source terms
        Cci = (Calpha / dt - 1 / dxc * p.c_t * np.diff(Calpha_pad) + Cbeta / dt  + 1 / dxc * p.c_t * np.diff(Cbeta_pad ) - 2 * Cvf) / 2
        Cci = Cci * p.Jc * p.rho
        CThreshold = p.C_thresProfile
        CollarFc = (abs(oc) < p.o_thres) * np.maximum(np.minimum(Cci, CThreshold), -CThreshold) + (abs(oc) >= p.o_thres) * np.sign(oc) * CThreshold * p.f_tRat
        
        # Riemann source terms
        S_P = Pvf + PipeFc / (p.Jp * p.rho)
        S_C = Cvf + CollarFc / (p.Jc * p.rho)
    
        ## 1st order Upwind
        Palpha = Palpha - dt / dxp * p.c_t * np.diff(Palpha_pad) - dt * S_P 
        Pbeta  = Pbeta  + dt / dxp * p.c_t * np.diff(Pbeta_pad)  - dt * S_P 
        Calpha = Calpha - dt / dxc * p.c_t * np.diff(Calpha_pad) - dt * S_C 
        Cbeta  = Cbeta  + dt / dxc * p.c_t * np.diff(Cbeta_pad)  - dt * S_C 

        ## Solve ODE for top drive
        tauTD = p.Jp * p.Gp / (2 * p.c_t) * (Palpha[0] - Pbeta[0])
        Otd = Otd + 1 / p.I_TD * dt * (u_tau_Motor - tauTD)

        t = t + dt
    
    ##
    op = 1/2         * (Pbeta + Palpha)
    fp = 1/(2*p.c_t) * (Palpha - Pbeta)
    oc = 1/2         * (Cbeta + Calpha)
    fc = 1/(2*p.c_t) * (Calpha - Cbeta)
    
    # Parse states
    x = np.zeros((2 * Pt + 1))
    x[0:p.Pp]               = fp   # Pipe Strain
    x[p.Pp:Pt]              = fc   # Collar Strain
    x[Pt:Pt+p.Pp]           = op   # Pipe Velocity
    x[Pt+p.Pp:Pt+p.Pp+p.Pc] = oc   # Collar Velocity
    x[2*Pt]                 = Otd  # Bit Angular position

    return x