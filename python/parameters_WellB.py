##
# Torsional Drillstring Model
#   Parameters for Well B
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
from scipy.interpolate import interp1d

class parameters_WellB:
    def __init__(self, MD, mu, fRat, oThres):
        # Physical parameters
        LcV = np.array([48.700, 114.800, 56.700, 38.880])
        ODV = np.array([0.150473, 0.187322,	0.209550, 0.247280])
        IDV = np.array([0.096207, 0.098912,	0.071438, 0.123391])
    
        # Averaged polar moment of inertia and cross sectional area for BHA
        self.Ac = np.sum(LcV * math.pi * ((ODV / 2) ** 2 - (IDV / 2) ** 2)) / np.sum(LcV)
        self.Jc = math.pi * np.sum(LcV * math.pi * (ODV ** 4 - IDV ** 4) / 32) / np.sum(LcV)
        self.Cro = np.sum(LcV * ODV / 2) / sum(LcV)

        # Pipe
        self.Pro = 0.146529 / 2     # [m] Drill string inner radius
        self.Pri = 0.123024 / 2     # [m] Drill string outer radius
        self.Jp = math.pi / 2 * (self.Pro ** 4 - self.Pri ** 4)     # [m^4] Drill string polar moment of inertia
        self.Ap = math.pi * (self.Pro ** 2 - self.Pri ** 2)         # [m^2] Drill string cross sectional area

        self.Gp = 61e9      # [m] Pipe shear modulus
        self.Gc = 67e9      # [m] Collar shear modulus
        self.rho = 7850     # [kg/m3] Density

        self.I_TD = 2900    # [kgm^2] Actual top drive measured mass moment of inertia

        self.kt = 50 * 1    # [-] Torsional damping

        # Computed quantities
        self.c_t = math.sqrt(self.Gp / self.rho)

        # Length parameters
        self.Lc = 230           # [m] Drill collar length
        self.Lp = MD - self.Lc  # [m] Drill pipe length

        # Numerics params
        Pt = 1000       # Number of cells in discretization
        self.Pp = round(Pt * self.Lp / (self.Lp + self.Lc))     # Pipe cells
        self.Pc = Pt - self.Pp                                  # Collar cells

        self.xp = np.linspace(0, self.Lp, self.Pp)
        self.xc = self.Lp + np.linspace(0, self.Lc, self.Pc)

        # Inclination vectors
        MD_prof = np.array([0, 500, 1600, 4000])
        inc_prof = np.array([0, 0, 60, 60])

        self.PthetaVec = interp1d(MD_prof, inc_prof)(self.xp) * math.pi / 180
        self.CthetaVec = interp1d(MD_prof, inc_prof)(self.xc) * math.pi / 180

        # Coulomb friction params
        self.f_thres = mu
        self.f_tRat = fRat
        self.o_thres = oThres

        g = 9.81    # Acceleration of gravity

        self.P_thresProfile = g * np.sin(self.PthetaVec) * self.rho * self.Ap * self.Pro * self.f_thres  # Torque per meter [Nm/m] 
        self.C_thresProfile = g * np.sin(self.CthetaVec) * self.rho * self.Ac * self.Cro * self.f_thres  # Torque per meter [Nm/m]
