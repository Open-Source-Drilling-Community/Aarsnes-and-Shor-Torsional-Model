%%
% Torsional Drillstring Model
%   Parameters for Well A
%
% @Authors: Ulf Jakob Aarsnes and Roman Shor
% 
% Copyright 2021 Open Drilling Project
%
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software 
% and associated documentation files (the "Software"), to deal in the Software without restriction, 
% including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, 
% and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, 
% subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all copies or substantial 
% portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT 
% NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
% IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
% WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
% SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
%
%

function p = parameters_WellA(MD,mu,fRat,oThres)

    %%
    % Physical parameters
    LpV = [19.080   10.240      57.240      9.750       95.400      38.750];
    ODV = [0.150167 0.206375    0.150167    0.206248    0.150167    0.222392];
    IDV = [0.096799 0.076200    0.096799    0.069850    0.096799    0.095964];

    % Averaged polar moment of inertia and cross sectional area for BHA
    p.Ac = sum(LpV.*pi.*( (ODV/2).^2 - (IDV/2).^2 )) / sum(LpV);
    p.Jc = pi*sum( LpV.*pi.*(ODV.^4 -  IDV.^4)/32 ) / sum(LpV);
    p.Cro = sum(LpV.*ODV/2 ) / sum(LpV);

    % Pipe
    p.Pro = .146529/2;     %[m] Drill string inner radius
    p.Pri = .123024/2;     %[m] Drill string outer radius
    p.Jp = pi/2*(p.Pro^4-p.Pri^4);  %[m^4] Drill string polar moment of inertia
    p.Ap = pi*(p.Pro^2-p.Pri^2);    %[m^2] Drill string cross sectional area

    p.Gp = 61e9;         % [m] Pipe shear modulus
    p.Gc = 67e9;         % [m] Collar shear modulus
    p.rho = 7850;        % [kg/m3] Density

    p.I_TD = 2900;  % [kgm^2] Actual top drive measured mass moment of inertia

    p.kt   = 50;                      %[-] Torsional damping

    % computed quantities
    p.c_t = sqrt(p.Gp/p.rho);      %[m/s] Torsional wave velocity

    %% Length parameters
    p.Lc    = 230;         %[m] Drill collar length
    p.Lp    = MD - p.Lc;   %[m] Drill pipe length

    %% Numerics params
    Pt = 1500;   % Number of cells in discretization
    p.Pp = round(Pt * p.Lp/(p.Lp+p.Lc));    % Pipe cells
    p.Pc = Pt - p.Pp;                       % Collar cells

    p.xp    = linspace(0,p.Lp,p.Pp).';
    p.xc    = p.Lp + linspace(0,p.Lc,p.Pc).';

    %% Inclination vectors
    load datasets/wellProfile
    p.PthetaVec  = interp1(wellProfile.MD,...
        wellProfile.inc,p.xp)*pi/180;
    p.CthetaVec  = interp1(wellProfile.MD,...
        wellProfile.inc,p.xc)*pi/180;
    %% Coulomb friction params
    p.f_thres   = mu;
    p.f_tRat    = fRat;
    p.o_thres   = oThres;

    g = 9.81; % Acceleration of gravity

    p.P_thresProfile   = g*sin(p.PthetaVec)*p.rho*p.Ap*p.Pro *p.f_thres; % Torque per meter [Nm/m]
    p.C_thresProfile   = g*sin(p.CthetaVec)*p.rho*p.Ac*p.Cro *p.f_thres; % Torque per meter [Nm/m]
