%%
% Torsional Drillstring Model
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


% Clear the Matlab Console
clc
clear;

%% Parameters
mu   = .3;                                          % Static Friction
fRat = 0.65;                                        % Ratio between Static and Dynamic Friction
oThres = 1.5;                                       % Omega Threshold
MD = 2500;                                          % Current Measure Depth
p = parameters_WellB(MD,mu,fRat,oThres);            % Well Parameters
omega_stat = 40 /60*2*pi;                           % Omega Setpoint

figure(11); clf;
plotyy([p.xp;p.xc],[p.PthetaVec;p.CthetaVec], ...
    [p.xp;p.xc],[p.P_thresProfile;p.C_thresProfile]/1e3)
title('Iinclination (rad) and Coulomb torque per meter (kNm/m)')

%% Basic analytics
% Max stored torque
tauMax = cumtrapz(flipud([p.xp;p.xc]),...
    flipud(-[p.P_thresProfile;p.C_thresProfile]));
tauMax = flipud(tauMax);

% Twist
phiMax = cumtrapz( flipud([p.xp;p.xc]), ...
    flipud([tauMax(1:p.Pp)/(p.Jp*p.Gp); tauMax(p.Pp+1:p.Pp+p.Pc)/(p.Jc*p.Gc)] ) );
phiMax = flipud(phiMax);

phiMax_TD = phiMax(1)

% Amplitude
tauMin = tauMax(1)*(1-2*(1-fRat));

%% Sim stuff
p.dt = .01;
T0 = 0;
T1 = 60*1;
t = T0:p.dt:T1;
Nt = numel(t);

% Init vals
op0 = ones(p.Pp,1)*omega_stat *0;
fp0 = zeros(p.Pp,1);
oc0 = ones(p.Pc,1)*omega_stat *0;
fc0 = zeros(p.Pc,1);
Otd0 = 0;
% fc0 = linspace(2.0e-3*p.f_tRat,0,p.Pc);
% fp0 = ones(p.Pp,1)*fc0(1)*p.Jc/p.Jp;
% Otd0 = omega_stat;


% Parse init states
Pt = p.Pp+p.Pc;     % Total number of cells
x0(1:p.Pp)                   = fp0;   % Pipe Strain
x0(p.Pp+1:Pt)                = fc0;   % Collar Strain
x0(Pt+1:Pt+p.Pp)             = op0;   % Pipe Velocity
x0(Pt+p.Pp+1:Pt+p.Pp+p.Pc)   = oc0;   % Collar Velocity
x0(2*Pt+1)                   = Otd0;   % Bit Angular position

x = zeros(2*(p.Pp+p.Pc)+1,Nt);
y = zeros(1,Nt);
x(:,1) = x0;    % Total number of cells

u.wb = 0;
%% Plot
p.t = t;
[hAx,hBx,hFig] = setFigs(p,MD);
%%
Phi_TD = zeros(Nt,12);
omega_sp = 0;
II = 0;
kk = 1; % Video capture increment
plotPeriod = 250;
% F = floor(Nt/plotPeriod);


for k=1:Nt
    omega_sp = omega_sp + (omega_stat-omega_sp)/(1/p.dt);

    % Compute Omega Surface
    Otd = x(2*Pt+1,k);
    u.tau_Motor = (omega_sp-Otd)*38e3 + II *100e3 ;
    II = II + p.dt*(omega_sp-Otd);
    
    %Increment the Model
    [x(:,k+1)] = torsionalWaveStep(p,x(:,k),u);
    
    y(k) = u.tau_Motor;
    Phi_TD(k+1) = Phi_TD(k) + p.dt * x(Pt+1,k);
    
    % Plot Data
    if mod(k,plotPeriod)==0
        yTorque = [x(1:p.Pp,k)*p.Jp*p.Gp; x(1+p.Pp:Pt,k) *p.Jc*p.Gc]/1e3;
        yVelocity = x(Pt+1:2*Pt,k) *60/2/pi;
        
        set(hAx(1),'XData',yTorque,'YData',[p.xp; p.xc]);
        set(hAx(2),'XData',yVelocity,'YData',[p.xp; p.xc]);

        set(hBx(1),'XData',t,'YData',x(2*Pt+1,:)*60/2/pi);
        set(hBx(3),'XData',t,'YData',y(:)/1e3);

        drawnow
        kk = kk+1;
    end
end

%% Finalize Plots
figure(20);
clf;
plot(t,y(:)/1e3)
hold on
plot(p.t,p.t*0 + tauMax(1)/1e3, '--k');
plot(p.t,p.t*0 + tauMin/1e3, '--k');

l = title('Motor torque, $\tau_m$'); set(l,'interpreter','latex');
ylabel('(kNm)')
xlabel('Time (s)');















