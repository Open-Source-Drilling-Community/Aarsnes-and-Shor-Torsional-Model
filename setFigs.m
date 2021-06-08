%%
% Torsional Drillstring Model
%   Set up Figures for data display
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

function [hAx,hBx,hFig] = setFigs(p, MD)
    %%

    hFig(1) = figure(1);
    clf;
    subplot(121);
    hAx(1) = plot(nan(size(p.xp)),p.xp);
    hold on;

    tauMax = cumtrapz(flipud([p.xp;p.xc]),...
        flipud(-[p.P_thresProfile;p.C_thresProfile]));
    tauMax = flipud(tauMax);
    tauMax = tauMax/1e3;
    plot(tauMax,[p.xp;p.xc])

    ylim([0 p.Lp+p.Lc]);
    xlim([-3 max(tauMax)*1.1]); 
    title('Torque'); xlabel('(kNm)')
    ylabel('MD (m)');
    set (gca,'Ydir','reverse')

    subplot(122); 
    hAx(2) = plot(nan(size(p.xp)),p.xp);
    hold on
    plot( [0 0]+p.o_thres*60/2/pi ,[0 p.Lc+p.Lp],'--k',...
        [0 0]*0-p.o_thres*60/2/pi ,[0 p.Lc+p.Lp],'--k');
    ylim([0 p.Lp+p.Lc]); 
    xlim([-10 120]);
    title('Angular Velocity'); xlabel('(RPM)')
    ylabel('MD (m)');
    set (gca,'Ydir','reverse')

    %%%%%%%%%%%%%%%%%%%%
    hFig(2) = figure(2);
    clf;
    subplot(211); 
    hBx(1) = plot(nan(size(p.xp)),p.xp);
    hold all
    hBx(2) = plot(nan(size(p.xp)),p.xp);
    hBx(5) = plot(nan(size(p.xp)),p.xp);
    hBx(6) = plot(nan(size(p.xp)),p.xp);
    hBx(7) = plot(nan(size(p.xp)),p.xp);
    l = title(sprintf('Top drive angular velocity, $\\omega_{TD}$, bit depth $= %i$m', round(MD))); 
    set(l,'interpreter','latex');
    ylabel('(RPM)')
    xlabel('Time (s)');


    subplot(212); 
    hBx(3) = plot(nan(size(p.xp)),p.xp);
    hold on
    hBx(4) = plot(nan(size(p.xp)),p.xp);
    plot(p.t,p.t*0 + tauMax(1), '--k');

    l = title('Motor torque, $\tau_m$'); set(l,'interpreter','latex');
    ylabel('(kNm)')
    xlabel('Time (s)');







