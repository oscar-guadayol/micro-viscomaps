function [figAbsVisc, figRelVisc] = MVM_viscosity_map(fileName,...
    empiricalViscosity, referenceViscosity, N, R,...
    binLength, frameHeight, frameWidth, scaling,...
    NThreshold, RThreshold)

% Creates and exports two maps, one for absolute dynamic viscosities and the other for
% viscosities relative to the background viscosity.
%
% input variables:
%
%       fileName: is a char variable with the prefix for the exported
%           figures.
%
%       empiricalViscosity is an MxP matrix of the estimated dynamic
%           viscosities (in [kg/m/s]), where M an P are the number of
%           partitions of binLength side in the vertical and the horizontal
%           respectively.
%
%       referenceViscosity is the reference dynamic viscosity (in [kg/m/s]).
%
%       N is an MxP matrix of the number of tracks registred in each
%           partition;
%
%       R is an MxP matrix of the R squared of the least square fitting to
%           least squares fitting the data to MSD(τ)=2dDτ.
%
%       binLength is the side of the squared bin in which the frames are
%           partitioned in microns
%
%       frameHeight is the height in pixels of the frame;
%
%       frameWidth is the height in pixels of the frame;
%
%       beadsDiameter is the diameter in microm of the microspheres.
%
%       scaling in pixels/microm;
%
%       NThreshold and RThreshold (optional) are values for the minimum
%       number of tracks and the minimum R squared for the value in the bin
%       two be plotted. Defauls are NThreshold= 20 and RThreshold = 0.95;
%
% output variables:
%
%       figAbsVisc and figRelVisc are the handles for the absolute and the
%           relative viscosity figures respectively.
%
% EXAMPLE:
%
% [figAbsVisc, figRelVisc] = MVM_viscosity_map(fileName,...
%   empiricalViscosity, theoreticalViscosity, N, R, binLength, frameHeight,
%   frameWidth, scaling);
% 
%  Copyright (C) 2020,  Oscar Guadayol
%  oscar_at_guadayol.cat
%
%
%  This program is free software; you can redistribute it and/or modify it
%  under the terms of the GNU General Public License, version 3.0, as
%  published by the Free Software Foundation.
% 
%  This program is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%  General Public License for more details.
% 
%  You should have received a copy of the GNU General Public License along
%  with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
%  This file is part of micro-viscomaps, a collection of matlab scripts to
%  generate maps of viscosity from a video of microspheres displaying
%  Brownian motion, using Multiple Particle Tracking microrhology (MPTM).

if ~exist('NThreshold','var') || isempty('NThreshold')
    NThreshold = 20;
end
if ~exist('RThreshold','var') || isempty('RThreshold')
    RThreshold = 0.95;
end

xBins = (0:binLength*1e-6:frameWidth/scaling*1e-6)+binLength*1e-6/2;
yBins = (0:binLength*1e-6:frameHeight/scaling*1e-6)+binLength*1e-6/2;
empiricalViscosity(N<NThreshold) = nan;
empiricalViscosity(R<RThreshold) = nan;

%% Plot style
plot_style.Version= '1';
plot_style.Format= 'png';
plot_style.Preview= 'none';
plot_style.Width= '11';
plot_style.Height= 'auto';
plot_style.Units= 'centimeters';
plot_style.Color= 'RGB';
plot_style.Background= 'w';
plot_style.FixedLineWidth= '1';
plot_style.ScaledLineWidth= 'auto';
plot_style.LineMode= 'none';
plot_style.LineWidthMin= '1';
plot_style.PSLevel= '2';
plot_style.Renderer= 'painters';
plot_style.Resolution= '600';
plot_style.LineStyleMap= 'none';
plot_style.ApplyStyle= '0';
plot_style.Bounds= 'loose';
plot_style.LockAxes= 'on';
plot_style.ShowUI= 'on';
plot_style.SeparateText= 'off';

%% Absolute viscosity figure
figAbsVisc = figure;

ax(1) = pcolor(repmat(xBins(:),1,length(yBins)),...
    repmat(yBins,length(xBins),1), empiricalViscosity);

axis equal
shading flat
cbar = colorbar;

ax(1).Parent.XLabel.String = '\mum';
ax(1).Parent.YLabel.String = '\mum';
ax(1).Parent.XLim = [0 xBins(end)];
ax(1).Parent.YLim = [0 yBins(end)];
ax(1).Parent.YDir = 'reverse';
ax(1).Parent.XTickLabel = ax(1).Parent.XTick*1e6;
ax(1).Parent.YTickLabel = flip(ax(1).Parent.YTick*1e6);
ax(1).Parent.ColorScale = 'log';
cbarLimits(1) = 10^floor(log10(prctile(empiricalViscosity(:),1)));
cbarLimits(2) = 10^ceil(log10(prctile(empiricalViscosity(:),99)));
ax(1).Parent.CLim = cbarLimits;
cbar.Ticks(1) = cbarLimits(1);
cbar.Ticks(end) = cbarLimits(2);
cbar.TickLabels{1} = ['\leq' cbar.TickLabels{1}];
cbar.TickLabels{end} = ['\geq' cbar.TickLabels{end}];
ylabel(cbar,'Dynamic viscosity [kg/m/s]')
hgexport(gcf,[fileName '_absolute.png'], plot_style)

%% Relative viscosity figure

figRelVisc = figure;
relativeViscosities = empiricalViscosity./referenceViscosity;
ax(2) = pcolor(repmat(xBins(:),1,length(yBins)), repmat(yBins,length(xBins),1),relativeViscosities);
axis equal
shading flat
cbar(2) = colorbar;

ax(2).Parent.XTickLabel = ax(2).Parent.XTick*1e6;
ax(2).Parent.YTickLabel = flip(ax(2).Parent.YTick*1e6);
ax(2).Parent.XLabel.String = '\mum';
ax(2).Parent.YLabel.String = '\mum';
ax(2).Parent.XLim = [0 xBins(end)];
ax(2).Parent.YLim = [0 yBins(end)];
ax(2).Parent.YDir = 'reverse';
cbarLimits(3) = 1;
cbarLimits(4) = round(prctile(relativeViscosities(:),99)*10/5)*5/10;
ax(2).Parent.CLim = cbarLimits(3:4);
cbar(2).Ticks(1) = cbarLimits(3);
cbar(2).Ticks(end) = cbarLimits(4);
ax(2).Parent.ColorScale = 'log';
cbar(2).TickLabels{1} = ['\leq' cbar(2).TickLabels{1}];
cbar(2).TickLabels{end} = ['\geq' cbar(2).TickLabels{end}];
ylabel(cbar(2),'Relative viscosity')

hgexport(figRelVisc,[fileName '_relative.png'], plot_style)
