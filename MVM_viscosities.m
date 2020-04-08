function [empiricalViscosity, theoreticalViscosity, N, R] =...
    MVM_viscosities(tracks, scaling, frameRate, frameHeight, frameWidth,...
    temperature, salinity,...
    beadsDiameter, binLength, maxLag, microrheologyBias)

% Divides the frame of the video into squares of binLength side, and
% calculates average empirical dynamic viscosities within each one of the
% squares. To estimate viscosities, it first computes the squared
% displacements with time lags (tau) between 1/frameRate and maxLag
% seconds. The squared displacements are then grouped into the partitions
% according to their mean centroid, and the ensemble MSDs is computed for
% each partition. The average diffusivity D of the microspheres within each
% partition is calculated by least squares fitting the data to
% MSD(tau)=2dDtau, where d=2 is the dimensionality of the tracks. The
% average dynamic viscosity for each partition is then calculated from the
% Stokes-Einstein equation for diffusion of a spherical particle. Details
% of the method are given in Guadayol et al 2020.
%
% input variables:
%
%       tracks: is MX3XP matrix the tracks of microspheres within the
%           video. M is the maximum track length in all tracks, and P is
%           the number of tracks. Columns are: 
%              1: x position of the centroid in microns
%              2: y position of the centroid in microns
%              3: projected area of the microsphere in microns^2;
%
%       scaling in pixels/microm;
%
%       frameRate is the video frame rate in frames/s;
%
%       frameHeight is the height in pixels of the frame;
%
%       frameWidth is the height in pixels of the frame;
%
%       temperature in Celsius;
%
%       salinity (no units);
%
%       beadsDiameter is the diameter in microm of the microspheres.
%
%       binLength is the side of the squared bin in which the frames are
%           partitioned in microns
%
%       maxLag is the maximum lag, in seconds, of the squared displacements
%
%       microrheologyBias (optional) is an empirical bias that can be found
%           comparing the outputs from this function with theoretical
%           viscosities or measurements performed with a bulk-phase
%           viscometer. It will depend on the experimental setup. Default
%           is 0.
%
% output variables:
%
%       empiricalViscosity is an MxP matrix of the estimated dynamic
%           viscosities (in [kg/m/s]), where M an P are the number of
%           partitions of binLength side in the vertical and the horizontal
%           respectively.
%
%       theoreticalViscosity is the dynamic viscosity (in [kg/m/s]) of
%           seawater at the experimental temperature and salinity
%           calculated with SW_Viscosity.m, from the SEAWATER
%           THERMOPHYSICAL PROPERTIES LIBRARY (Sharqawy et al 2010).
%
%       N is an MxP matrix of the number of tracks registred in each
%           partition;
%
%       R is an MxP matrix of the R squared of the least square fitting to
%           least squares fitting the data to MSD(τ)=2dDτ.
%
%
% Dependencies: Econometrics Toolbox
%
% REFERENCES:
%   Guadayol, O., Mendonca, T., Segura-Noguera, M., Wright, A., Tassieri,
%    M., Humphries, S. (2020) Microrheology reveals microscale viscosity
%    gradients in planktonic systems. biorxv.
%   Sharqawy,M.H., Lienhard, J.H., Zubair S.M. (2010). Thermophysical properties of
%    seawater: a review of existing correlations and data, Desalination and
%    Water Treatment, 16 (2010) 354–380. doi:10.5004/dwt.2010.1079.
%
% EXAMPLE:
% [empiricalViscosity, theoreticalViscosity, N, R] =...
%     MVM_viscosities(tracks, scaling,...
%     frameRate, frameHeight, frameWidth,...
%     temperature, salinity,...
%     beadsDiameter, binLength, maxLag, microrheologyBias)
%
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

warning('off','stats:statrobustfit:IterationLimit')
if iscell(tracks)
    X = cellfun(@(x) (squeeze(x(:,1,:))*1e-6), tracks, 'UniformOutput', false);
    X = [X{:}];
    Y = cellfun(@(x) (squeeze(x(:,2,:))*1e-6), tracks, 'UniformOutput', false);
    Y = [Y{:}];
else
    X = squeeze(tracks(:,end-2,:))*1e-6;
    Y = squeeze(tracks(:,end-1,:))*1e-6;
end

if ~exist('microrheologyBias','var') || isempty('microrheologyBias')
    microrheologyBias = 0;
%     microrheologyBias = 2.337e-04;
end

maxLagFrames = floor(maxLag*frameRate);
time = (1:maxLagFrames)/frameRate;

%% Partitioning algorithm
XLAG = lagmatrix(X,1:maxLagFrames);
YLAG = lagmatrix(Y,1:maxLagFrames);
XLAG = reshape(XLAG,size(X,1),size(X,2),maxLagFrames);
YLAG = reshape(YLAG,size(Y,1),size(Y,2),maxLagFrames);
SD = (X-XLAG).^2 + (Y-YLAG).^2;
SD = reshape(SD,size(SD,1)*size(SD,2),size(SD,3));
X = X(:);
Y = Y(:);
isValid = ~any(isnan(SD),2);
SD = SD(isValid,:);
X = X(isValid,:);
Y = Y(isValid,:);

%% Diffusivity per xy bin
xbins = 0:binLength*1e-6:frameWidth/scaling*1e-6;
ybins = 0:binLength*1e-6:frameHeight/scaling*1e-6;
D = nan(length(xbins), length(ybins));
N = D;
R = D;
for xx = 1:length(xbins)-1
    for yy = 1:length(ybins)-1
        iParticles = find(X>xbins(xx) & X<=xbins(xx+1) & Y>ybins(yy) & Y<=ybins(yy+1));
        msd_z = mean(SD(iParticles,:),1);
        N(xx,yy) = length(iParticles);
        try
            if sum(~isnan(msd_z))==floor(maxLag*frameRate)
                fitted_model = fitlm(4*([0; time(:)]),[0; msd_z(:)],'RobustOpts','on','Intercept',false);
                D(xx,yy) = fitted_model.Coefficients.Estimate;
                R(xx,yy) = fitted_model.Rsquared.Adjusted;
            end
        end
    end
end

K = temperature+273.15; % temperature in Kelvins
k_b = 1.38e-23; % Boltzmann constant (J/K)
a = beadsDiameter*1e-6/2;
empiricalViscosity = k_b*K/6/pi./D./a;
empiricalViscosity = empiricalViscosity-microrheologyBias; % removes methodological bias
theoreticalViscosity = SW_Viscosity(temperature, 'C', salinity, 'ppt');