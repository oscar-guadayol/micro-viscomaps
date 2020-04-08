%  'micro-viscomaps' is a set of matlab scripts that uses a Multiple
%  Particle Tracking microrhology (MPTM) approach to generate maps of
%  viscosity from a video of microspheres displaying Brownian motion.
%
%  Copyright (C) 2020,  Oscar Guadayol <oscar_at_guadayol.cat>
%
% LICENSE:
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
%  with this program. If not, see <http://www.gnu.org/licenses/>.
%
% FUNCTIONS:
%
%   MVM_beads_detection.m: detects and characterizes fluorescent beads displaying
%       Brownian motion in a grayscale video.
%   MVM_beadss_tracker.m: tracks the beads in the video using the 'Boundaries'
%       cell array from beads_detection.
%   MVM_viscosities.m: partitions the frame into squares and computes the
%       viscosity within each square.
%   MVM_viscosity_map.m: creates and exports maps of absolute and relative
%       dynamic viscosities.
%
%   
% DEPENDENCIES:
%
%   Econometrics Toolbox,
%   Image Processing Toolbox,
%   SW_Viscosity.m (from SEAWATER THERMOPHYSICAL PROPERTIES
%       LIBRARY, Sharqawy et al 2010)
%
% REFERENCES:
%
%	Crocker, J.C. and Grier, D.G. (1996). Methods of digital video
%     microscopy for colloidal studies", Journal of Colloid and Interface
%     Science 179, 298-310 (1996). https://doi.org/10.1006/jcis.1996.0217.
%   Guadayol, O. 2016 trackbac v1.1 (Version v1.1). Zenodo.
%     http://doi.org/10.5281/zenodo.3452818.
%   Guadayol, O., Thornton, K.L., Humphries, S. 2017. Cell morphology
%     governs directional control in swimming bacteria. Scientific Reports
%     7:2061. http://doi.org/10.1038/s41598-017-01565-y.
%   Guadayol, O.,Mendonca, T., Segura-Noguera, M., Wright, A., Tassieri,
%     M., Humphries, S. (2020) Microrheology reveals microscale viscosity
%     gradients in planktonic systems. biorxv.
%   Sharqawy,M.H., Lienhard, J.H., Zubair S.M. (2010). Thermophysical
%     properties of seawater: a review of existing correlations and data,
%     Desalination and Water Treatment, 16 (2010) 354–380.
%     doi:10.5004/dwt.2010.1079.







