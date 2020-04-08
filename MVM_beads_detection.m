function [Centroids, Boundaries, Area] = MVM_beads_detection(FRAMES, scaling, beadsDiameter)
%
% Detects moving fluorescent microspheres displaying Brownian motion in a
% grayscale video and records geometrical parameters of each particle
% in each frame. The particle detection algorithm is based on the algorithm
% develped by Crocker and Grier (1996). More details can be found at
% Guadayol et al 2020.
%
% input variables:
%
%       FRAMES is a MxNxP stack of uint8 images, where M and N are the
%           height and width of the frames respectively and P is the number
%           of frames.
%
%       scaling in pixels/microm;
%
%       beadsDiameter is the diameter in microm of the microspheres.
%
% output variables:
%
%       Centroids is a cell array of size (NumberofFramesX1), in which each
%          cell is a NX2 matrix, N being the number of particles detected within
%          the frame, and the 2 columns represent position x and y in pixels of
%          the centroids of the particles.
%
%       Boundaries is a cell array of size (NumberofFramesX1) in which each
%           cell is another array of 2-column matrices, one for each particle in
%           the frame, with x,y points of the boundary of the particle.
%
%       Area is a cell array of size (NumberofFramesX1), in which
%          each cell is an N long vector of areas (in pixels^2) for each
%          particle in the frame.
%
% DEPENDENCIES: image processing toolbox
%
% REFERENCES:
%    Crocker, J.C. and Grier, D.G. (1996). Methods of digital video
%     microscopy for colloidal studies", Journal of Colloid and
%     Interface Science 179, 298-310 (1996).
%     https://doi.org/10.1006/jcis.1996.0217
%
%
% ex.:
%     [Centroids, Boundaries, Area]  =...
%       MVM_beads_detection(FRAMES, scaling, beadsDiameter);
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
%  This files is part of micro-viscomaps, a collection of matlab scripts to
%  generate maps of viscosity from a video of microspheres displaying
%  Brownian motion, using Multiple Particle Tracking microrhology (MPTM).

highpass_window = ceil(beadsDiameter*scaling*2);
lowpass_window = 1;
normalize = @(x) x/sum(x);
gaussian_kernel = normalize(...
    exp(-((-ceil(5*lowpass_window):ceil(5*lowpass_window))/(2*lowpass_window)).^2));
boxcar_kernel = normalize(...
    ones(1,length(-round(highpass_window):round(highpass_window))));
minimum_area = round((pi*(beadsDiameter*scaling)^2));
maximum_area = round((pi*(beadsDiameter*2*scaling)^2));

background =  median(single(FRAMES),3);
parfor ff = 1:size(FRAMES,3)
    frame = uint8(single(FRAMES(:,:,ff))-background + min(background(:)));
    gconv = convn(frame,gaussian_kernel,'same');
    gconv = convn(gconv,gaussian_kernel','same');
    bconv = convn(frame,boxcar_kernel,'same');
    bconv = convn(bconv,boxcar_kernel','same');
    frame = gconv - bconv;
    frame = (uint8(frame*255/max(frame(:))));
    frame(:,1:highpass_window) = 0;
    frame(1:highpass_window,:) = 0;
    frame(size(frame,1)-highpass_window:end,:) = 0;
    frame(:,size(frame,2)-highpass_window:end) = 0;
    
    bw = imbinarize(frame);
    bw = imclearborder(bw); % removes particles on the edges
    bw = bwareaopen(bw, highpass_window); % removes objects smaller than n pixels
    bw = bwareafilt(bw,[minimum_area maximum_area]);
    [boundaries,L] = bwboundaries(bw,'noholes');
    Boundaries{ff} = cellfun(@uint16,boundaries, 'UniformOutput', false);
    stats = regionprops(L, 'Area', 'Centroid');
    Area{ff} = uint16([stats.Area]');
    Centroids{ff} = single(reshape([stats.Centroid],2,length(stats))');
end
