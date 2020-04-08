function tracks = MVM_beads_tracker(...
    Centroids, Boundaries, Area, scaling, frameRate, minimumTrackLength)

% Finds tracks that are longer than the minimumTrackLength and creates a
% MX3XP matrix with all the good tracks. M is the maximum track length in
% all tracks, and P is the number of good tracks detected. The code is
% based on trackbac toolbox.
% Columns are:
%
%       1: x position of the centroid in microns
%       2: y position of the centroid in microns
%       3: projected area of the microsphere in microns^2
%
% input variables:
%
%       Centroids is a cell array of size (NumberofFramesX1), in
%           which each cell is a NX2 matrix, N being the number of
%           particles detected within the frame, and the 2 columns
%           represent position x and y in pixels of the centroids of the
%           particles.
%
%       Area is a cell array of size (NumberofFramesX1), in which
%          each cell is an N long vector of areas (in pixels^2) for each
%          particle in the frame.
%
%       scaling is pixel/microm;
%
%       frameRate is the video frame rate in frames/s;
%
%       minimumTrackLength is the minimum length in seconds of the track to be
%           considered as such.
%
% REFERENCES:
%  Guadayol, O. 2016 trackbac v1.1 (Version v1.1). Zenodo.
%   http://doi.org/10.5281/zenodo.3452818
%  Guadayol, O., Thornton, K.L., Humphries, S. 2017. Cell morphology
%   governs directional control in swimming bacteria. Scientific Reports
%   7:2061. http://doi.org/10.1038/s41598-017-01565-y
%     
% 
%
% Ex.: tracks = MVM_bead_tracker(...
%      Centroids, Boundaries, Area, scaling, frameRate, minimumTrackLength)
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

warning('off','MATLAB:inpolygon:ModelingWorldLower')


%% tracks construction
nFrames = size(Centroids,2);
mm = 0;
dummy = Boundaries;
t{1} = [];
for ff = 1:nFrames-1 % frame
    for pp = 1:length(Boundaries{ff}) % particle
        
        if isempty(dummy{ff}{pp})
            continue
        end
        mm = mm+1;
        e = pp;
        t{mm} = [ff,pp]; % [frame, particle]
        for kk = ff:nFrames-1 % frame
            if isempty(Boundaries{kk+1})
                break
            end
            near = find(sqrt(sum((Centroids{kk+1}-repmat(Centroids{kk}(e,:),size(Centroids{kk+1},1),1)).^2,2))<=scaling*5); % finds particles in the following frame that are near the current particle (i,e, less than 5 microns away)
            E = false(length(near),1);
            for bb = 1:length(near)
                E(bb) = any([inpolygon(Boundaries{kk+1}{near(bb)}(:,1),Boundaries{kk+1}{near(bb)}(:,2),Boundaries{kk}{e}(:,1),Boundaries{kk}{e}(:,2));...
                    inpolygon(Boundaries{kk}{e}(:,1),Boundaries{kk}{e}(:,2),Boundaries{kk+1}{near(bb)}(:,1),Boundaries{kk+1}{near(bb)}(:,2))]);
            end
            e = near(E);
            
            if numel(e)==1
                t{mm} = [t{mm};kk+1,e]; % first column of each track is the
                % frame, second is the particle
                dummy{kk+1}{e} = []; % to remove particles that cross each
                % and to avoid retracking already
                % tracked particles
            elseif numel(e)>1
                for ll = 1:numel(e)
                    dummy{kk+1}{e(ll)} = []; % to remove particles that cross
                    % each and to avoid retracking
                    % already tracked particles
                end
                break
            else
                break
            end
        end
    end
end
t = cellfun(@uint16,t,'UniformOutput', 0);


%% tracks selection

[Cnrows,~] = cellfun(@size, t);
t(Cnrows<floor(frameRate*minimumTrackLength)) = [];
tracks = single(nan(nFrames, 3, size(t,2)));

warning('off', 'MATLAB:inpolygon:ModelingWorldLower')
for jj = 1:size(t,2)
    tracks(1:size(t{jj}),1:2,jj) = t{jj}; % frame number and particle ID
    for ii = 1:size(t{jj},1)
        tracks(ii,3:4,jj) = Centroids{t{jj}(ii,1)}(t{jj}(ii,2),:)/scaling; % position in microm
        tracks(ii,5,jj) = single(Area{t{jj}(ii,1)}(t{jj}(ii,2)))/scaling^2; % area in microm^2
    end
end

minimumArea = pi()*2*1;
dividedTracks = nan(size(tracks,1),size(tracks,2),0);
for jj = 1:size(tracks,3)
    kk = find(abs(diff(tracks(:,5,jj)))>minimumArea);% breaks the track when there is a change in area larger than 25%
    if ~isempty(kk)
        if kk(1)~=1
            kk = [1;kk];
        end
        k = find(~isnan(tracks(:,5,jj)),1,'last');
        if kk(end)~= k
            kk = [kk; k];
        end
        FF = flip(kk(2:end-1));
        for ii = 1:length(FF)
            if sum(~isnan(tracks(FF(ii):kk(end),3,jj)))>=floor(frameRate*minimumTrackLength)
                dividedTracks(:,:,end+1) = nan;
                dividedTracks(1:kk(end)-(FF(ii))+1,:,end) = tracks(FF(ii):kk(end),:,jj);
            end
            tracks(FF(ii):kk(end),:,jj) = nan;
        end
    end
end

tracks = cat(3,tracks, dividedTracks);
%% removes short tracks
goodTracks = squeeze(sum(~isnan(tracks(:,3,:))))>=floor(frameRate*minimumTrackLength);
tracks = tracks(:,:,goodTracks);