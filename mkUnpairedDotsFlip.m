% function s = mkDotsC(stimSz, stimSz, dotDirection, dotSpeed, stimLoc
%                      stim.density, dotCoherence, dotContrast, dotLuminance,
%                      dotSz, densityStyle, sampleFactor, interpMethod)
%
% MKDOTSC makes a circular aperture drifting dot stimulus.
%
% Required arguments:
% stimSz            The dimensions of the entire stimulus, in [Y X T] coordinates;
% maskRad           The dimensions of the stimulus aperture (scalar).
% dotDirection      The direction of movement, in radians, with 0 = rightward.
%                   If dotDirection is just one number, the dots will move
%                   in that direction at all times. If dotDirection is a
%                   vector of the same length as the number of frames in
%                   the stimulus, the direction of movement in each frame
%                   will be specified by the corresponding element of
%                   dotDirection.
% dotSpeed          The speed of movement, in frames/second.
%                   If dotSpeed is just one number, the dots will move
%                   with that speed at all times. If dotSpeed is a
%                   vector of the same length as the number of frames in
%                   the stimulus, the speed of movement in each frame
%                   will be specified by the corresponding element of
%                   dotSpeed.
%
% Optional arguments:
% stimLoc           The location of the stimulus within the surround.
%                   Options are: N, NE, E, SE, S, SW, W, NW, DEFAULT = C
% stimDist          The distance from stimulus center to aperture center. 0
%                   if stimLoc == C.
% stim.density        The density of the dots, which can be from 0 to 1. DEFAULT = .1
% dotCoherence      The coherence of dot motion. DEFAULT = 1.
% dotContrast       The Michelson contrast ratio of the dots with the background.
% dotLuminance      The luminance of the dots, which can be 0 to 1.
%                   DEFAULT = 1 (i.e. 100% Luminance).
% dotSz         The radius of the dots. If dotSz < 0, the dots will
%                   be single pixels. If dotSz > 0, the dots will be
%                   Gaussian blobs with sigma = .5 * dotSz. DEFAULT = -1                    
% dotPlacementStyle The number of dots is calculated by multiplying the
%                   stim.density by the size of the image window. If the dots
%                   are placed randomly, some dots will overlap, so the
%                   actual dot density will be lower than stim.density. If,
%                   however, the dots are placed exactly so that they never
%                   overlap, all the dots will move in register. Both of
%                   these are problems that only occur when you use pixel
%                   dots rather than Gaussian dots. Your choices for
%                   dotPlacementStyle are 'random' (for the first problem)
%                   or 'exact' (for the second problem). DEFAULT = 'random'

function stim = mkUnpairedDotsFlip(stim)


% Parse arguments out of varargin and set default values
if ~isfield(stim, 'paired');     stim.paired = 1;                end
if ~isfield(stim, 'loc');        stim.loc = [0 0];               end
if ~isfield(stim, 'density');    stim.density = 0.01;            end
if ~isfield(stim, 'coherence');  stim.coherence = 1;             end
if ~isfield(stim, 'contrast');   stim.contrast = 0.89;           end
if ~isfield(stim, 'dotSz');      stim.dotSz = 2;                 end
if ~isfield(stim, 'luminance');  stim.luminance = -1;            end
if ~isfield(stim, 'pathlength'); stim.pathlength = 7;            end
if ~isfield(stim, 'color');      stim.color = [1 1 1];           end


if ~isfield(stim, 'size') || ~isfield(stim, 'direction') || ~isfield(stim, 'speed')
    error('stim needs fields size, direction, and speed');
end
if ~isfield(stim, 'maskRad')
    stim.maskRad = floor(min(stim.size(1:2))./2);
    stim.loc = [0 0];
end


% resize stim.size, dotDirection and dotSpeed if necessary
if length(stim.direction) == 1
    stim.direction = repmat(stim.direction, stim.size(3), 1);
end
if length(stim.speed) == 1
    stim.speed = repmat(stim.speed, stim.size(3), 1);
end

if length(stim.direction) ~= stim.size(3)
    error('If stim.direction is a vector, it must have the same number of entries as there are frames in the stimulus.');
end
if length(stim.speed) ~= stim.size(3)
    error('If stim.speed is a vector, it must have the same number of entries as there are frames in the stimulus.');
end
stim.density = stim.density;

% Make sure stim size is odd to center the circular aperture
if ~mod(stim.size(1),2) || ~mod(stim.size(2),2)
    error('Stimulus size must be odd to center circular aperture.');
end

%%%%%%%%%%%%%%%%% NOW, ON WITH THE CODE!!! %%%%%%%%%%%%%%%%%%%%%%%%%%

lum(1) = 4.45;
if (stim.luminance == -1)
    lum(2) = lum(1).*(stim.contrast+1)./(1-stim.contrast);
    lum = lum./80.1;
else
    lum(2) = stim.luminance;
end
loLum = 0;%lum(1);
hiLum = lum(2);
color = reshape(stim.color(1,:),1,1,3).*hiLum;

% Other dots
pathlength = stim.pathlength;         %in degrees
lifetime = pathlength./stim.speed(1); %in frames

% There is a buffer area around the image so that we don't have to worry
% about getting wraparounds exactly right.  This buffer is twice the size
% of a dot diameter.
bufferSize = max(stim.speed(:))+stim.dotSz+2*pathlength;
rad = floor(stim.dotSz./2);

% the 'frame' is the field across which the dots drift. The final ouput of
% this function will consist of 'snapshots' of this frame without buffer.
% We store the size of the frame and a coordinate system for it here:
frameSzX = stim.size(2) + 2.*bufferSize;
frameSzY = stim.size(1) + 2.*bufferSize;
frameSz = [frameSzY frameSzX];
[xFrame, yFrame] = meshgrid([1:frameSzX], [1:frameSzY]);
yFrame = flipud(yFrame);


% set up mask
cFrame = ceil(([frameSzY frameSzX]./2));
cStim = cFrame+[-1*stim.loc(1), stim.loc(2)];

% nDots is the number of coherently moving dots in the stimulus.
% nDotsNonCoherent is the number of noncoherently moving dots.
pos = ((xFrame-cStim(2)).^2+(yFrame-cStim(1)).^2)<=stim.maskRad^2;
pos = find(pos == 1);

nDots = max([1, round(stim.density.*length(pos))]);

nDotsNonCoherent = round((1-stim.coherence).*stim.density.*length(pos));

% Set the initial dot positions.
% dotPositions is a matrix of positions of the coherently moving dots in
% [y, x] coordinates; each row in dotPositions stores the position of one
% dot.
dotInds = randperm(length(pos), nDots);
dotInds = pos(dotInds);
[I,J] = ind2sub(frameSz, dotInds);
dotPositions = [I,J];

a = max([tan(stim.direction(1)-pi/2), eps]);

nKill = nDots./lifetime;        %dots to kill per frame
nKill = [0:stim.size(3)].*nKill;
nKill = diff(floor(nKill));
dot2kill = 1;


% s will store the output. After looping over each frame, we will trim away
% the buffer from s to obtain the final result.
r = zeros(frameSzY, frameSzX, stim.size(3).*3);
for t = 1:stim.size(3)
    
    % move the positions of all the dots
    dotVelocity = [sin(stim.direction(t)), cos(stim.direction(t))];
    dotVelocity = dotVelocity*stim.speed(t);
    dotPositions = dotPositions + repmat(dotVelocity, size(dotPositions, 1), 1);
    tmpDotPositions = round(dotPositions);
    dotInds = sub2ind(frameSz, tmpDotPositions(:,1), tmpDotPositions(:,2));
    
    % FOR RANDOMIZED NON-COHERENCE IB_EDIT
    if nDotsNonCoherent>0
        dotIndsNonCoherent = randperm(length(pos), nDotsNonCoherent);
        dotIndsNonCoherent = pos(dotIndsNonCoherent);
        overlap = ismember(dotIndsNonCoherent,dotInds);
        while sum(overlap)>0
            newInds = randperm(length(pos), sum(overlap));
            dotIndsNonCoherent(overlap) = pos(newInds);
            overlap = ismember(dotIndsNonCoherent,dotInds);
        end

        randDotIndex = randperm(nDots,nDotsNonCoherent);
        dotInds(randDotIndex) = dotIndsNonCoherent;
        [tmpDotPositions(randDotIndex,1), tmpDotPositions(randDotIndex,2)] = ind2sub(frameSz, dotInds(randDotIndex));
        dotPositions(randDotIndex,:) = tmpDotPositions(randDotIndex,:);
    end
    
    % FOR LIFETIME RESTRICTION
    if nKill(t)>0
        dotIndsRespawn = randperm(length(pos), nKill(t));
        dotIndsRespawn = pos(dotIndsRespawn);
        overlap = ismember(dotIndsRespawn,dotInds);
        while sum(overlap)>0
            newInds = randperm(length(pos), sum(overlap));
            dotIndsRespawn(overlap) = pos(newInds);
            overlap = ismember(dotIndsRespawn,dotInds);
        end

        kInds = [dot2kill:dot2kill+nKill(t)-1];
        kInds = 1+mod(kInds,nDots);
        dot2kill = dot2kill+nKill(t);
        dotInds(kInds) = dotIndsRespawn;
        [tmpDotPositions(kInds,1), tmpDotPositions(kInds,2)] = ind2sub(frameSz, dotInds(kInds));
        dotPositions(kInds,:) = tmpDotPositions(kInds,:);
    end
    
    dotInds = sub2ind(frameSz, tmpDotPositions(:,1), tmpDotPositions(:,2));
    
    % wrap around for all dots that have gone past the image borders
    w = ~ismember(dotInds,pos);
    if sum(w)>0
        tmpCStim = repmat(cStim,sum(w),1);
        tmpDots = tmpDotPositions(w,:)-tmpCStim;
        d = (tmpDots(:,2)+a*tmpDots(:,1))/(1+a^2);
        tmpDots = [2*d*a-tmpDots(:,1), 2*d-tmpDots(:,2)];
        tmpDotPositions(w,:) = round(tmpDots+cStim+repmat(dotVelocity, sum(w), 1));
        dotPositions(w,:) = tmpDotPositions(w,:);
    end
    
    % prepare a matrix of zeros for this frame
    thisFrame = zeros([frameSz,3])+loLum;
    
%     size(thisFrame)
    if stim.dotSz>1
        if mod(stim.dotSz,2)
            for j = 1:3
                for i = 1:size(tmpDotPositions,1)
                    thisFrame(tmpDotPositions(i,1)-rad:tmpDotPositions(i,1)+rad, ...
                        tmpDotPositions(i,2)-rad:tmpDotPositions(i,2)+rad,j) = color(j);
                end
            end
        else
            for j = 1:3
                for i = 1:size(tmpDotPositions,1)
                    thisFrame(tmpDotPositions(i,1)-rad+1:tmpDotPositions(i,1)+rad, ...
                        tmpDotPositions(i,2)-rad+1:tmpDotPositions(i,2)+rad,j) = color(j);
                end
            end
        end
    else
        for j = 1:3
            thisFrame(tmpDotPositions,j) = color(j);
        end
    end
%     size(thisFrame)
    thisFrame = thisFrame(1:frameSzY,1:frameSzX,1:3);
%     size(thisFrame)
    r(:,:,(t-1)*3+1:(t-1)*3+3) = flipud(thisFrame);%.*mask;

end
% Now trim away the buff
s = r(bufferSize+1:end-bufferSize, bufferSize+1:end-bufferSize, :);

% Set the luminances
% s = s.*hiLum;
% s(s==0) = loLum;
s(s>max(hiLum)) = max(hiLum);
stim.s = s;
% flipBook(s);
end
