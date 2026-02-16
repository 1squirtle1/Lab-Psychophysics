function [ maskedStim ] = patchMask( varargin )
%Overlays a circular aperture over a stimulus.

patchSz = 'default';
patchDist = 'default';

stim = varargin{1};
patchLoc = varargin{2};
if nargin >= 3;  patchSz = varargin{3};         end
if nargin >= 4;  patchDist = varargin{4};         end

stimSz = size(stim);

if (strcmp(patchSz, 'default') || strcmp(patchSz, 'full')); patchSz = floor(min(stimSz(1:2))); end
if strcmp(patchDist, 'default'); patchDist = 9;                     end

maskRad = floor(patchSz./2);
cStim = floor(stimSz(1:2)./2)+1;
stimDistLeg = ceil(patchDist./sqrt(2));
switch patchLoc
    case 'C'
        cPatch = cStim;
    case 'N'
        cPatch = [cStim(1)-patchDist, cStim(2)];
    case 'NE'
        cPatch = [cStim(1)-stimDistLeg, cStim(2)+stimDistLeg];
    case 'E'
        cPatch = [cStim(1), cStim(2)+patchDist];
    case 'SE'
        cPatch = [cStim(1)+stimDistLeg, cStim(2)+stimDistLeg];
    case 'S'
        cPatch = [cStim(1)+patchDist, cStim(2)];
    case 'SW'
        cPatch = [cStim(1)+stimDistLeg, cStim(2)-stimDistLeg];
    case 'W'
        cPatch = [cStim(1), cStim(2)-patchDist];
    case 'NW'
        cPatch = [cStim(1)-stimDistLeg, cStim(2)-stimDistLeg];
    otherwise
        disp('invalid stimulus location');
end

mask = bsxfun(@plus, (((1:stimSz(2))-(cPatch(2))).^2),...
    (transpose(1:stimSz(1))-(cPatch(1))).^2) < maskRad.^2;
maskedStim = stim.*(repmat(mask, 1, 1, stimSz(3)));
maskedStim(maskedStim == 0) = min(stim(:));
end

