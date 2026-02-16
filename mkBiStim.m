function s = mkBiStim( varargin )
% mkBiDots Makes a bidirectional drifting dot stimulus from
%          two component sets of dot field parameters.

% Parse arguments
if nargin >= 2
    s1 = varargin{1};
    s2 = varargin{2};
    loLum = max([min(s1(:)), min(s2(:))]);
    % If arrays given, just append
    s = (s1+s2)-loLum;
else
    disp('two stimuli must be entered');
    return
end

if nargin ==3
    contrast = varargin{3};
elseif nargin > 3
    disp('too many arguments');
    return
end

lum = max([max(s1(:)), max(s2(:))]);
bck = min([min(s1(:)), min(s2(:))]);

s(s>lum) = lum;
s(s<bck) = bck;

end