function shMovie(varargin)
stim = varargin{1};
name = varargin{2};
if nargin>2;    color = varargin;    end


if isfield(stim,'s')
    matrix = stim.s;
else
    matrix = stim;
end

if exist('color','var')
    matSz = stim.size;
    colorLength = length(color);
else
    matSz = size(matrix);
    colorLength = 1;
end

%v = VideoWriter([name '.avi']);%'Archival');
v = VideoWriter([name '.mp4'],'MPEG-4');
% v.Quality = 100;
% v.Path = path;
v.FrameRate=100;

open(v);
for i = 0:matSz(3)-1
   writeVideo(v,matrix(:,:,i*colorLength+1:i*colorLength+colorLength));
   %v.FrameRate=100;
   %v.FrameRate
end
close(v);