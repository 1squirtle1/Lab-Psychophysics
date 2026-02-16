% Paired / Unpaired Stimulus Generator
clear; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set desired parameters here, this should be the only place you edit
ppdeg = 25;          % pixels per degree
refresh = 100;       % screen refresh rate in Hz
patchDiam = 7.5;     % patch diameter in degrees
stimDur = 750;       % duration in ms
%stimDur = 2000;
ddensity = 2;        % dot density in dots/degree^2
dspeed = [7.5];        % dot speed in deg/sec %previously 7.5

s.contrast = [0.89 0.89];   %[CW;CCW]
s.dotSz = 2;                % dot size in px
s.colors = [1 0 1 ; 0 1 0]; % RGB color triplets for two components ([CW;CCW])

dirVec = [0:pi/6:11*pi/6];       % Directions in radians
%dirVec = [0];       % Directions in radians
% dsVec = [pi/6 pi/3 pi/2 2*pi/3 pi]; % Degrees of separation in radians
dsVec = [pi/2]; % Degrees of separation in radians
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These lines convert parameters from degrees,seconds to pixels,frames
s.size = round([patchDiam*ppdeg patchDiam*ppdeg stimDur*refresh/1000]);
s.size(1:2) = s.size(1:2) + mod(s.size(1:2),2) - 1;
s.speed = round(dspeed*ppdeg/refresh);
s.density = ddensity/ppdeg^2;
s.ppdeg = ppdeg;

% %thisDir = which('pairedDotGenerator.m');
% %thisDir = thisDir(1:end-length('pairedDotGenerator.m'));
% %addpath(thisDir);
% %componentDir = [thisDir,'\components'];
% %mkdir(componentDir);

if ~isequal(s.colors,ones(2,3))
    for j = dsVec
        dirName = ([thisDir '\' sprintf('%1iDS',j/pi*180)]);
        mkdir(dirName);
        cd(dirName);
        for i = dirVec
            s.direction = i;
            s.ds = j;
            P = mkPairedDotsFlip(s);

            r = s;
            r.color = s.colors(1,:);
            r.contrast = s.contrast(1);
            r.direction = i-s.ds./2;
            n1 = mkUnpairedDotsFlip(r);

            t = s;
            t.color = s.colors(2,:);
            t.contrast = s.contrast(2);
            t.direction = i+s.ds./2;
            n2 = mkUnpairedDotsFlip(t);

            N = s;
            N.s = mkBiStim(n1.s,n2.s);
            
            if(j==j(1))
                cd(componentDir);
                shMovie(n1,sprintf('componentColor%1iDIR',(i-j/2)/pi*180),r.color);
                cd(dirName);
            end
            shMovie(P,sprintf('pairedColor%1iDS%1iVA',j/pi*180,i/pi*180),r.color);
            shMovie(N,sprintf('unpairedColor%1iDS%1iVA',j/pi*180,i/pi*180),t.color);
        end
    end
end
