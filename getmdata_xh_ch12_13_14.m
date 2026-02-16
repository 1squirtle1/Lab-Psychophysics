% getmdata_xh_ch12_13.m
% to get response from channel 12 and 13

% based on getmdata_xh_111.m
%
% modified on 11/1/2011 by XH
%         by: justin
%       date: 03/01/00
%    purpose: get multiple data for trials with matching names and
%             put data into an array of data structures.
%      usage: d = getdata(trialtype,locutoff,hicutoff)
%       e.g.: d = getdata('b,type1:type2',200,450);
%             --This will create an array of d structures with
%             two elements, one for type1 trials and one
%             for type2 trials. It will look, for the 'b'
%             subtype, i.e. trials with names pk000811b.????
%             --If no trials match, getdata will print
%             out a list of all trialtypes found in directory
%             --if the type is just ':' it will get all trials
%             and sort them by trialname into an array of
%             structures.
%
%      notes: (1) most of what has been commented out is related to
%                 unnecessary processing
%             (2) added error checking if statements, the
%                 .choice field and the .targets field

function [d dTname lenTname] = getmdata_xh_ch12_13_14(trialsdir, trialtype,trialnums,cut,oth)

dTname = zeros(10, 30);
% check input arguments
if (nargin == 0)
  help getmdata_xh_ch12_13
  return
end

if(nargin<4 | isempty(cut))
  cut = true;
end

if(nargin<5 | isempty(oth))
  oth = true;
end

% see if we have a trial subcode
subcodepos = findstr(trialtype,',');
if (~isempty(subcodepos))
  subcode = trialtype(1:subcodepos(1)-1);
  if (length(trialtype) ~= subcodepos)
    trialtype = trialtype(subcodepos(1)+1:length(trialtype));
  else
    trialtype = [];
  end
else
  subcode = [];
end

% get match trialnames
matchtypes = [];
if (~strcmp(trialtype,':'))
  typepos = findstr(trialtype,':');
  typepos = [0 typepos length(trialtype)+1];
  numtrialtypes = length(typepos)-1;
  for i = 1:numtrialtypes
    % save all the trial names in our array
    matchtype = trialtype(typepos(i)+1:typepos(i+1)-1);
    % init paramaters for the trial type
    d(i).trialtype = matchtype;
    d(i).trialname = matchtype;
    d(i).n = 0;
    matchtypes = savenames(matchtypes,1,matchtype,1);
  end
else
  numtrialtypes = 0;
end

% calibration factors for traces
velcal = 10.8826;
poscal = 40;
dvelcal = velcal/2.84;

% only eight calibrations available. add more if there are more traces
cal = [poscal poscal velcal velcal dvelcal ...
       poscal dvelcal dvelcal velcal];
%use this for Mark's data
%cal = [poscal velcal poscal poscal poscal ...
%       poscal dvelcal dvelcal velcal];

% directory information
% trialsdir = 'C:\hx-ws\work\Physiology\Maestro_data\Puck\P32307d\';

% modified this on 5/4/2010
% cd C:\workspace\work\physiology\Maestro\MaestroData\Jennifer\41610\GI41610stt32\

% trialsdir = 'C:\workspace\work\physiology\Maestro\MaestroData\Jennifer\50510\GI50510deg60coh100-40\';
%  cd C:\workspace\work\physiology\Maestro\MaestroData\Jennifer\50510\GI50510deg60coh100-40

% cd C:\workspace\work\physiology\Maestro\MaestroData\Jennifer\42910\GI42910coh100-45
% trialsdir = 'C:\workspace\work\physiology\Maestro\MaestroData\Jennifer\42910\GI42910coh100-45\';

% cd C:\workspace\work\physiology\Maestro\MaestroData\Jennifer\50610\GI50610coh100-40
% trialsdir = 'C:\workspace\work\physiology\Maestro\MaestroData\Jennifer\50610\GI50610coh100-40\'

% cd C:\workspace\work\physiology\Maestro\MaestroData\Jennifer\50710\GI50710coh100-40
% trialsdir = 'C:\workspace\work\physiology\Maestro\MaestroData\Jennifer\50710\GI50710coh100-40\'

% cd C:\workspace\work\physiology\Maestro\MaestroData\Jennifer\51010\GI51010coh100-60
% trialsdir = 'C:\workspace\work\physiology\Maestro\MaestroData\Jennifer\51010\GI51010coh100-60\'

% cd C:\workspace\work\physiology\Maestro\MaestroData\Xin\0710\70910\E270910GT30c34\
% trialsdir = 'C:\workspace\work\physiology\Maestro\MaestroData\Xin\0710\70910\E270910GT30c34\'

% working directory must be the same as trialsdir

% [fname,trialsdir] = uigetfile('pick a maestro from the target directory')

experiment = [getlastdir(pwd) subcode];
if isunix
  numtrials = length(findstr(ls([trialsdir experiment '.????']),trialsdir));
else
  %doesn't seem to work if you're using unix dirs mounted over samba
  [numtrials, filelen] = size(ls([experiment '.*']));
end

%if trialnums not specified, just do 1:numtrials
if(nargin==2 | isempty(trialnums))
  trialnums = 1:numtrials;
end

trialnums
% pause;

% init some variables
names = []; counts = [];
saccadescut = 0;

% figure out how many spike traces there are
temp = readcxdata(strcat(trialsdir,experiment,'.0001'));

numunits = 0;
for i=1:13
  if ~isempty(temp.sortedSpikes{i})
    numunits = numunits+1;
  end
end

d.totaltrialn = length(trialnums);

for i = 1:length(trialnums)
  % load up file
  filename = strcat(experiment,sprintf('.%04d', trialnums(i)));
  mark = [];
  tdata = readcxdata(strcat(trialsdir,filename));
  tdata.key.flags
  
  chan12_13_14_temp = tdata.data;
%   d(i).chan6 = chan6_14_temp(1,:)';
%   d(i).chan7 = chan6_14_temp(2,:)';
%   d(i).chan8 = chan6_14_temp(3,:)';
%   d(i).chan9 = chan6_14_temp(4,:)';
%   d(i).chan10 = chan6_14_temp(5,:)';
%   d(i).chan11 = chan6_14_temp(6,:)';
  d(i).chan12 = chan12_13_14_temp(3,:)';
  d(i).chan13 = chan12_13_14_temp(4,:)';
  d(i).chan14 = chan12_13_14_temp(5,:)';
  
  astr = tdata.trialname
  dTname(1:length(astr),i) = astr;
  lenTname(i) = length(astr);
  
  tdata
  i
%   pause;
  
  if (~isempty(tdata.data))
    triallen = length(tdata.data(1,:));

    % find out how many traces there are
    numtraces = size(tdata.data,1);
    usetrial = 1;
  else
    triallen = 0;
    usetrial = 1;
    numtraces = 0;
  end

  % error checking for rejected trials
  %if ~isempty(tdata.mark1)
  %  if tdata.mark1(1) < 0
  %    usetrial = 0;
  %  end
  %end
  %switching to marked field
  %legacy feature below (old style of marking trials for discard)
  if((~isempty(tdata.marked) & tdata.marked==1) | ...
     (~isempty(tdata.mark1) & tdata.mark1(1)<0))
    usetrial = 0;
  end
  
  % if using SelbyFix2 and the computer chose the target
  if (tdata.key.flags == 108) | (tdata.key.flags == 92)
    usetrial = 0;
  end

  % remember trialnames
  [names counts] = savenames(names, counts, tdata.trialname, usetrial);

  % see if we are collecting all trial types
  if (strcmp(trialtype,':'))
    trialmatchnum = findname(matchtypes, tdata.trialname);
    % add trialtype if we haven't seen it before
    if (trialmatchnum == 0)
      matchtypes = savenames(matchtypes,1,tdata.trialname,1);
      numtrialtypes = numtrialtypes + 1;
      trialmatchnum = numtrialtypes;
      d(numtrialtypes).trialtype = tdata.trialname;
      d(numtrialtypes).trialname = tdata.trialname;
      d(numtrialtypes).n = 0;
    end
  else
    % just look up trials
    trialmatchnum = findname(matchtypes, tdata.trialname);
  end

%   % save data if it fits our search criteria
%   if (trialmatchnum & usetrial)
%     if (usetrial)
%       d(trialmatchnum).n = d(trialmatchnum).n + 1;
%       matchtrials = d(trialmatchnum).n;
%       d(trialmatchnum).trialnums(matchtrials) = trialnums(i);
% 
%       scaleddata = tdata.data ./ repmat(cal(1:numtraces)',1, triallen);
%       d(trialmatchnum).data{matchtrials} = scaleddata;
% 
%       % remember first mark1 and mark2
%       if (~isempty(tdata.mark1))
%         d(trialmatchnum).mark1(matchtrials,1:size(tdata.mark1,2)) = ...
%             tdata.mark1';
%       else
%         d(trialmatchnum).mark1(matchtrials,1) = 0;
%       end
% 
%       if (~isempty(tdata.mark2))
%         d(trialmatchnum).mark2(matchtrials,1:size(tdata.mark2,2)) = ...
%             tdata.mark2';
%       else
%         d(trialmatchnum).mark2(matchtrials,1) = 0;
%       end
% 
%       % if there are n events in tdata.other, the first n/2 will be the
%       % channel, and the second n/2 the time of the pulse
%       %d(trialmatchnum).others(matchtrials,k:) = tdata.other
%       if(oth)
%         %number of douts on this trial
%         numoth = size(tdata.other,1);
%         for k=1:numoth
%           d(trialmatchnum).othern(matchtrials,k) = tdata.other(k,1);
%           d(trialmatchnum).othert(matchtrials,k) = tdata.other(k,2);
%         end
%         if(isfield(d(trialmatchnum),'othern') & ...
%            numoth<size(d(trialmatchnum).othern,2))
%           d(trialmatchnum).othern(matchtrials,numoth+1:end) = NaN;
%           d(trialmatchnum).othert(matchtrials,numoth+1:end) = NaN;
%         end
%       end
% 
%       % process the .targets field.  we need to get rid of the data
%       % for targnum 0, which is RedLED1, as it is 0
%       index = find(tdata.targets.targnums ~= 0);
%       for j=1:length(index)
%         d(trialmatchnum).targets(matchtrials).hpos(j,:) = ...
%           tdata.targets.hpos(index(j),1:triallen);
%         d(trialmatchnum).targets(matchtrials).vpos(j,:) = ...
%           tdata.targets.vpos(index(j),1:triallen);
% 
%         d(trialmatchnum).targets(matchtrials).hvel(j,:) = ...
%             tdata.targets.hvel(j,:);
%         d(trialmatchnum).targets(matchtrials).vvel(j,:) = ...
%             tdata.targets.vvel(j,:);
%         
%         d(trialmatchnum).targets(matchtrials).patvelH(j,:) = ...
%             tdata.targets.patvelH(j,:);
%         d(trialmatchnum).targets(matchtrials).patvelV(j,:) = ...
%             tdata.targets.patvelV(j,:);
%       end
% 
% 
%     end
%   end
end

% unitColorMat = ['u','t','a','b','v'];
% if numunits
%   [unitNames dataTypes] = nex_info([getlastdir(pwd) '.' subcode '.nex']);
%   unitMarkers = find(~dataTypes);
%   for i=1:numtrialtypes
%     for m=1:numunits
%       d(i).color(m) = unitColorMat(str2num(unitNames(unitMarkers(m),4:6)));
%       d(i).unitname(m,:) = strcat(getlastdir(pwd),subcode,'.',unitNames(unitMarkers(m),1:7));
%     end
%   end
% end
% 
% % count the total number of found trials
% totaln = 0;
% for i = 1:numtrialtypes
%   totaln = totaln + d(i).n;
% end
% 
% % report trials found
% if (saccadescut == 0)
%   %disp('SACCADES NOT CUT');
% end
% if (totaln == 0)
%   if (isempty(trialtype))
%     disp('Found trials:');
%   else
%     disp(sprintf('Could not find trials of type: %s',trialtype))
%   end
%   printnames(names, counts);
%   d = [];
% else
%   if (nargout < 5)
%     disp(sprintf('Found %i trials', totaln));
%   end
% 
%   if 0;
%     for i = 1:numtrialtypes
%       if (usewaitbar & ~mod(i,usewaitbar)),waitbar(j/numtraces * (i/numtrialtypes));,end
%       %calibrate data if we have it.
%       if (d(i).n >= 1)
%         triallen = size(d(i).data,3);
%         %calibrate offsets using first 50 ms of traces
%         for j = 1:numtraces
%           offset = squeeze(nanmedian(squeeze(nanmedian(squeeze(d(i).data(:,j,1:50))))));
%           d(i).data(:,j,1:triallen) = d(i).data(:,j,1:triallen)-offset;
%         end
%       else
%         %put in one naned out trial, if there is no trials
%         d(i).n = 1;
%         d(i).data = nan*ones(1,8,2000);
%       end
%       d(i).median = squeeze(median(d(i).data,1));
%       d(i).mean = squeeze(mean(d(i).data,1));
%       d(i).std = squeeze(std(d(i).data,1));
%       %d(i).saclimits = [locutoff hicutoff];
%       %d(i).basename = experiment;
%     end
%   end
% end
% 
% %filter vel again
% %d = redovel(d);
% %interp saccades
% if(cut)
%   d = cutsac(d);
%   %cut vertical as well
%   d = cutsac(d,4);
% end
% 
%%%%%%%%%%%%%%%%%%%%%%%%

% remember unique name
%%%%%%%%%%%%%%%%%%%%%%%%
function [names, counts] = savenames(names, counts, trialname, usetrial)

namesize = size(names);
numnames = namesize(1);
match = 0;

for i = 1:numnames
  if (names(i,1) == length(trialname))
    if (strcmp(char(names(i,2:length(trialname)+1)),trialname))
      match = 1;
      counts(i,1) = counts(i,1) + 1;
      if (usetrial)
        counts(i,2) = counts(i,2) + 1;
      end
    end
  end
end

if (match == 0)
  names(numnames+1,1) = length(trialname);
  names(numnames+1,2:length(trialname)+1) = trialname;
  counts(numnames+1,1) = 1;
  if (usetrial)
    counts(numnames+1,2) = 1;
  else
    counts(numnames+1,2) = 0;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find if there is a matching name
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function match = findname(names, trialname)

namesize = size(names);
numnames = namesize(1);
match = 0;

for i = 1:numnames
  if (names(i,1) == length(trialname))
    if (strcmp(char(names(i,2:length(trialname)+1)),trialname))
      match = i;
    end
  end
end

%%%%%%%%%%%%%%%
% print names
%%%%%%%%%%%%%%%
function printnames(names, counts)

namesize = size(names);
numnames = namesize(1);

for i = 1:numnames
  disp(sprintf('%s\t: %i (%i)',char(names(i,2:names(i,1)+1)),counts(i,1),counts(i,2)));
end


