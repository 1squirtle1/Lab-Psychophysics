% getlastdir
%
%   usage: getlastdir(dir)
%      by: justin gardner
%    date: 2/04/99
% purpose: returns last directory in path
function retval = getlastdir(d)

if (isunix)
  slashes = findstr('/', d);
else
  slashes = findstr('\', d);
end

if (slashes(length(slashes)) == length(d))
  s = slashes(length(slashes)-1)+1;
  e = length(d)-1;
else
  s = slashes(length(slashes))+1;
  e = length(d);
end

retval = d(s:e);