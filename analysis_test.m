% Analyze_SpeedPerception_Transparency.m
% 8/29/2017

% Maestro filename: HuangLab_Psycho_JX.cxe
% Maestro trialset: 1or2speed_2T4T40coh
 
clear all;
close all;
  
% cd C:\LabFolderNew\Xin\SpeedPsychophysics\MaestroData
%cd C:\Workspace\Workspace\Work\Projects\Psychophysics\Speed_Perception\Emily\MaestroData\SpeedSegmentationData\Practice\OR
cd P:\labFolderNew\Bikalpa\matlab_base\project_human_psych\speed_psychophysics\
[fname,pname] = uigetfile('*.*','pick a maestro data file');
 
trialsdir = pname;
cd(trialsdir);
[d dTname lenTname] = getmdata_xh_ch12_13(trialsdir, 'b,');
 
 
total_trial_n = d(1).totaltrialn;
 
% The data length for each trial may be slightly different
 
chan12_len = zeros(total_trial_n,1);%edit ib
chan13_len = zeros(total_trial_n,1);
 
for i = 1:total_trial_n,
    temp12 = d(i).chan12;
    chan12_len(i) = length(temp12);
    
    temp13 = d(i).chan13;
    chan13_len(i) = length(temp13);
           
end;    
 
minCh12_len = min(chan12_len);
minCh13_len = min(chan13_len);
 
START_TIME = 100;  % examine respone after 600ms
 
chan12 = zeros(minCh12_len'-START_TIME, total_trial_n);
chan13 = zeros(minCh13_len'-START_TIME, total_trial_n);
 
 
% parse trialname
 
total_cdnum = 44;
 
 
% d(1).totaltrialn
for i = 1:total_trial_n,
    temp12 = d(i).chan12;
    chan12(:,i) = temp12(START_TIME+1:minCh12_len);
    
    temp13 = d(i).chan13;
    chan13(:,i) = temp13(START_TIME+1:minCh13_len);
end;    
 
figure(1);
clf(1);
subplot(2,2,1);
for i = 1:total_trial_n,
  plot(chan12(:,i), 'b');
  title('chan12-single');
  hold on;
end;  
 
subplot(2,2,2);
for i = 1:total_trial_n,
  plot(chan13(:,i), 'r');
  title('chan13-double');
  hold on;
end;  
 
% superimpose Chan12 and Chan13 for each trial
 
% logic OR of chan12 and 13
chan12OR13 = chan12 + chan13;
 
% difference chan12OR13, finding the rising phase
diff_chan12OR13 = diff(chan12OR13);
 
% difference of chan12, finding the rising phase
diff_chan12 = diff(chan12);
 
% difference of chan13, finding the rising phase
diff_chan13 = diff(chan13);
 
 
% set threshold as 800, if voltage difference > 800, consider an event
% the first event has to be after 1000 ms 
% the interval between the first and second intervals need > 200 ms
 
% find location of the rising phase
 
THRESHOLD = 600;
BASELINE = 200;
% FIRST_DELAY = 300;    % ms %changed from 500 to 300ms on 12/28/14
FIRST_DELAY = 300;    % ms %changed from 300 to 0 ms on 7/25/16
% INTERVAL = 200;       % ms
 
Evnt_1 = zeros(total_trial_n, 1);
 
missed_trial_ct = 0;
for i = 1:total_trial_n,
 
    % because it is from difference, +1 (move one step next) the get the
    % real peak
    if (max(diff_chan12OR13(:,i)) > THRESHOLD),
        first_rise_loc_chan12OR13(i) = min(find(diff_chan12OR13(:,i) > THRESHOLD))+1;
    else
        first_rise_loc_chan12OR13(i) = -999;
    end
    
           
    % make sure the timing of the first event and the interval between 
    % the first and second events meet criteria
    
    if first_rise_loc_chan12OR13(i) < FIRST_DELAY,
        sprintf('First EVENT OCCURRED TOO SOON');
        Evnt_1(i) = -3;
        Evnt_2(i) = -3;
        missed_trial_ct = missed_trial_ct + 1;
        missed_trial(missed_trial_ct) = i; 
    end;    
                
    % meet the criteria    
    
    if Evnt_1(i) == 0;
    
        if     (chan12(first_rise_loc_chan12OR13(i),i) > THRESHOLD) && ...
             (chan13(first_rise_loc_chan12OR13(i),i) < BASELINE)
            Evnt_1(i) = 12;
        elseif (chan13(first_rise_loc_chan12OR13(i),i) > THRESHOLD) && ...
             (chan12(first_rise_loc_chan12OR13(i),i) < BASELINE)
            Evnt_1(i) = 13;  
        end;    
    
    end;
 
             
               
  figure(2);
  clf(2);
  
  subplot(2,2,1);
  plot(chan12(:,i), 'b');
  hold on;
  plot(chan13(:,i), 'r');
  events = sprintf('channel %d',Evnt_1(i));
  title(events);
  
%   subplot(2,2,2);
%   plot(chan12OR13(:,i), 'k');
%   hold on;
%   plot(first_rise_loc_chan12OR13(i), chan12OR13(first_rise_loc_chan12OR13(i),i), 'go');
    
  subplot(2,2,4);
  plot(diff_chan12OR13(:,i), 'k');
         
  
end;    
 
[Evnt_1]
 
for i = 1:total_trial_n;
    ptname(i, 1:lenTname(i)) = char(dTname(1:lenTname(i),i)');
    
    switch (ptname(i,1:lenTname(i)))
        case '1'
            cond_num(i) = 1;
        case '110'
            cond_num(i) = 2;
        case '2'
            cond_num(i) = 3;        
        case '23'
            cond_num(i) = 4;  
                        
        case '3'
            cond_num(i) = 5;  
        case '31'
            cond_num(i) = 6;    
        case '4'
            cond_num(i) = 7;             
        case '41'
            cond_num(i) = 8; 
                        
        case '5'
            cond_num(i) = 9;             
        case '51'
            cond_num(i) = 10;              
        case '6'
            cond_num(i) = 11;            
        case '61'
            cond_num(i) = 12;
                        
        case '7'
            cond_num(i) = 13;            
        case '71'
            cond_num(i) = 14;              
        case '8'
            cond_num(i) = 15;            
        case '81'
            cond_num(i) = 16;
                        
        case '9'
            cond_num(i) = 17;            
        case '91'
            cond_num(i) = 18;              
        case '10'
            cond_num(i) = 19;            
        case '101'
            cond_num(i) = 20;
                        
        case '11'
            cond_num(i) = 21;             
        case '111'
            cond_num(i) = 22; 

        case '12'
            cond_num(i) = 23;
        case '121'
            cond_num(i) = 24;

        case '13'
            cond_num(i) = 25;        
        case '131'
            cond_num(i) = 26;  
                        
        case '14'
            cond_num(i) = 27;  
        case '141'
            cond_num(i) = 28;    
        case '15'
            cond_num(i) = 29;             
        case '151'
            cond_num(i) = 30; 
                        
        case '16'
            cond_num(i) = 31;             
        case '161'
            cond_num(i) = 32;              
        case '17'
            cond_num(i) = 33;            
        case '171'
            cond_num(i) = 34;
                        
        case '18'
            cond_num(i) = 35;            
        case '181'
            cond_num(i) = 36;              
        case '19'
            cond_num(i) = 37;            
        case '191'
            cond_num(i) = 38;
                        
        case '20'
            cond_num(i) = 39;            
        case '201'
            cond_num(i) = 40;              
        case '21'
            cond_num(i) = 41;            
        case '211'
            cond_num(i) = 42;
                        
        case '22'
            cond_num(i) = 43;             
        case '221'
            cond_num(i) = 44;

            
        
                        
        otherwise   
        display('Incorrect trialname');  
    end
    
end;

 
for mtct = 1:missed_trial_ct
    missed_cond(mtct) = cond_num(missed_trial(mtct))       
end;     
    
strength = zeros(1,total_trial_n);
 
figure(3);
clf(3);
plot(1:total_trial_n,cond_num,'ro-');
xlabel('trial number');
ylabel('Condition #');
 
%edit ib
listResp(:,1) = 1:total_trial_n;
listResp(:,2) = cond_num;
listResp(:,3) = Evnt_1;

Resp = zeros(1,total_cdnum); 
Count = zeros(1,total_cdnum); 

for i = 1:total_trial_n
    if listResp(i,3) == 13
       listResp(i, 4) = 1;
       Resp(cond_num(i)) = Resp(cond_num(i)) + listResp(i, 4);
       Count(cond_num(i)) = Count(cond_num(i)) + 1;
        
    elseif listResp(i,3) == 12
       listResp(i, 4) = 0;
       Resp(cond_num(i)) = Resp(cond_num(i)) + listResp(i, 4);
       Count(cond_num(i)) = Count(cond_num(i)) + 1;
        
    elseif listResp(i,3) == -3
        listResp(i, 4) = NaN;
    else disp 'error';           
    end
    
    
end

listResp;

Ratio = Resp./Count

% [condition number, total times reporting chan 12 , trial count, ratio]
[[1:total_cdnum]' Resp' Count' Ratio']

% x2 Speed difference
T2_condnum     = [12:17];
T2_equal       = [1:6];

% x4 Speed difference
T4_condnum     = [18:22];
T4_equal       = [7:11];

%%%%%%%%%%%%%%%%%%% For x2 speed difference %%%%%%%%%%%%%%%%%%% %%%%%%%%%%
% 2T, 40% Coh
T2_bothx2 = [T2_condnum' Resp(T2_condnum)' Count(T2_condnum)' Ratio(T2_condnum)']

% 2T, 40% Coh Equal
T2_singlex2 = [T2_equal' Resp(T2_equal)' Count(T2_equal)' Ratio(T2_equal)']


%%%%%%%%%%%%%%%%%%% For x4 speed difference %%%%%%%%%%%%%%%%%%% %%%%%%%%%%
% 4T, 40% Coh
T4_bothx4 = [T4_condnum' Resp(T4_condnum)' Count(T4_condnum)' Ratio(T4_condnum)']

% 4T, 40% Coh Equal
T4_singlex4 = [T4_equal' Resp(T4_equal)' Count(T4_equal)' Ratio(T4_equal)']

figure; 
subplot(1,2,1)
x = [1 2 3 4 5];
vals = [Ratio(T4_condnum); Ratio(T4_equal)];
b = bar(x,vals);
title('x4 HR & FA')

subplot(1,2,2)
x = [1 2 3 4 5 6];
vals = [Ratio(T2_condnum); Ratio(T2_equal)];
b = bar(x,vals);
title('x2')



% if hit(i) >= 101/102 || fa(i) <= 1/102 || hit(i) <= 1/102 || fa(i) >= 101/102
%     dprime(i)= norminv(((100 * hit(i)) + 1)/102) - norminv(((100 * fa(i)) + 1)/102);
% % elseif fa(i) > 1/102 || hit(i) < 101/102
% else
%     dprime(i)= norminv(hit(i)) - norminv(fa(i));
% end


