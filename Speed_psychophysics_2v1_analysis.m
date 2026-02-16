% Analyze_2vs1_Speed_Segmentation.m

% made by Steven on 2/25/15
% call getmdata_xh_ch12_13

% Maestro filename: 
% Maestro trialset: 

clear all;
close all;


% cd C:\LabFolder\MaestroData\Chris\5vs2_spd5_stair\0815\081915
% cd Y:\MaestroData\Xin
cd P:\labFolderNew\Emily\'Human Psychophysics'
[fname,pname] = uigetfile('*.*','pick a maestro data file');

trialsdir = pname;
cd(trialsdir);
[d dTname lenTname] = getmdata_xh_ch12_13(trialsdir, 'a,:');


total_trial_n = d(1).totaltrialn;

% The data length for each trial may be slightly different

chan12_len = zeros(total_trial_n);
chan13_len = zeros(total_trial_n);

for i = 1:total_trial_n,
    temp12 = d(i).chan12;
    chan12_len(i) = length(temp12);
    
    temp13 = d(i).chan13;
    chan13_len(i) = length(temp13);
           
end;    

minCh12_len = min(chan12_len);
minCh13_len = min(chan13_len);

START_TIME = 100;  % examine respone after 600ms

% chan12 = zeros(minCh12_len'-START_TIME, total_trial_n);
% chan13 = zeros(minCh13_len'-START_TIME, total_trial_n);


% parse trialname

total_cdnum = 56;


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

THRESHOLD = 800;
BASELINE = 200;
FIRST_DELAY = 300;    % ms %changed from 500 to 300ms on 12/28/14
INTERVAL = 200;       % ms

Evnt_1 = zeros(total_trial_n, 1);


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
        Evnt_1(i) = -1;
        Evnt_2(i) = -1;
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
% cond_num = [];
% TwoDirTrials = cond_num(1:2:56);
TwoDirTrials = (1:2:56);

cond_num = zeros(1,total_trial_n);




for i = 1:total_trial_n;
    ptname(i, 1:lenTname(i)) = char(dTname(1:lenTname(i),i)');
    switch (ptname(i,1:lenTname(i)))
        case '1_2dir_2.5deg'
            cond_num(i) = 1;
        case '1_1dir_2.5deg'
            cond_num(i) = 2;
        case '2_2dir_5deg'
            cond_num(i) = 3;
        case '2_1dir_5deg'
            cond_num(i) = 4;
        case '3_2dir_7.5deg'
            cond_num(i) = 5;
        case '3_1dir_7.5deg'
            cond_num(i) = 6;
        case '4_2dir_10deg'
            cond_num(i) = 7;
        case '4_1dir_10deg'
            cond_num(i) = 8;
        case '5_2dir_12.5deg'
            cond_num(i) = 9;
        case '5_1dir_12.5deg'
            cond_num(i) = 10;            
        case '6_2dir_15deg'
            cond_num(i) = 11;
        case '6_1dir_15deg'
            cond_num(i) = 12;
        case '7_2dir_17.5deg'
            cond_num(i) = 13;
        case '7_1dir_17.5deg'
            cond_num(i) = 14;
        case '8_2dir_20deg'
            cond_num(i) = 15;
        case '8_1dir_20deg'
            cond_num(i) = 16;
        case '9_2dir_22.5deg'
            cond_num(i) = 17;
        case '9_1dir_22.5deg'
            cond_num(i) = 18;
        case '10_2dir_25deg'
            cond_num(i) = 19;
        case '10_1dir_25deg'
            cond_num(i) = 20;    
        case '11_2dir_27.5deg'
            cond_num(i) = 21;
        case '11_1dir_27.5deg'
            cond_num(i) = 22;
        case '12_2dir_30deg'
            cond_num(i) = 23;
        case '12_1dir_30deg'
            cond_num(i) = 24;
        case '13_2dir_32.5deg'
            cond_num(i) = 25;
        case '13_1dir_32.5deg'
            cond_num(i) = 26;
        case '14_2dir_35deg'
            cond_num(i) = 27;
        case '14_1dir_35deg'
            cond_num(i) = 28;
        case '15_2dir_37.5deg'
            cond_num(i) = 29;
        case '15_1dir_37.5deg'
            cond_num(i) = 30;            
        case '16_2dir_40deg'
            cond_num(i) = 31;
        case '16_1dir_40deg'
            cond_num(i) = 32;
        case '17_2dir_42.5deg'
            cond_num(i) = 33;
        case '17_1dir_42.5deg'
            cond_num(i) = 34;
        case '18_2dir_45deg'
            cond_num(i) = 35;
        case '18_1dir_45deg'
            cond_num(i) = 36;
        case '19_2dir_47.5deg'
            cond_num(i) = 37;
        case '19_1dir_47.5deg'
            cond_num(i) = 38;
        case '20_2dir_50deg'
            cond_num(i) = 39;
        case '20_1dir_50deg'
            cond_num(i) = 40;
        case '21_2dir_52.5deg'
            cond_num(i) = 41;
        case '21_1dir_52.5deg'
            cond_num(i) = 42;    
        case '22_2dir_55deg'
            cond_num(i) = 43;
        case '22_1dir_55deg'
            cond_num(i) = 44;
        case '23_2dir_57.5deg'
            cond_num(i) = 45;
        case '23_1dir_57.5deg'
            cond_num(i) = 46;
        case '24_2dir_60deg'
            cond_num(i) = 47;
        case '24_1dir_60deg'
            cond_num(i) = 48;
        case '25_2dir_62.5deg'
            cond_num(i) = 49;
        case '25_1dir_62.5deg'
            cond_num(i) = 50;
        case '26_2dir_65deg'
            cond_num(i) = 51;
        case '26_1dir_65deg'
            cond_num(i) = 52;
        case '27_2dir_67.5deg'
            cond_num(i) = 53;
        case '27_1dir_67.5deg'
            cond_num(i) = 54;
        case '28_2dir_70deg'
            cond_num(i) = 55;
        case '28_1dir_70deg'
            cond_num(i) = 56;    

            
            
        otherwise   
        display('Incorrect trialname');  
    end
    
end;



% Creating a variable name for pulling out the strength from trialname

for i = 1:total_trial_n,
stname(i,1:2) = char(dTname(1:2,i))';
end


%finding the strength and plotting the results

strength = zeros(1,length(Test_Trials));

for i = 1:length(Test_Trials),
 if stname(i,2) == '_',
    strength(i) = str2num(stname(i,1));
 else strength(i) = str2num(stname(i,1:2));
 end
end


% converting the strength from Maestro into actual % speed increase

act_strength = zeros(1,total_trial_n);
% for i = 1:total_trial_n,
for i = 1:Test_Trials,
 act_strength(i) = strength(i) * 2.5 ;
end;

% Pull out the 2 direction trials only
Two_Dir = rem(cond_num,2) == 1;
Test_Trials = cond_num(Two_Dir);

figure(3);
clf(3);
plot(1:total_trial_n,act_strength,'ro-');
xlabel('trial number');
ylabel('Direction Separation (deg)');
axis([0 total_trial_n+5 min(act_strength - 3) max(act_strength + 3)]);

%%finding reversals and averaging the last 4

reversalpos = 1;
reversalcount = 0;
reversals = zeros(7,1);

for trials = 2:total_trial_n,
if (act_strength(trials-1) > act_strength(reversalpos)) ...
     && act_strength(trials) < act_strength(trials-1),
        reversalpos = trials-1;
        reversalcount = reversalcount + 1;
        reversals(reversalcount) = act_strength(trials-1);
    elseif  (act_strength(trials-1) < act_strength(reversalpos)) ...
     && (act_strength(trials) > act_strength(trials-1)),  
        reversalpos = trials-1;
        reversalcount = reversalcount + 1;
        reversals(reversalcount) = act_strength(trials-1);
end;
 
end;

% last input where the paradigm stopped is the last reversal

reversals(reversalcount+1) = act_strength(length(cond_num));

display('reversals');
reversals'
reversalpos;

n = length(find(reversals ~= -99));
if (n >= 7)
    display ('staircase completed');
    reversals(find(reversals ~= -99))
    display('Mean of the last 4 reversals:');
    meanRev = mean(reversals(n-3:n))
    display('Total trial numbers') 
    trialNum = total_trial_n
end;



