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
[d dTname lenTname] = getmdata_xh_ch12_13_14(trialsdir, 'b,');
 
 
total_trial_n = d(1).totaltrialn;
 
% The data length for each trial may be slightly different
 
chan12_len = zeros(total_trial_n,1);%edit ib
chan13_len = zeros(total_trial_n,1);
chan14_len = zeros(total_trial_n,1);
 
for i = 1:total_trial_n
    temp12 = d(i).chan12;
    chan12_len(i) = length(temp12);
    
    temp13 = d(i).chan13;
    chan13_len(i) = length(temp13);

    temp14 = d(i).chan14;
    chan14_len(i) = length(temp14);
           
end    
 
minCh12_len = min(chan12_len);
minCh13_len = min(chan13_len);
minCh14_len = min(chan14_len);
 
START_TIME = 1500;  % examine respone after 600ms
% Resp_window_len=1500;%Response window length=1500 ms

chan12 = zeros(minCh12_len-START_TIME, total_trial_n);
chan13 = zeros(minCh13_len-START_TIME, total_trial_n);
chan14 = zeros(minCh14_len-START_TIME, total_trial_n);
 
 
% parse trialname
 
total_cdnum = 40;
first_rise_loc_diff_chan121314=ones(1,total_trial_n);

 
% d(1).totaltrialn
for i = 1:total_trial_n
    temp12 = d(i).chan12;
    chan12(:,i) = temp12(START_TIME+1:minCh12_len);
    
    temp13 = d(i).chan13;
    chan13(:,i) = temp13(START_TIME+1:minCh13_len);

    temp14 = d(i).chan14;
    chan14(:,i) = temp14(START_TIME+1:minCh14_len);
end    
 
figure(1);
clf(1);
subplot(1,3,1);
for i = 1:total_trial_n
  plot(chan12(:,i), 'b');
  title('chan12-first_two');
  hold on;
end  
 
subplot(1,3,2);
for i = 1:total_trial_n
  plot(chan13(:,i), 'r');
  title('chan13-second_two');
  hold on;
end

subplot(1,3,3);
for i = 1:total_trial_n
  plot(chan14(:,i), 'c');
  title('chan14-notsure');
  hold on;
end
 
% superimpose Chan12 and Chan13 for each trial
 
% % logic OR of chan12 and 13
% chan12OR13 = chan12 + chan13;
%  
% % difference chan12OR13, finding the rising phase
% diff_chan12OR13 = diff(chan12OR13);
%  
% % difference of chan12, finding the rising phase
% diff_chan12 = diff(chan12);
%  
% % difference of chan13, finding the rising phase
% diff_chan13 = diff(chan13);
 
chanOR121314 = chan12 + chan13 + chan14;

% difference chan6-14, finding the rising phase
diff_chan121314= diff(chanOR121314);

% difference of chan6-14, finding the rising phase

diff_chan12 = diff(chan12);
diff_chan13 = diff(chan13);
diff_chan14 = diff(chan14);


 
% set threshold as 800, if voltage difference > 800, consider an event
% the first event has to be after 1000 ms 
% the interval between the first and second intervals need > 200 ms
 
% find location of the rising phase
 
THRESHOLD = 600;
BASELINE = 200;
FIRST_DELAY = 300;    % ms %changed from 500 to 300ms on 12/28/14
% FIRST_DELAY = 0;    % ms %changed from 300 to 0 ms on 7/25/16
% INTERVAL = 200;       % ms
 
Evnt_1 = zeros(total_trial_n, 1);
 
missed_trial_ct = 0;
for i = 1:total_trial_n
 
    % because it is from difference, +1 (move one step next) the get the
    % real peak
% % %     if (max(diff_chan12OR13(:,i)) > THRESHOLD),
% % %         first_rise_loc_chan12OR13(i) = min(find(diff_chan12OR13(:,i) > THRESHOLD))+1;
% % %     else
% % %         first_rise_loc_chan12OR13(i) = -999;
% % %     end
    
           
    % make sure the timing of the first event and the interval between 
    % the first and second events meet criteria
    
% % %     if first_rise_loc_chan12OR13(i) < FIRST_DELAY,
% % %         sprintf('First EVENT OCCURRED TOO SOON');
% % %         Evnt_1(i) = -3;
% % %         Evnt_2(i) = -3;
% % %         missed_trial_ct = missed_trial_ct + 1;
% % %         missed_trial(missed_trial_ct) = i; 
% % %     end;    
                
    % meet the criteria    
    
% % %     if Evnt_1(i) == 0;
% % %     
% % %         if     (chan12(first_rise_loc_chan12OR13(i),i) > THRESHOLD) && ...
% % %              (chan13(first_rise_loc_chan12OR13(i),i) < BASELINE)
% % %             Evnt_1(i) = 12;
% % %         elseif (chan13(first_rise_loc_chan12OR13(i),i) > THRESHOLD) && ...
% % %              (chan12(first_rise_loc_chan12OR13(i),i) < BASELINE)
% % %             Evnt_1(i) = 13;  
% % %         end;    
% % %     
% % %     end;
 
%       if (max(diff_chan121314(:,i)) > THRESHOLD )
%             first_rise_loc_diff_chan121314(i) = min(find(diff_chan121314(:,i) > THRESHOLD))+1;
%          if (chanOR121314(first_rise_loc_diff_chan121314(i),i)<0 && max(diff_chan14(:,i)> THRESHOLD))
%             first_rise_loc_diff_chan121314(i) = min(find(diff_chan14(:,i)> THRESHOLD))+1; 
%          end
%          
%      
%      else
%          
%         Evnt_1(i)=0;       
%       end    

      if (max(diff_chan121314(:,i)) > THRESHOLD)
        first_rise_loc_diff_chan121314(i) = find(diff_chan121314(:,i) > THRESHOLD, 1 )+1;
    else
        first_rise_loc_diff_chan121314(i) = -999;
    end
    
           
    % make sure the timing of the first event and the interval between 
    % the first and second events meet criteria
    
    if first_rise_loc_diff_chan121314(i) < FIRST_DELAY,
        sprintf('First EVENT OCCURRED TOO SOON');
        Evnt_1(i) = -3;
        Evnt_2(i) = -3;
        missed_trial_ct = missed_trial_ct + 1;
        missed_trial(missed_trial_ct) = i; 
    end;    

%     if (max(diff_chan121314(:,i)) > THRESHOLD)
%         first_rise_loc_diff_chan121314(i) = find(first_rise_loc_diff_chan121314(:,i) > THRESHOLD, 1 )+1;
%     else
%         first_rise_loc_diff_chan121314(i) = -999;
%     end
% % %   figure(2);
% % %   clf(2);
% % %   
% % %   subplot(2,2,1);
% % %   plot(chan12(:,i), 'b');
% % %   hold on;
% % %   plot(chan13(:,i), 'r');
% % %   events = sprintf('channel %d',Evnt_1(i));
% % %   title(events);
% % %   
% % % %   subplot(2,2,2);
% % % %   plot(chan12OR13(:,i), 'k');
% % % %   hold on;
% % % %   plot(first_rise_loc_chan12OR13(i), chan12OR13(first_rise_loc_chan12OR13(i),i), 'go');
% % %     
% % %   subplot(2,2,4);
% % %   plot(diff_chan12OR13(:,i), 'k');
%      switch(first_rise_loc_diff_chan121314)
%         case (chan12(first_rise_loc_diff_chan121314(i),i) > THRESHOLD) 
%             Evnt_1(i) = 12;  
%         case (chan13(first_rise_loc_diff_chan121314(i),i) > THRESHOLD)
%             Evnt_1(i) = 13; 
%         case (chan14(first_rise_loc_diff_chan121314(i),i) > THRESHOLD)
%             Evnt_1(i) = 14;  
%      end   


     if Evnt_1(i) == 0
    
        if     (chan12(first_rise_loc_diff_chan121314(i),i) > THRESHOLD) && ...
             (chan13(first_rise_loc_diff_chan121314(i),i) < BASELINE)&&...
             (chan14(first_rise_loc_diff_chan121314(i),i) < BASELINE)
            Evnt_1(i) = 12;
        elseif (chan13(first_rise_loc_diff_chan121314(i),i) > THRESHOLD) && ...
             (chan12(first_rise_loc_diff_chan121314(i),i) < BASELINE)&&...
             (chan14(first_rise_loc_diff_chan121314(i),i) < BASELINE)
            Evnt_1(i) = 13; 
        elseif (chan14(first_rise_loc_diff_chan121314(i),i) > THRESHOLD) && ...
             (chan12(first_rise_loc_diff_chan121314(i),i) < BASELINE)&&...
             (chan13(first_rise_loc_diff_chan121314(i),i) < BASELINE)
            Evnt_1(i) = 14; 
        end
     end 
end    
 
[Evnt_1]
 
for i = 1:total_trial_n;
    ptname(i, 1:lenTname(i)) = char(dTname(1:lenTname(i),i)');
    
    switch (ptname(i,1:lenTname(i)))
        case '1'
            cond_num(i) = 1;
        case '2'
            cond_num(i) = 2;
        case '3'
            cond_num(i) = 3;        
        case '4'
            cond_num(i) = 4;           
        case '5'
            cond_num(i) = 5;  

        case '6'
            cond_num(i) = 6;    
        case '7'
            cond_num(i) = 7;             
        case '8'
            cond_num(i) = 8;                  
        case '9'
            cond_num(i) = 9;             
        case '10'
            cond_num(i) = 10; 

        case '11'
            cond_num(i) = 11;            
        case '12'
            cond_num(i) = 12;            
        case '13'
            cond_num(i) = 13;            
        case '14'
            cond_num(i) = 14;              
        case '15'
            cond_num(i) = 15; 

        case '16'
            cond_num(i) = 16;
        case '17'
            cond_num(i) = 17;            
        case '18'
            cond_num(i) = 18;              
        case '19'
            cond_num(i) = 19;            
        case '20'
            cond_num(i) = 20;

        case '21'
            cond_num(i) = 21;
        case '22'
            cond_num(i) = 22;
        case '23'
            cond_num(i) = 23;        
        case '24'
            cond_num(i) = 24;           
        case '25'
            cond_num(i) = 25;  

        case '26'
            cond_num(i) = 26;    
        case '27'
            cond_num(i) = 27;             
        case '28'
            cond_num(i) = 28;                  
        case '29'
            cond_num(i) = 29;             
        case '30'
            cond_num(i) = 30; 

        case '31'
            cond_num(i) = 31;            
        case '32'
            cond_num(i) = 32;            
        case '33'
            cond_num(i) = 33;            
        case '34'
            cond_num(i) = 34;              
        case '35'
            cond_num(i) = 35; 

        case '36'
            cond_num(i) = 36;
        case '37'
            cond_num(i) = 37;            
        case '38'
            cond_num(i) = 38;              
        case '39'
            cond_num(i) = 39;            
        case '40'
            cond_num(i) = 40;
                        
        
                        
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

Resp12 = zeros(1,total_cdnum); 
Resp13 = zeros(1,total_cdnum); 
Resp14 = zeros(1,total_cdnum); 
Count = zeros(1,total_cdnum); 

for i = 1:total_trial_n
    if listResp(i,3) == 12
       listResp(i, 4) = 0;
       listResp(i, 5) = 1;
%        Resp(cond_num(i)) = Resp(cond_num(i)) + listResp(i, 4);
        Resp12(cond_num(i)) = Resp12(cond_num(i)) + 1;
       Count(cond_num(i)) = Count(cond_num(i)) + 1;
        
    elseif listResp(i,3) == 13
       listResp(i, 4) = 1;
        listResp(i, 6) = 1;
%        Resp(cond_num(i)) = Resp(cond_num(i)) + listResp(i, 4);
        Resp13(cond_num(i)) = Resp13(cond_num(i)) + 1;
       Count(cond_num(i)) = Count(cond_num(i)) + 1;

   elseif listResp(i,3) == 14
   listResp(i, 4) = 2;
   listResp(i, 7) = 1;
%    Resp(cond_num(i)) = Resp(cond_num(i)) + listResp(i, 4);
    Resp14(cond_num(i)) = Resp14(cond_num(i)) + 1;
   Count(cond_num(i)) = Count(cond_num(i)) + 1;
        
    elseif listResp(i,3) == -3
        listResp(i, 4) = NaN;
    else disp 'error';           
    end
    
end
x2_single_double     = (1:5);
x2_double_single     = (6:10);
x2_Nsingle_Ndouble     = (11:15);
x2_Ndouble_Nsingle     = (16:20);

x4_single_double     = (21:25);
x4_double_single     = (26:30);
x4_Nsingle_Ndouble     = (31:35);
x4_Ndouble_Nsingle     = (36:40);



rest_together = [Resp12',Resp13',Resp14'];
rest_together(:,4) = sum(rest_together,2);

%extraction
x4_single_double_resp = rest_together(x4_single_double,:)+rest_together(x4_Nsingle_Ndouble,:);
x4_double_single_resp1 = rest_together(x4_double_single,:)+rest_together(x4_Ndouble_Nsingle,:);
x4_double_single_resp = x4_double_single_resp1(:,[2 1 3 4]); 
x4_resp = x4_single_double_resp + x4_double_single_resp; 

x2_single_double_resp = rest_together(x2_single_double,:)+rest_together(x2_Nsingle_Ndouble,:);
x2_double_single_resp1 = rest_together(x2_double_single,:)+rest_together(x2_Ndouble_Nsingle,:);
x2_double_single_resp = x2_double_single_resp1(:,[2 1 3 4]); 
x2_resp = x2_single_double_resp + x2_double_single_resp;

%calculation
% x4_resp_hits = x4_resp(:,2)./x4_resp(:,4);
% x4_resp_falseR = x4_resp(:,3)./x4_resp(:,4);
% x4_resp_fa = x4_resp(:,1)./x4_resp(:,4); 

x4_resp_hits = x4_resp(:,2)./(x4_resp(:,1)+x4_resp(:,2));
% x4_resp_falseR = x4_resp(:,3)./x4_resp(:,4);
x4_resp_fa = x4_resp(:,1)./(x4_resp(:,1)+x4_resp(:,2)); 

% x4_resp_hits = (x4_resp(:,2)+(x4_resp(:,3)/2))./x4_resp(:,4);
% % x4_resp_falseR = x4_resp(:,3)./x4_resp(:,4);
% x4_resp_fa = (x4_resp(:,1)+(x4_resp(:,3)/2))./x4_resp(:,4);

% x2_resp_hits = x2_resp(:,2)./x2_resp(:,4);
% x2_resp_falseR = x2_resp(:,3)./x2_resp(:,4);
% x2_resp_fa = x2_resp(:,1)./x2_resp(:,4); 

x2_resp_hits = x2_resp(:,2)./(x2_resp(:,1)+x2_resp(:,2));
% x2_resp_falseR = x2_resp(:,3)./x2_resp(:,4);
x2_resp_fa = x2_resp(:,1)./(x2_resp(:,1)+x2_resp(:,2)); 

% x2_resp_hits = (x2_resp(:,2)+(x2_resp(:,3)/2))./x2_resp(:,4);
% % x2_resp_falseR = x2_resp(:,3)./x2_resp(:,4);
% x2_resp_fa = (x2_resp(:,1)+(x2_resp(:,3)/2))./x2_resp(:,4);

x4_resp_hits(isnan(x4_resp_hits))=0; x4_resp_fa(isnan(x4_resp_fa))=0;
x2_resp_hits(isnan(x2_resp_hits))=0; x2_resp_fa(isnan(x2_resp_fa))=0;



figure; 
subplot(1,2,1)
x = [1 2 3 4 5];
vals = [x4_resp_hits,x4_resp_fa];
b = bar(x,vals);
title('x4 hitRate & falseAlarm')

subplot(1,2,2)
x = [1 2 3 4 5];
vals1 = [x2_resp_hits,x2_resp_fa];
b = bar(x,vals1);
title('x2 hitRate & falseAlarm')

[x4_resp_hits' x4_resp_fa']
[x2_resp_hits' x2_resp_fa']

% dprime calc
for i = 1:length(x4_resp_hits)

    if x4_resp_hits(i) >= 101/102 || x4_resp_fa(i) <= 1/102 || x4_resp_hits(i) <= 1/102 || x4_resp_fa(i) >= 101/102
        dprime_x4(i)= norminv(((100 * x4_resp_hits(i)) + 1)/102) - norminv(((100 * x4_resp_fa(i)) + 1)/102);
    % elseif fa(i) > 1/102 || hit(i) < 101/102
    else
        dprime_x4(i)= norminv(x4_resp_hits(i)) - norminv(x4_resp_fa(i));
    end

    if x2_resp_hits(i) >= 101/102 || x2_resp_fa(i) <= 1/102 || x2_resp_hits(i) <= 1/102 || x2_resp_fa(i) >= 101/102
        dprime_x2(i)= norminv(((100 * x2_resp_hits(i)) + 1)/102) - norminv(((100 * x2_resp_fa(i)) + 1)/102);
    % elseif fa(i) > 1/102 || hit(i) < 101/102
    else
        dprime_x2(i)= norminv(x2_resp_hits(i)) - norminv(x2_resp_fa(i));
    end
end

figure; 
subplot(1,2,1)
x = [1 2 3 4 5];
vals = [dprime_x4];
b = bar(x,vals);
title('x4 dprime')

subplot(1,2,2)
x = [1 2 3 4 5];
vals = [dprime_x2];
b = bar(x,vals);
title('x2 dprime')

out(1,:) = [x4_resp_hits' x4_resp_fa']
out(2,:) = [x2_resp_hits' x2_resp_fa']
out(3,:) = [dprime_x4 dprime_x2]
open out
