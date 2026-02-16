close all; clear all; 

%%   ANOVA to compare d' of x4 and x2 speed differences across 5 speed
% pairs 3AFC and 2AFC

%import data from excel file

% AFC_3 = xlsread('speed_psycho_excel.xlsx','stat','B4:K7');
% AFC_2 = xlsread('speed_psycho_excel.xlsx','stat','N4:W7');
AFC_3 = xlsread('speed_psycho_excel.xlsx','stat','B37:K40');
% AFC_2 = xlsread('speed_psycho_excel.xlsx','stat','N4:W7');

% 1. 3AFC comparision  ANOVA
for i = 1:5
    [p,tbl,stats] = anova1([AFC_3(:,i), AFC_3(:,i+5)]);
    out_3AFC(i,:) = [tbl{2,3},tbl{3,3},round(tbl{2,5},2),tbl{2,6}];
    clear p tbl stats
end 

% 2. 2AFC comparision  ANOVA
for i = 1:5
    [p,tbl,stats] = anova1([AFC_2(:,i), AFC_2(:,i+5)]);
    out_2AFC(i,:) = [tbl{2,3},tbl{3,3},round(tbl{2,5},2),tbl{2,6}];
    clear p tbl stats
end 

%% 3. one way anova within 3AFC

close all

%%%% 3AFC x4 combo
speed_combo_x4 = {'1.25-5'	'2.5-10' '5-20'	'10-40' '20-80'};
% does d' mean from any of the speed combo differs?
% one way analysis of variance
[p, tbl, stats] = anova1(AFC_3(:,1:5), speed_combo_x4);
% p= 6.2441e-09. Hence, strong evidence that at least one of the d' means differ from the other speed combo.
% F(4,15) = 57.72

% post hoc test: there is difference but where?
[c,m,h,gnames] = multcompare(stats)
tbl1 = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])
%The means of groups _ and _are significantly different if p<=0.05
%'1.25-5'vs'20-80', p= 7.1293e-08***
% '2.5-10'vs'20-80', p= 2.0172e-08***
% '5-20'vs'20-80', p= 2.1808e-08***	
% '10-40'vs'20-80', p= 9.816e-08***
% other combo n.s

[p, tbl, stats] = anova1(AFC_3(:,1:5), speed_combo_x4);

close all;
%%%% 3AFC x2 combo
speed_combo_x2 = {'1.25-2.5' '2.5-5' '5-10'	'10-20' '20-40'};
% does d' mean from any of the speed combo differs?
% one way analysis of variance
[p, tbl, stats] = anova1(AFC_3(:,6:10), speed_combo_x2);
% p= 0.0001. Hence, strong evidence that at least one of the d' means differ from the other speed combo.
% F(4,15) = 12.17

% post hoc test: there is difference but where?
[c,m,h,gnames] = multcompare(stats)
tbl1 = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])
%The means of groups _ and _are significantly different if p<=0.05
%'1.25-2.5'vs'5-10', p=0.013255*
%'2.5-5'vs'20-40', p=0.0033823**
%'5-10'vs'20-40', p=0.00016951***
%'10-20'vs'20-40',p=0.00080898***
%other combo n.s

% Assuming you have loaded your data into MATLAB as a matrix named 'data'

% Extract the data for speed combinations 10&20 and 20&40 for x2
data_10_20_x2 = AFC_3(:,9); % Adjust the column index as per your data structure
data_20_40_x2 = AFC_3(:,10); % Adjust the column index as per your data structure

% Combine the data for ANOVA
grouped_data = [data_10_20_x2; data_20_40_x2];

% Create a grouping variable
groups = [repmat({'10&20'}, size(data_10_20_x2)); repmat({'20&40'}, size(data_20_40_x2))];

% Perform the one-way ANOVA
[p, tbl, stats] = anova1(grouped_data, groups, 'off');  % 'off' suppresses the ANOVA table display

% Display the results
F_statistic = tbl{2,5};
p_value = p;

fprintf('F(1,%d) = %.2f, p = %.5f\n', stats.df, F_statistic, p_value);



%% 4. one way anova within 2AFC

close all

%%%% 2AFC x4 combo
speed_combo_x4 = {'1.25-5'	'2.5-10' '5-20'	'10-40' '20-80'};
% does d' mean from any of the speed combo differ?
% one way analysis of variance
[p, tbl, stats] = anova1(AFC_2(:,1:5), speed_combo_x4);
% p= 1.0250e-08. Hence, strong evidence that at least one of the d' means differ from the other speed combo.
% F(4,15) = 53.76

% post hoc test: there is difference but where?
[c,m,h,gnames] = multcompare(stats)
tbl1 = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])
%The means of groups _ and _are significantly different if p<=0.05
%'1.25-5'vs'20-80', p= 5.3363e-08***
% '2.5-10'vs'20-80', p= 9.8462e-08***
% '5-20'vs'20-80', p= 5.3363e-08***	
% '10-40'vs'20-80', p= 5.3363e-08***
% other combo n.s


close all;
%%%% 2AFC x2 combo
speed_combo_x2 = {'1.25-2.5' '2.5-5' '5-10'	'10-20' '20-40'};
% does d' mean from any of the speed combo differs?
% one way analysis of variance
[p, tbl, stats] = anova1(AFC_2(:,6:10), speed_combo_x2);
% p= 0.0005. Hence, strong evidence that at least one of the d' means differ from the other speed combo.
% F(4,15) = 9.65

% post hoc test: there is difference but where?
[c,m,h,gnames] = multcompare(stats)
tbl1 = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])
%The means of groups _ and _are significantly different if p<=0.05
%'1.25-2.5'vs'10-20', p=0.035894*
%'2.5-5'vs'20-40', p=0.0064474**
%'5-10'vs'20-40', p=0.0018956**
%'10-20'vs'20-40',p=0.0008561***
%other combo n.s


%% 5. Two- way anova to test
% doing a two-way ANOVA, in which the two factors are 
% the speed separation (x2 and x4) and the stimulus speeds (from low to high)

% To test if the main effect of speed separation on d' is significant, and
% if d' depends on the stimulus speed.

%%%% 3AFC x4 combo
speed_combo_x4 = {'1.25-5'	'2.5-10' '5-20'	'10-40' '20-80'};
% [p, tbl, stats] = anova1(AFC_3(:,1:5), speed_combo_x4);

close all; 
data_restru = AFC_3(:,1:5);
data_restru = [data_restru; AFC_3(:,6:10)]

[p,tbl,stats] = anova2(data_restru,4);

% Source          SS      df     MS        F     Prob>F
% -----------------------------------------------------
% Columns        62.417    4   15.6042   46.81   0     
% Rows           31.407    1   31.4068   94.22   0     
% Interaction     7.401    4    1.8503    5.55   0.0018
% Error          10       30    0.3333                 
% Total         111.225   39                           

% p-values for the stimulus speed (from low to high combo) p = 0, F=46.81
% the speed seperation (x2 and x4) is p = 0, F = 94.22
% and the interaction between stimulus speed and speed seperation, p=0.0018

%These values indicate that stimulus speed and speed seperation affect d'
%Also, evidence of an interaction effect of the two. 


% %multicompare
% [~,~,stats]  = anova2(data_restru,4,"off");
% c1 = multcompare(stats);
% c2 = multcompare(stats,"Estimate","row");




%%
%The discrimination performance was significantly better at x4 than at x2 speed separation 


x4 = AFC_3(:,1:5);
x2 = AFC_3(:,6:10);

[p,stats] = anova1([x4(:,1) x2(:,1)]);
%F(1,6)=51.84, p = 0.0004
[p,stats] = anova1([x4(:,2) x2(:,2)]);
%F(1,6)=45.33, p = 0.0005
[p,stats] = anova1([x4(:,3) x2(:,3)]);
%F(1,6)=10.99, p = 0.0161
[p,stats] = anova1([x4(:,4) x2(:,4)]);
%F(1,6) = 12.6, p = 0.0121
[p,stats] = anova1([x4(:,5) x2(:,5)]);
%F(1,6) = 1.5, p = 0.266

n = 5;
alpha = 0.05;
alpha_bonf = alpha / n;


[h,p,ci,stats] = ttest(x4(:,1) ,x2(:,1))
%p = 0.0015 tval= 11.1616
[h,p,ci,stats] = ttest(x4(:,2), x2(:,2))
 %p = 0.0093 tval= 5.9915
[h,p,ci,stats] = ttest(x4(:,3), x2(:,3))
 %p = 0.0405 tval= 3.4654
[h,p,ci,stats] = ttest(x4(:,4), x2(:,4))
 %p = 0.0447 tval= 3.3291
[h,p,ci,stats] = ttest(x4(:,5), x2(:,5))
% p=0.2818 tval= 1.3090


%The discrimination performance was poor at fast speeds. d' dropped significantly as the stimulus speeds increased from 10 and 40°/s to 20 and 80°/s for x4 separation (one-way ANOVA, F(1,6) = 67.6, p = 1.75x10-4) (Fig. 1E), and as the stimulus speeds increased from 10 and 20°/s to 20 and 40°/s  for x2 separation (one-way ANOVA, F(1,6) = 47.1, p = 4.71x10-4 ) (Fig. 1F).  

[h,p,ci,stats] = ttest(x4(:,4) ,x4(:,5))
%9.7062e-04,tstat= 13.0549
[h,p,ci,stats] = ttest(x2(:,4) ,x2(:,5))
%p=0.0072, tstat = 6.5583

% Perform Bonferroni correction
c = multcompare(stats, 'CType', 'bonferroni');













