clear all; 
close all; 

% % %NP 1.25_2.5
% data(:,1) =[0.0499    0.2231    0.3964    0.5697    0.7430    0.9163    1.0896    1.2629]; % 1p25_2p5_control
% data(:,2) = [0.00	0.05	0.05	0.14	0.64	0.79	0.96	1.00]; %NP 1p25_2p5_control
% 
% %CO 1.25_2.5
data(:,1) =[0.0499    0.2231    0.3964    0.5697    0.7430    0.9163    1.0896    1.2629]; % 1p25_2p5_control
data(:,2) = [0.00	0.08	0.25	0.21	0.71	0.75	0.92	1.00]; 

%IN 20_80
% data(:,1) =[2.6400    2.9900    3.3400    3.6800    4.0300    4.3800    4.7200]; %20_80_control
% data(:,2) = [0.0800    0.2900    0.4200    0.5200    0.8300    0.9200    1.0000]; %IN 20_80_control

% % %NP_20_80
% data(:,1) =[2.6400    2.9900    3.3400    3.6800    4.0300    4.3800    4.7200]; %20_80_control
% data(:,2) = [0.00	0.13	0.42	0.95	0.95	0.96	1.00]; 


data(:,3) = repelem(24,length(data));


options             = struct;   % initialize as an empty struct
options.sigmoidName = 'norm';   % choose a cumulative Gaussian as the sigmoid
% options.expType     = '2AFC';   % choose 2-AFC as the paradigm of the experiment
                                % this sets the guessing rate to .5 and
%                                 % fits the rest of the parameters
% options.sigmoidName = 'weibull';

result = psignifit(data,options);



result.Fit
result.conf_Intervals


plotPsych(result);
%% Weibull fitting

% dataStims = [10, 14, 25, 40, 55, 66, 70]; % these are the same for all participants
% dataProp = [1, 1, 0.777777777777778, 0.722222222222222, 0.722222222222222, 0.500000000000000, 0.555555555555556];
%dataProp = [0.94444, 1, 0.77778, 0.88889, 0.55556, 0.44444, 0.66667];
%f = @(F) (1 - exp(-1*(dataStims./F(1)).^F(2))) - dataProp;
% figure; 
% gamma = 0; lambda = 0.01;
% f = @(F) gamma + (1 - gamma - lambda).*(1 - exp(-1*(X./F(1)).^F(2))) - Y;
% p0 = [1 2];
% p = lsqnonlin(f,p0);
% plot(X,f(p)+Y)
% hold on
% plot(X,Y,'+')


h = findobj(gca,'Type','line');
x=get(h,'Xdata') ;
y=get(h,'Ydata') ;


%%
%Fitting Cumulative Normal Distribution Function to Data
%https://www.mathworks.com/matlabcentral/answers/345595-fitting-cumulative-normal-distribution-function-to-data
fcn = @(b,x) normcdf(x, b(1), b(2));                    % Objective Function
NRCF = @(b) norm(Y - fcn(b,X));                     % Norm Residual Cost Function
B = fminsearch(NRCF, [0; 10]);                          % Estimate Parameters
Xplot = linspace(min(X), max(X));
figure(1)
plot(X, Y, 'pg')
hold on
plot(Xplot, fcn(B,Xplot))
hold off
grid
text(0.2, 0.5, sprintf('\\mu = %.1f\n\\sigma = %.1f', B))
%%