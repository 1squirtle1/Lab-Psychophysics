clear all; 
close all;

% Define the path to the Excel file
filePath = 'CollectedData.xlsx';
% Get sheet names
[~, sheetNames] = xlsfinfo(filePath);
pathlength= (0.1:0.1:1);
% % pathlength = pathlength.*sin(90/2);
pathlength = pathlength.*(sqrt(2)/2);

% Ask the user to input the participant's name
participant_name = input('Enter the participant name (e.g., rb, ap, bg, sbw): ', 's');

% Define dg_start based on participant_name
switch participant_name
    case 'rb'
        dg_start = {(1:12), (21:33), (42:55), (64:76), (85:98)}; % RB
        participant_no = 1;
    case 'ap'
        dg_start = {(1:13), (23:35), (45:59), (68:82), (92:102)}; % AB
        participant_no =2; 
    case 'bg'
        dg_start = {(1:9), (18:25), (34:43), (52:60), (69:77)}; % BG
        participant_no = 3; 
    case 'sbw'
        dg_start = {(1:10), (17:26), (33:42), (49:58), (65:74)}; % SBW
        participant_no = 4; 
    otherwise
        error('Invalid participant name. Please enter "rb", "ab", or "bg".');
end



% % Define Weibull function
% % weibull_func = @(p, x) 1 - exp(-((x./p(1)).^p(2)));
% weibull_func = @(p, x) p(1) * (1 - exp(-(x / p(2)).^p(3))); 
% logistic_func = @(p, x) p(1) ./ (1 + exp(-p(2) * (x - p(3))));
% logistic_func = @(p, x) p1 ./ (1 + exp(-p(1) * (x - p(2))));

% Initialize cell arrays for storing different sizes
reval_Dir = cell(1, 5);
reval_Trials = cell(1, 5);
reval_Data_Motion = cell(1, 5);
reval_Data_NoMotion = cell(1, 5);



% Loop through participants
for i = participant_no:participant_no

    % Read the sheet data
    sheetData = readtable(filePath, 'Sheet', sheetNames{i});
    colors = lines(5);
%     figure; 
    for dg = 1:5
        % Extract the direction (character) column from specified rows
        reval_Dir{dg} = table2cell(sheetData(dg_start{dg}, 5)); % Store in cell array
        
        % Extract the trials (numeric) column from specified rows
        reval_Trials{dg} = table2array(sheetData(dg_start{dg}, 6)); % Store in cell array
        
        % Extract the motion data (numeric) columns from specified rows
        reval_Data_Motion_Seg{dg} = table2array(sheetData(dg_start{dg}, 8:17)); % Store in cell array
        
        % Extract the no-motion data (numeric) columns from specified rows
        reval_Data_NoMotion{dg} = table2array(sheetData(dg_start{dg}, 19:28)); % Store in cell array
        
        merge_index = 1;
        sessions = 1;
            while sessions <= length(reval_Trials{dg})
                if reval_Trials{dg}(sessions) == 5; 
                    % Merge current session with the next one
                    reval_Trials_merge{dg}(merge_index) = reval_Trials{dg}(sessions)+reval_Trials{dg}(sessions+1);
                    reval_Data_Motion_Seg_merge{dg}(merge_index,:) = reval_Data_Motion_Seg{dg}(sessions,:)+reval_Data_Motion_Seg{dg}(sessions+1,:);
                    reval_Data_NoMotion_merge{dg}(merge_index,:) = reval_Data_NoMotion{dg}(sessions,:)+reval_Data_NoMotion{dg}(sessions+1,:);
                    reval_Data_Motion_Int_merge{dg}(merge_index,:) = 10-(reval_Data_Motion_Seg_merge{dg}(merge_index,:)+reval_Data_NoMotion_merge{dg}(merge_index,:));

                    % Skip to the session after the next one, as we've merged two
                    sessions = sessions+2;
                else
                    % If current value is not 5, just add it as is to the merged array
                    reval_Trials_merge{dg}(merge_index) = reval_Trials{dg}(sessions);
                    reval_Data_Motion_Seg_merge{dg}(merge_index,:) = reval_Data_Motion_Seg{dg}(sessions,:);
                    reval_Data_NoMotion_merge{dg}(merge_index,:) = reval_Data_NoMotion{dg}(sessions,:);
                    reval_Data_Motion_Int_merge{dg}(merge_index,:) = 10-(reval_Data_Motion_Seg_merge{dg}(merge_index,:)+reval_Data_NoMotion_merge{dg}(merge_index,:));
                    % Move to the next session
                    sessions = sessions + 1;
                end 
                merge_index = merge_index + 1;
            end 

        % Compute seg_percent based on different lengths
%         seg_percent = zeros(size(reval_Data_Motion{dg}, 1), size(reval_Data_Motion{dg}, 2)); % Initialize seg_percent
     
        for trials = 1:length(reval_Trials_merge{dg})
            % seg_percent(trials, :) = reval_Data_Motion_Seg_merge{dg}(trials, :) ./ (reval_Data_Motion_Seg_merge{dg}(trials, :)+reval_Data_Motion_Int_merge{dg}(trials, :));
            seg_percent_all(trials, :) = reval_Data_Motion_Seg_merge{dg}(trials, :) ./ (reval_Data_Motion_Seg_merge{dg}(trials, :)+reval_Data_Motion_Int_merge{dg}(trials, :)+reval_Data_NoMotion_merge{dg}(trials, :));
           
            
            % int_percent(trials, :) = reval_Data_Motion_Int_merge{dg}(trials, :) ./ (reval_Data_Motion_Seg_merge{dg}(trials, :)+reval_Data_Motion_Int_merge{dg}(trials, :));
            int_percent_all(trials, :) = reval_Data_Motion_Int_merge{dg}(trials, :) ./ (reval_Data_Motion_Seg_merge{dg}(trials, :)+reval_Data_Motion_Int_merge{dg}(trials, :)+reval_Data_NoMotion_merge{dg}(trials, :));

            noMotion_percent_all(trials, :) = reval_Data_NoMotion_merge{dg}(trials, :) ./ (reval_Data_Motion_Seg_merge{dg}(trials, :)+reval_Data_Motion_Int_merge{dg}(trials, :)+reval_Data_NoMotion_merge{dg}(trials, :));
%             seg_percent(trials, :) = reval_Data_Motion{dg}(trials, :) ./ (reval_Trials{dg}(trials)-reval_Data_NoMotion{dg}(trials, :));

%             seg_percent(trials, :) = (reval_Trials{dg}(trials) - reval_Data_Motion{dg}(trials, :) - reval_Data_NoMotion{dg}(trials, :))./...
%                 reval_Trials{dg}(trials);
%             seg_percent(trials, :) = reval_Data_NoMotion{dg}(trials, :) ./ reval_Trials{dg}(trials);

            % coherentMotion_percent_all(trials, :) = (reval_Data_Motion_Int_merge{dg}(trials, :) + reval_Data_Motion_Seg_merge{dg}(trials, :))./...
                % (reval_Data_Motion_Int_merge{dg}(trials, :) + reval_Data_Motion_Seg_merge{dg}(trials, :) + reval_Data_NoMotion_merge{dg}(trials, :));

        end 
        
        % seg_percent_all_mean(:,dg) = (nanmean(seg_percent_all*100))';
        % seg_percent_all_sem(:,dg) = (nanstd(seg_percent_all*100) ./ sqrt(sum(~isnan(seg_percent_all*100))))';

        % seg_percent_all_mean(:,dg) = (nanmean(seg_percent_all*100))';
        % seg_percent_all_sem(:,dg) = (nanstd(seg_percent_all*100) ./ sqrt(sum(~isnan(seg_percent_all*100))))';

        seg_percent_all_mean(:,dg) = (nanmean(seg_percent_all))';
        seg_percent_all_sem(:,dg) = (nanstd(seg_percent_all) ./ sqrt(sum(~isnan(seg_percent_all))))';

        int_percent_all_mean(:,dg) = (nanmean(int_percent_all))';
        int_percent_all_sem(:,dg) = (nanstd(int_percent_all) ./ sqrt(sum(~isnan(int_percent_all))))';

        noMotion_percent_all_mean(:,dg) = (nanmean(noMotion_percent_all))';
        noMotion_percent_all_sem(:,dg) = (nanstd(noMotion_percent_all) ./ sqrt(sum(~isnan(noMotion_percent_all))))';

        % coherentMotion_percent_all_mean(:,dg) = (nanmean(coherentMotion_percent_all))';
        % coherentMotion_percent_all_sem(:,dg) = (nanstd(coherentMotion_percent_all) ./ sqrt(sum(~isnan(coherentMotion_percent_all))))';


        current_color = colors(dg, :);
        figure(6)
        subplot(1,4,1)
        errorbar(pathlength,seg_percent_all_mean(:,dg),seg_percent_all_sem(:,dg) ,'o-','Color', current_color); hold on;
        title("segmentation_percentage")
        
        subplot(1,4,2)
        errorbar(pathlength, int_percent_all_mean(:,dg), int_percent_all_sem(:,dg),'o-','Color', current_color); hold on;
        title("integration_percentage")

        subplot(1,4,3)
        errorbar(pathlength, noMotion_percent_all_mean(:,dg), noMotion_percent_all_sem(:,dg),'o-','Color', current_color); hold on;
        title("integration_percentage")

        % subplot(1,4,4)
        % errorbar(pathlength, nanmean(coherentMotion_percent_all), nanstd(coherentMotion_percent_all) ./ sqrt(sum(~isnan(coherentMotion_percent_all))),'o-','Color', current_color); hold on;
        % title("coherent motion")

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% Fit logistic trial by trial and at each  eccentricity
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(7)
        for t = 1:size(seg_percent_all,1)  
            x_values = pathlength;
            y_values = seg_percent_all(t,:); % Normalize the y-values
           
             % 1) SET UP LOGISTIC
            p1 = max(y_values);  
            logistic_func = @(p,x) p1 ./ (1 + exp(-p(1)*(x - p(2))));
            p0 = [1, mean(x_values)];          
            lb = [0, min(x_values)];
            ub = [50, max(x_values)];

            % 2) FIT
            p_fit = lsqcurvefit(logistic_func, p0, x_values, y_values, lb, ub);
    
            % 3) STORE FITTED CURVE
            x_fit = linspace(min(x_values), max(x_values), 100);
            y_fit(t,:,dg) = logistic_func(p_fit, x_fit);

            % 4) EXTRACT PARAMETERS
            p2(dg,t) = p_fit(1);    % slope
            p3(dg,t) = p_fit(2);    % inflection = 50% point
            y_50      = 0.5 * p1;
            threshold_50(dg,t) = p3(dg,t) - (1/p2(dg,t))*log((p1/y_50)-1);

            % 5) PLOT
            plot(x_fit, y_fit(t,:,dg), '--', 'Color', current_color);
            
        end 
        hold on;

    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Bootstrap values trial by traial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set bootstrap parameters
rng(12345,'twister');
nBoot = 1000;       % Number of bootstrap iterations
ci = 95;            % Confidence interval percentage
lowerPerc = (100 - ci) / 2;
upperPerc = 100 - lowerPerc;

nStimuli = size(threshold_50, 1);

% Initialize result vectors
bootMean    = zeros(nStimuli, 1);
bootCI_lower = zeros(nStimuli, 1);
bootCI_upper = zeros(nStimuli, 1);

for i = 1:nStimuli
    % Remove 0's from the data of stimulus i
    data = threshold_50(i, :);
    data = data(data ~= 0);
    reg_mean(i) = mean(data);
    reg_std(i) = std(data);
    reg_sem(i) = std(data)./sqrt(length(data));
    if isempty(data)
        bootMean(i)    = NaN;
        bootCI_lower(i)= NaN;
        bootCI_upper(i)= NaN;
    else
        nData = numel(data);
        bootMeans = zeros(nBoot, 1);
        
        % Perform bootstrap resampling
        for j = 1:nBoot
            sampleIndices = randi(nData, [1, nData]);  % sample with replacement
            sample = data(sampleIndices);
            bootMeans(j) = mean(sample);
        end
        
        % Compute bootstrap mean and CI percentiles
        bootSamples(:,i) = bootMeans;
        bootMean(i)     = mean(bootMeans);
        bootCI_lower(i) = prctile(bootMeans, lowerPerc);
        bootCI_upper(i) = prctile(bootMeans, upperPerc);
        % error(i) = (bootCI_upper - bootCI_lower) / 2;

    end
    bootMeans_compile(:,i) = bootMeans;
end

[reg_mean', reg_std', reg_sem']
% ——————————————————————————————————————————————
% PAIRWISE BOOTSTRAP SIGNIFICANCE TESTS
% ——————————————————————————————————————————————
pVals = nan(nStimuli);  % symmetric matrix of p-values

for i = 1:nStimuli
    for k = i+1:nStimuli
        diffs = bootSamples(:,i) - bootSamples(:,k);
        % two‐tailed p-value
        p = 2 * min( mean(diffs>0), mean(diffs<0) );
        pVals(i,k) = p;
        pVals(k,i) = p;
    end
end

% ——————————————————————————————————————————————
% DISPLAY RESULTS
% ——————————————————————————————————————————————
fprintf('Ecc\tMean\t95%% CI Lower\t95%% CI Upper\n');
for i = 1:nStimuli
    fprintf('%2d:\t%.3f\t%.3f\t\t%.3f\n', ...
        i, bootMean(i), bootCI_lower(i), bootCI_upper(i));
end

disp('Pairwise bootstrap p-values (rows vs columns):');
disp(pVals);


error_vals = (bootCI_upper - bootCI_lower) / 2;

% Plot error bars
figure;
x = 1:nStimuli;
errorbar(x, bootMean, bootMean - bootCI_lower, bootCI_upper - bootMean, 'o', 'LineWidth', 1.5);
xlabel('Stimulus');
ylabel('Mean Threshold');
title('Bootstrap Mean and 95% CI for Each Stimulus');
grid on;


% Print results for each stimulus
for i = 1:nStimuli
    fprintf('Stimulus %d: Mean = %.3f, 95%% CI = [%.3f, %.3f], Error = %.3f\n', ...
        i, bootMean(i), bootCI_lower(i), bootCI_upper(i), error_vals(i));
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Fit logistic to the mean response at each eccentricity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%——— Prepare two figures ———
figure(8); clf; hold on
title('Logistic Fits'); 
xlabel('Pathlength (deg)'); 
ylabel('Segmentation (%)');

figure(9); clf; hold on
title('Cumulative-Gaussian Fits'); 
xlabel('Pathlength (deg)'); 
ylabel('Segmentation (%)');

nEcc = size(seg_percent_all_mean,2);
R2_log   = nan(1,nEcc);
R2_gauss = nan(1,nEcc);

for t = 1:nEcc   % looping over different eccentricities  
    x_values = pathlength';
    y_values = seg_percent_all_mean(:,t);
    y_sem    = seg_percent_all_sem(:,t);

    %——— 1) LOGISTIC FIT ———
    p1 = max(y_values);  
    logistic_func = @(p,x) p1 ./ (1 + exp(-p(1) * (x - p(2))));
    p0 = [1, mean(x_values)];  
    lower_bounds = [0, min(x_values)];
    upper_bounds = [50, max(x_values)];
    p_fit = lsqcurvefit(logistic_func, p0, x_values, y_values, lower_bounds, upper_bounds);

    % generate logistic curve
    x_fit_mean = linspace(min(x_values), max(x_values), 100);
    y_fit_log  = logistic_func(p_fit, x_fit_mean);

    % extract logistic params & threshold
    p2(t) = p_fit(1);  
    p3(t) = p_fit(2);  
    y_50  = 0.5 * p1;
    threshold_logistic(t) = p3(t) - (1/p2(t)) * log((p1/y_50) - 1);
    %——— plot into Fig. 8 ———
    figure(8);
    errorbar(x_values, y_values, y_sem, 'o','Color',current_color);
    plot(   x_fit_mean, y_fit_log,   '--','Color',current_color);

    %——— 2) CUMULATIVE-GAUSSIAN FIT ———
    gauss_func = @(q,x) p1 * 0.5 .* (1 + erf((x - q(2)) ./ (q(1)*sqrt(2))));
    q0     = [std(x_values), mean(x_values)];  
    lb_g   = [eps,          min(x_values)];
    ub_g   = [max(x_values), max(x_values)];
    q_fit  = lsqcurvefit(gauss_func, q0, x_values, y_values, lb_g, ub_g);

    % generate Gaussian curve
    y_fit_gauss = gauss_func(q_fit, x_fit_mean);

    % store 50% point for Gaussian
    threshold_gauss(t) = q_fit(2);

    %——— plot into Fig. 9 ———
    figure(9);
    errorbar(x_values, y_values, y_sem, 'o','Color',current_color);
    plot(   x_fit_mean, y_fit_gauss,':','Color',current_color);

    %——— 3) GOODNESS-OF-FIT ———
    yhat_log   = logistic_func(p_fit,   x_values);
    yhat_gauss = gauss_func(  q_fit,    x_values);
    SSE_log    = sum((y_values - yhat_log).^2);
    SSE_gauss  = sum((y_values - yhat_gauss).^2);
    SST        = sum((y_values - mean(y_values)).^2);
    R2_log(t)   = 1 - SSE_log./SST;
    R2_gauss(t) = 1 - SSE_gauss./SST;

    % %——— 4) PLOT DATA + FITS ———
    % errorbar(x_values, y_values, y_sem, 'o', 'Color', current_color);
    % plot(   x_fit_mean, y_fit_log,   '--', 'Color', current_color);
    % plot(   x_fit_mean, y_fit_gauss, ':',  'Color', current_color);
end

% after the loop you can inspect
disp(table((1:nEcc)', R2_log', R2_gauss', ...
    'VariableNames',{'EccIdx','R2_Logistic','R2_Gaussian'}));