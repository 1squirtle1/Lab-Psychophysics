clear; clc; 
% close all;

%% 1. SETUP
%--------------------------------------------------------------------------
% Define file, sheet, and experimental parameters
fileName    = 'CollectedData.xlsx';
sheetName   = 'threshold_vals'; 
participants = [1, 2, 3]; % Participant IDs
eccAngles   = [3, 4, 5, 6, 7]; % Eccentricities in degrees

% Define the row ranges for each participant in the Excel sheet
% p_rows = {1:5, 7:11, 13:17,19:23}; 
p_rows = {1:5, 7:11, 13:17}; 

%% 2. LOAD AND PARSE RAW SESSION DATA
%--------------------------------------------------------------------------
% Load the entire sheet of raw session data
rawData = readmatrix(fileName, 'Sheet', sheetName);

% Initialize a table to store data in long format (one row per session)
T_long = table();

% Loop through each participant and eccentricity to build the long table
for p_idx = 1:numel(participants)
    p_id = participants(p_idx);
    rows = p_rows{p_idx};
    
    for e_idx = 1:numel(eccAngles)
        ecc = eccAngles(e_idx);
        
        % Get threshold values for all sessions for this condition
        session_vals = rawData(rows(e_idx), 1:11);
        
        % Remove missing sessions (which read as NaN)
        session_vals = session_vals(~isnan(session_vals));
        
        % Create a temporary table for this condition's data
        n_sessions = numel(session_vals);
        temp_tbl = table( ...
            repmat(p_id, n_sessions, 1), ...
            repmat(ecc, n_sessions, 1), ...
            session_vals', ...
            'VariableNames', {'ParticipantID', 'Eccentricity', 'Threshold'});
            
        % Append this data to our main long-format table
        T_long = [T_long; temp_tbl];
    end
end

fprintf('Successfully loaded and parsed data for %d total sessions.\n', height(T_long));

%% 3. FIT LINEAR MIXED-EFFECTS MODEL
%--------------------------------------------------------------------------
% This model finds the overall trend (fixed effect of eccentricity) while
% allowing each participant to have their own baseline (random intercept)
% AND their own unique trend (random slope). This is a powerful way to
% model individual differences.
lme_formula = 'Threshold ~ 1 + Eccentricity + (1 + Eccentricity | ParticipantID)';
lme = fitlme(T_long, lme_formula);

% Extract key results from the model
fixed_effects = lme.Coefficients;
slope      = fixed_effects.Estimate(2);
intercept  = fixed_effects.Estimate(1);
slope_CI   = coefCI(lme, 'alpha', 0.05); % 95% CI
slope_CI   = slope_CI(2,:);
adjR2      = lme.Rsquared.Adjusted;

fprintf('\n*** Linear Mixed-Effects Model Results ***\n');
disp(lme);
fprintf('Fixed Effect of Eccentricity (Slope): β = %.4f, 95%% CI [%.4f, %.4f]\n', ...
    slope, slope_CI(1), slope_CI(2));
fprintf('Model Fit: Adjusted R² = %.3f\n\n', adjR2);

% Add this after Section 3
figure;
subplot(1,2,1);
plotResiduals(lme, 'fitted'); % Check for homoscedasticity (should look like a random cloud)
title('Residuals vs Fitted');
subplot(1,2,2);
plotResiduals(lme, 'probability'); % Check for normality (should follow the straight line)
title('Normal Probability Plot of Residuals');

%% 4. CALCULATE SUMMARY STATISTICS FOR PLOTTING
%--------------------------------------------------------------------------
% For plotting, we calculate the mean and SEM for each condition across sessions.
[G, plot_tbl] = findgroups(T_long(:, {'ParticipantID', 'Eccentricity'}));
plot_tbl.mean_thresh = splitapply(@mean, T_long.Threshold, G);
plot_tbl.sem_thresh  = splitapply(@(x) std(x)/sqrt(numel(x)), T_long.Threshold, G);

%% 5. GENERATE LME PREDICTIONS FOR PLOT
%--------------------------------------------------------------------------
% Create a smooth range of eccentricities for the fit line
xFit = linspace(min(eccAngles), max(eccAngles), 100)';
newTbl = table(xFit, repmat(participants(1), numel(xFit), 1), ...
    'VariableNames', {'Eccentricity', 'ParticipantID'});

% Get fixed-effects predictions only
yFit = predict(lme, newTbl, 'Conditional', false);
% yFit = feval(lme, xFit); % Get fixed-effects predictions


%% 6. PLOT DATA AND MODEL
%--------------------------------------------------------------------------

hold off;


figure('Color', 'w');
hold on;
box on;
ax = gca;
ax.LineWidth = 1.5;
ax.FontSize = 12;

colors = [0.93, 0.69, 0.13; ... % Gold for P1
          0,    0.45, 0.74; ... % Blue for P2
%           0.47, 0.67, 0.19;...
          0.17, 0.27, 0.19];   % Green for P3

% Plot individual participant means with error bars
for p_idx = 1:numel(participants)
    p_id = participants(p_idx);
    p_data = plot_tbl(plot_tbl.ParticipantID == p_id, :);
    
    errorbar(p_data.Eccentricity, p_data.mean_thresh, p_data.sem_thresh, ...
        'o-', 'Color', colors(p_idx,:), 'LineWidth', 1.5, ...
        'MarkerSize', 8, 'MarkerFaceColor', colors(p_idx,:), 'CapSize', 0);
end

% Plot the LME model fit line
plot(xFit, yFit, 'k-', 'LineWidth', 2.5);

% --- ADD THIS BLOCK BELOW ---
% Plot participant-specific fits (random effects included)
for p = participants
    newTbl = table(xFit, repmat(p, numel(xFit), 1), ...
        'VariableNames', {'Eccentricity', 'ParticipantID'});
    yFit_p = predict(lme, newTbl, 'Conditional', true);
    plot(xFit, yFit_p, '--', 'Color', colors(p == participants,:), 'LineWidth', 1);
end
% --- END ADD ---

% Add labels, title, and legend
xlabel('Eccentricity (deg)');
ylabel('Critical Spatial Separation (deg)');
title('Perceptual Thresholds vs. Eccentricity');
legend_labels = arrayfun(@(p) sprintf('Participant %d', p), participants, 'uni', 0);
legend([legend_labels, {'LME Fit'}], 'Location', 'NorthWest', 'Box', 'off');

hold off;

%--------------------------------------------------------------------------
%% 7. ROBUST EXTRACTION OF FIXED & RANDOM EFFECTS
%--------------------------------------------------------------------------
re_raw = randomEffects(lme);
fixed_raw = lme.Coefficients;

fprintf('--- Diagnostics ---\n');
fprintf('randomEffects class: %s\n', class(re_raw));
fprintf('lme.Coefficients class: %s\n', class(fixed_raw));
fprintf('-------------------\n');

% Handle fixed effects (works for titleddataset or table)
if istable(fixed_raw)
    fixed_tbl = fixed_raw;
elseif isa(fixed_raw, 'classreg.regr.lmeutils.titleddataset')
    fixed_tbl = dataset2table(fixed_raw);
elseif isstruct(fixed_raw)
    fixed_tbl = struct2table(fixed_raw);
else
    error(['Unexpected type for lme.Coefficients: ' class(fixed_raw)]);
end

fixed_est   = fixed_tbl.Estimate;
fixed_names = fixed_tbl.Name;

fixed_int_idx = strcmp(fixed_names,'(Intercept)');
fixed_slope_idx = strcmp(fixed_names,'Eccentricity');

fixed_int = fixed_est(fixed_int_idx);
fixed_slope = fixed_est(fixed_slope_idx);

% Handle random effects (double for R2020b)
if istable(re_raw)
    re_tbl = re_raw;
elseif isstruct(re_raw)
    re_tbl = struct2table(re_raw);
elseif isnumeric(re_raw)
    nParams = 2; % (Intercept), Eccentricity
    nParticipants = numel(participants);
    if numel(re_raw) ~= nParams * nParticipants
        error('Unexpected size of randomEffects output.');
    end
    re_mat = reshape(re_raw, [nParams, nParticipants])';
    re_tbl = table();
    re_tbl.Level = participants(:);
    re_tbl.Intercept = re_mat(:,1);
    re_tbl.Slope = re_mat(:,2);
else
    error(['randomEffects returned unknown type: ' class(re_raw)]);
end

%% 8. BUILD PARTICIPANT SUMMARY TABLE
%--------------------------------------------------------------------------
participant_summary = table();

for i = 1:numel(participants)
    pid = participants(i);
    if all(ismember({'Intercept','Slope'}, re_tbl.Properties.VariableNames))
        rows_int = re_tbl.Intercept(re_tbl.Level == pid);
        rows_slope = re_tbl.Slope(re_tbl.Level == pid);
    else
        rows_int = 0; rows_slope = 0;
    end
    
    intercept_p = fixed_int + rows_int;
    slope_p = fixed_slope + rows_slope;
    
    subj_data = T_long(T_long.ParticipantID == pid,:);
    mean_thresh = mean(subj_data.Threshold);
    sem_thresh  = std(subj_data.Threshold)/sqrt(height(subj_data));
    
    participant_summary = [participant_summary;
        table(pid, intercept_p, slope_p, mean_thresh, sem_thresh, ...
        'VariableNames', {'ParticipantID','PredIntercept','PredSlope','MeanThreshold','SEMThreshold'})];
end

disp('Participant summary:');
disp(participant_summary);
writetable(participant_summary, 'Participant_LME_Summary.xlsx');




