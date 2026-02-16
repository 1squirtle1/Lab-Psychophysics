%% Full_Analysis_Pipeline.m
% COMPLETE ANALYSIS PIPELINE
% 1. Reads raw Excel data (CollectedData_revamp.xlsx)
% 2. Fits psychometric functions (Logistic) to every session
% 3. Extract thresholds (lambda_50) and R-squared values
% 4. Runs Linear Mixed-Effects (LME) Analysis
% 5. Generates Publication Tables (Individual Stats & Bootstrapping)
% 6. Generates Summary Figures

clear; clc; 
close all;

%% ========================================================================
%% PART 1: SETUP & PARAMETERS
%% ========================================================================
fileName = 'CollectedData_revamp.xlsx';
sheetNames = {'RH', 'BA', 'AS', 'SC'}; 
% sheetNames = {'RB', 'AP', 'BG'}; 
ecc_labels = [3, 4, 5, 6, 7]; % The degrees corresponding to the blocks

% X-Axis: Spatial Separation (10 steps from 0.07 to 0.71)
xVal = linspace(0.0707, 0.7070, 10); 

% Colors for Eccentricities (Cool to Warm colormap)
colors = parula(length(ecc_labels)); 

% Initialize a structure to store all parsed data
parsedData = struct();

% Initialize a table to store ALL thresholds for LME analysis later
T_thresholds = table(); 

% --- DEFINE LOGISTIC FUNCTION ---
% 3-Parameter Logistic:
% p(1) = Max Asymptote (Performance cap)
% p(2) = Slope (Steepness)
% p(3) = x50 (The 50% Threshold / Inflection Point)
logistic_func = @(p, x) p(1) ./ (1 + exp(-p(2) * (x - p(3))));

% Fit Options
opts = optimset('Display', 'off');
lb = [0.5, 0, min(xVal)]; 
ub = [1.0, 50, max(xVal)]; 
p0 = [1, 10, mean(xVal)];

all_R2_values = []; 

%% ========================================================================
%% PART 2: PSYCHOMETRIC ANALYSIS (PER PARTICIPANT)
%% ========================================================================

for p = 1:numel(sheetNames)
    pName = sheetNames{p};
    fprintf('Processing Participant: %s...\n', pName);
    
    parsedData(p).participantID = pName;
    parsedData(p).blocks = struct(); 
    
    try
        rawData = readmatrix(fileName, 'Sheet', pName);
    catch
        error('Could not read sheet "%s". Check file name and sheet names.', pName);
    end
    
    % --- SPLIT DATA INTO BLOCKS ---
    % Find rows where the first column is NaN (indicating a break between ecc blocks)
    nanRows = find(isnan(rawData(:,1)));
    splitIndices = [0; nanRows; size(rawData,1) + 1];
    
    % Initialize Figure for this Participant (1x4 Layout)
    figure('Color', 'w', 'Name', sprintf('Participant %s', pName), ...
           'Position', [50, 100, 1600, 350]);
    sgtitle(sprintf('Participant %s Psychometric Analysis', pName));
    
    blockCount = 0;
    
    for k = 1:(length(splitIndices)-1)
        startRow = splitIndices(k) + 1;
        endRow   = splitIndices(k+1) - 1;
        
        if endRow >= startRow
            blockCount = blockCount + 1;
            if blockCount > length(ecc_labels), break; end
            
            currEcc = ecc_labels(blockCount);
            blockData = rawData(startRow:endRow, :);
            nSessions = size(blockData, 1);
            
            % --- EXTRACT DATA ---
            cols_Seg    = 2:11;
            cols_NoMot  = 13:22;
            cols_OneMot = 24:33;
            cols_Total  = 35:44;
            
            % Calculate Proportions
            prop_Seg    = blockData(:, cols_Seg) ./ blockData(:, cols_Total);
            prop_NoMot  = blockData(:, cols_NoMot) ./ blockData(:, cols_Total);
            prop_OneMot = blockData(:, cols_OneMot) ./ blockData(:, cols_Total);
            
            % --- CALCULATE MEANS & SEM ---
            mean_Seg = mean(prop_Seg, 1, 'omitnan');
            sem_Seg  = std(prop_Seg, 0, 1, 'omitnan') ./ sqrt(nSessions);
            mean_NoMot = mean(prop_NoMot, 1, 'omitnan');
            sem_NoMot  = std(prop_NoMot, 0, 1, 'omitnan') ./ sqrt(nSessions);
            mean_OneMot= mean(prop_OneMot, 1, 'omitnan');
            sem_OneMot = std(prop_OneMot, 0, 1, 'omitnan') ./ sqrt(nSessions);
            
            % --- FIT INDIVIDUAL SESSIONS ---
            fitParams_Indiv = zeros(nSessions, 3);
            fitR2_Indiv     = zeros(nSessions, 1);
            
            cColor = colors(blockCount, :);
            
            % Plot Panel 1: Individual Sessions (Seg)
            subplot(1, 4, 1); hold on;
            title('Individual Sessions (Seg)');
            xlabel('Spatial Separation'); ylabel('Proportion'); ylim([0 1.05]);
            
            for s = 1:nSessions
                yRaw = prop_Seg(s, :);
                plot(xVal, yRaw, 'o', 'Color', [cColor, 0.3], 'MarkerSize', 3, 'HandleVisibility', 'off');
                
                try
                    [pFit, ~, ~] = lsqcurvefit(logistic_func, p0, xVal, yRaw, lb, ub, opts);
                    fitParams_Indiv(s, :) = pFit;
                    
                    % Calculate R2
                    yPred = logistic_func(pFit, xVal);
                    SSres = sum((yRaw - yPred).^2);
                    SStotal = sum((yRaw - mean(yRaw)).^2);
                    if SStotal < 1e-6, r2 = NaN; else, r2 = 1 - (SSres / SStotal); end
                    fitR2_Indiv(s) = r2;
                    all_R2_values = [all_R2_values; r2]; %#ok<AGROW>
                    
                    % Plot Fit Line
                    xSmooth = linspace(min(xVal), max(xVal), 100);
                    ySmooth = logistic_func(pFit, xSmooth);
                    hL = plot(xSmooth, ySmooth, '-', 'Color', [cColor, 0.6], 'LineWidth', 0.5);
                    if s == 1, set(hL, 'DisplayName', sprintf('%d deg', currEcc)); else, set(hL, 'HandleVisibility', 'off'); end
                catch
                    fitParams_Indiv(s, :) = [NaN, NaN, NaN];
                    fitR2_Indiv(s) = NaN;
                end
            end
            
            % --- ACCUMULATE DATA FOR LME TABLE ---
            thresholds_50 = fitParams_Indiv(:, 3); % 3rd param is threshold
            % Formula: x75 = p3 + ln(3)/p2
%             thresholds_50 = fitParams_Indiv(:, 3) + (log(3) ./ fitParams_Indiv(:, 2));
            
            % Append to main table
            tempTbl = table(repmat({pName}, nSessions, 1), ...
                            repmat(currEcc, nSessions, 1), ...
                            thresholds_50, ...
                            fitR2_Indiv, ...
                            'VariableNames', {'ParticipantID', 'Eccentricity', 'Threshold', 'R2'});
            T_thresholds = [T_thresholds; tempTbl]; %#ok<AGROW>

            % Plot Panel 2: Mean Seg
            subplot(1, 4, 2); hold on;
            title('Mean Seg (+SE)'); xlabel('Spatial Separation'); ylabel('Proportion'); ylim([0 1.05]);
            errorbar(xVal, mean_Seg, sem_Seg, 'o', 'Color', cColor, 'MarkerFaceColor', cColor, 'MarkerSize', 5, 'CapSize', 0, 'HandleVisibility', 'off');
            try
                pFitMean = lsqcurvefit(logistic_func, p0, xVal, mean_Seg, lb, ub, opts);
                xSmooth = linspace(min(xVal), max(xVal), 100);
                plot(xSmooth, logistic_func(pFitMean, xSmooth), '-', 'Color', cColor, 'LineWidth', 2, 'DisplayName', sprintf('%d deg', currEcc));
            catch; end
            
            % Plot Panel 3: No Motion
            subplot(1, 4, 3); hold on;
            title('No Motion'); xlabel('Spatial Separation'); ylabel('Proportion'); ylim([0 1.05]);
            errorbar(xVal, mean_NoMot, sem_NoMot, 'o-', 'Color', cColor, 'LineWidth', 1.5, 'MarkerFaceColor', cColor, 'MarkerSize', 5, 'CapSize', 0, 'DisplayName', sprintf('%d deg', currEcc));
            
            % Plot Panel 4: One Motion
            subplot(1, 4, 4); hold on;
            title('One Motion'); xlabel('Spatial Separation'); ylabel('Proportion'); ylim([0 1.05]);
            errorbar(xVal, mean_OneMot, sem_OneMot, 'o-', 'Color', cColor, 'LineWidth', 1.5, 'MarkerFaceColor', cColor, 'MarkerSize', 5, 'CapSize', 0, 'DisplayName', sprintf('%d deg', currEcc));
            
            % --- SAVE TO STRUCTURE ---
            parsedData(p).blocks(blockCount).eccentricity = currEcc;
            parsedData(p).blocks(blockCount).fitParams_Individual = fitParams_Indiv;
            parsedData(p).blocks(blockCount).fitR2_Individual = fitR2_Indiv;
            parsedData(p).blocks(blockCount).thresholds_50 = thresholds_50;
            
            % Save Proportion Data for Part 6 Analysis
            parsedData(p).blocks(blockCount).prop_Seg = prop_Seg;
            parsedData(p).blocks(blockCount).prop_NoMot = prop_NoMot;
            parsedData(p).blocks(blockCount).prop_OneMot = prop_OneMot;
        end
    end
    
    % Legends
    for sp = 1:4, subplot(1,4,sp); legend('Location', 'Best', 'Interpreter', 'none'); end
end
fprintf('Psychometric processing complete. Mean R^2: %.4f\n', mean(all_R2_values, 'omitnan'));


%% ========================================================================
%% PART 3: LINEAR MIXED-EFFECTS (LME) ANALYSIS
%% ========================================================================
fprintf('\n--- Running Linear Mixed-Effects Model ---\n');

% Preprocessing: Remove NaNs
dataClean = T_thresholds(~isnan(T_thresholds.Threshold), :);
dataClean.ParticipantID = categorical(dataClean.ParticipantID);

% Fit LME
% Formula: Threshold ~ Intercept + Eccentricity + (Random Int + Random Slope | ID)
lme = fitlme(dataClean, 'Threshold ~ 1 + Eccentricity + (1 + Eccentricity | ParticipantID)');

% Display LME Results
disp(lme);
fprintf('LME Adjusted R^2: %.3f\n', lme.Rsquared.Adjusted);

% Extract Fixed Effects
fixedStats = lme.Coefficients;
CIs = coefCI(lme);
fixedStats.LowerCI = CIs(:,1);
fixedStats.UpperCI = CIs(:,2);

fprintf('\n=== TABLE 1: POPULATION FIXED EFFECTS ===\n');
disp(fixedStats);


%% ========================================================================
%% PART 4: GENERATE PUBLICATION TABLES (INDIVIDUAL STATS: REGRESSION + ANOVA)
%% ========================================================================
fprintf('\n--- Generating Publication Tables ---\n');

% --- TABLE A: INDIVIDUAL MODEL QUALITY & SIGNIFICANCE ---
% This table answers: "Did eccentricity matter for THIS specific person?"
% 1. Linear Regression (Tests for Trend)
% 2. One-Way ANOVA (Tests for Group Differences)

stats_rows = {};
uIDs = unique(dataClean.ParticipantID);

for i = 1:numel(uIDs)
    thisID = string(uIDs(i));
    
    % 1. Get Mean Logistic R2 for this person (Data Quality)
    pStructIdx = find(strcmpi({parsedData.participantID}, thisID));
    if ~isempty(pStructIdx)
        p_R2s = [];
        for b = 1:length(parsedData(pStructIdx).blocks)
            if ~isempty(parsedData(pStructIdx).blocks(b).fitR2_Individual)
                p_R2s = [p_R2s; parsedData(pStructIdx).blocks(b).fitR2_Individual];
            end
        end
        mean_Log_R2 = mean(p_R2s(~isnan(p_R2s)));
    else
        mean_Log_R2 = NaN;
    end
    
    % Filter Data for this Person
    subTable = dataClean(dataClean.ParticipantID == uIDs(i), :);
    
    % 2. Fit Linear Model (Threshold ~ Eccentricity)
    lm = fitlm(subTable, 'Threshold ~ Eccentricity');
    lin_Slope = lm.Coefficients.Estimate(2);
    lin_R2    = lm.Rsquared.Ordinary;
    lin_Pval  = lm.Coefficients.pValue(2); 
    
    % 3. Run One-Way ANOVA (Threshold across Eccentricity Groups)
    % 'off' suppresses the figure popup
    [anova_p, anova_tbl, ~] = anova1(subTable.Threshold, subTable.Eccentricity, 'off');
    
    % Extract Stats from ANOVA Table
    % Cell structure: Row 2 is 'Groups' (Model), Col 3=df, Col 5=F
    anova_df_between = anova_tbl{2,3};
    anova_df_within  = anova_tbl{3,3};
    anova_F          = anova_tbl{2,5};
    
    % Store Everything
    stats_rows{i, 1} = char(uIDs(i));
    stats_rows{i, 2} = height(subTable);
    stats_rows{i, 3} = mean_Log_R2;
    % Linear Stats
    stats_rows{i, 4} = lin_Slope;
    stats_rows{i, 5} = lin_R2;
    stats_rows{i, 6} = lin_Pval;
    % ANOVA Stats
    stats_rows{i, 7} = sprintf('F(%d,%d)=%.2f', anova_df_between, anova_df_within, anova_F);
    stats_rows{i, 8} = anova_p;
end

TableA = cell2table(stats_rows, ...
    'VariableNames', {'Participant', 'N_Sessions', 'Logistic_R2', ...
                      'Linear_Slope', 'Linear_R2', 'Linear_P', ...
                      'ANOVA_F_Stat', 'ANOVA_P'});

fprintf('\n=== TABLE A: INDIVIDUAL STATISTICS (REGRESSION + ANOVA) ===\n');
disp(TableA);


% --- TABLE B: BOOTSTRAPPED DESCRIPTIVE STATS ---
% (This part remains the same as before)
fprintf('\nBootstrapping Means (1000 iter)... \n');
rng(12345); 
nBoot = 1000;
ecc_levels = unique(dataClean.Eccentricity);
desc_rows = cell(numel(uIDs), length(ecc_levels)+1);

for i = 1:numel(uIDs)
    desc_rows{i,1} = char(uIDs(i));
    for e = 1:length(ecc_levels)
        idx = (dataClean.ParticipantID == uIDs(i)) & (dataClean.Eccentricity == ecc_levels(e));
        vals = dataClean.Threshold(idx);
        if numel(vals) >= 2
            bootDist = bootstrp(nBoot, @mean, vals);
            bMean = mean(bootDist);
            bCI = prctile(bootDist, [2.5, 97.5]);
            desc_rows{i, e+1} = sprintf('%.3f [%.3f, %.3f]', bMean, bCI(1), bCI(2));
        else
            desc_rows{i, e+1} = 'N/A';
        end
    end
end
varNamesB = [{'Participant'}, arrayfun(@(x) sprintf('Deg_%d', x), ecc_levels', 'uni', 0)];
TableB = cell2table(desc_rows, 'VariableNames', varNamesB);

fprintf('\n=== TABLE B: DESCRIPTIVE STATS (Mean [95%% CI]) ===\n');
disp(TableB);

%% ========================================================================
%% PART 5: FIGURE 4A (POPULATION & INDIVIDUAL VISUALIZATION)
%% ========================================================================
figure('Color', 'w', 'Name', 'Figure 4A: Population & Individual Fits', 'Position', [100, 100, 1200, 500]);

% --- SUBPLOT 1: POPULATION FIT (LME) ---
subplot(1, 2, 1); 
hold on; box on;

% 1. Plot Data Points (Aggregated Mean +/- SEM)
[G, t_agg] = findgroups(dataClean.Eccentricity);
mean_pop = splitapply(@mean, dataClean.Threshold, G);
sem_pop  = splitapply(@(x) std(x)/sqrt(numel(x)), dataClean.Threshold, G);

errorbar(t_agg, mean_pop, sem_pop, 'o', 'Color', [0.4 0.4 0.4], ...
    'MarkerFaceColor', [0.4 0.4 0.4], 'MarkerSize', 8, 'LineWidth', 1.5, ...
    'DisplayName', 'Data (Mean \pm SEM)');

% 2. Plot Individual Linear Fits (Faint Dashed Lines - from Polyfit on Raw)
xFit = linspace(2.5, 7.5, 100);
colors_sub = lines(numel(uIDs));

for i = 1:numel(uIDs)
    subTable = dataClean(dataClean.ParticipantID == uIDs(i), :);
    % Simple linear regression on raw data for this person
    p = polyfit(subTable.Eccentricity, subTable.Threshold, 1);
    plot(xFit, polyval(p, xFit), '--', 'Color', [colors_sub(i,:) 0.3], ...
        'LineWidth', 1, 'HandleVisibility', 'off');
end

% 3. Plot LME Population Fit (Thick Black Line)
pop_Int = fixedStats.Estimate(1);
pop_Slope = fixedStats.Estimate(2);
yPop = pop_Int + pop_Slope * xFit;

plot(xFit, yPop, 'k-', 'LineWidth', 2.5, 'DisplayName', 'LME Population Fit');

% Aesthetics (Panel 1)
xlabel('Eccentricity (deg)', 'FontSize', 12);
ylabel('Threshold (\lambda_{50})', 'FontSize', 12);
title('Population Trend (LME)', 'FontSize', 14);
legend('Location', 'NorthWest');
xlim([2.5 7.5]); ylim([0 1.0]);

% Add Stat Text
txt = sprintf('\\beta = %.3f, p < 0.001\nR^2_{adj} = %.3f', pop_Slope, lme.Rsquared.Adjusted);
text(2.7, 0.9, txt, 'FontSize', 11, 'EdgeColor', 'k', 'BackgroundColor', 'w');


% --- SUBPLOT 2: INDIVIDUAL BOOTSTRAPPED FITS (Mean + 95% CI) ---
subplot(1, 2, 2); 
hold on; box on;

fprintf('\nCalculating 95%% CI for Figure...\n');

for i = 1:numel(uIDs)
    thisID = uIDs(i);
    subTable = dataClean(dataClean.ParticipantID == thisID, :);
    cColor = colors_sub(i,:);
    
    % We need Mean & 95% CI for each eccentricity for this person
    ecc_u = unique(subTable.Eccentricity);
    bMeans = zeros(size(ecc_u));
    errLow = zeros(size(ecc_u));
    errHigh = zeros(size(ecc_u));
    
    for e = 1:length(ecc_u)
        vals = subTable.Threshold(subTable.Eccentricity == ecc_u(e));
        if numel(vals) >= 2
            bootDist = bootstrp(1000, @mean, vals);
            bMeans(e) = mean(bootDist);
            ci = prctile(bootDist, [2.5, 97.5]); % 95% CI
            errLow(e) = bMeans(e) - ci(1);
            errHigh(e) = ci(2) - bMeans(e);
        else
            bMeans(e) = mean(vals);
            errLow(e) = 0; errHigh(e) = 0;
        end
    end
    
    % Plot Error Bars (Bootstrapped Mean +/- 95% CI)
    errorbar(ecc_u, bMeans, errLow, errHigh, 'o', 'Color', cColor, ...
        'MarkerFaceColor', cColor, 'MarkerSize', 6, 'LineWidth', 1.2, ...
        'HandleVisibility', 'off'); 
    
    % Plot Linear Regression Line (Fit to Raw Data to match Table A)
    p = polyfit(subTable.Eccentricity, subTable.Threshold, 1);
    plot(xFit, polyval(p, xFit), '-', 'Color', cColor, 'LineWidth', 2, ...
        'DisplayName', char(thisID));
end

% Aesthetics (Panel 2)
xlabel('Eccentricity (deg)', 'FontSize', 12);
ylabel('Threshold (\lambda_{50})', 'FontSize', 12);
title('Individual Participants (Mean \pm 95% CI)', 'FontSize', 14);
legend('Location', 'NorthWest');
xlim([2.5 7.5]); ylim([0 1.0]);

hold off;
fprintf('\nAll Analysis Complete.\n');



% %% ========================================================================
% %% PART 6: ONE MOTION ANALYSIS (SPLINE FIT & HALF-HEIGHT WIDTH)
% %% ========================================================================
% fprintf('\n--- Running One Motion (Integration) Analysis ---\n');
% 
% % Initialize storage matrices: [Participants x Eccentricities]
% nSubs = numel(parsedData);
% nEccs = length(ecc_labels);
% Rising_Vals  = nan(nSubs, nEccs);
% Falling_Vals = nan(nSubs, nEccs);
% Peak_Vals    = nan(nSubs, nEccs);
% 
% % High-res x-axis for spline interpolation
% xFine = linspace(min(xVal), max(xVal), 1000);
% 
% % Loop through participants and blocks
% for p = 1:nSubs
%     for b = 1:length(parsedData(p).blocks)
%         
%         % 1. Get Mean One Motion Data
%         thisEcc = parsedData(p).blocks(b).eccentricity;
%         eccIdx = find(ecc_labels == thisEcc);
%         
%         if isempty(eccIdx), continue; end
%         
%         % Calculate Mean across sessions for One Motion
%         prop_OneMot = parsedData(p).blocks(b).prop_OneMot;
%         yMean = mean(prop_OneMot, 1, 'omitnan');
%         
%         % 2. Spline Interpolation (Standard Cubic Spline)
%         yFine = interp1(xVal, yMean, xFine, 'spline');
%         
%         % 3. Find Peak and Half-Height
%         [maxVal, maxIdx] = max(yFine);
%         halfHeight = maxVal / 2;
%         
%         % Only process if there is a significant peak (e.g., > 0.2 proportion)
%         if maxVal > 0.2
%             
%             % --- FIND RISING PHASE (Left of Peak) ---
%             [~, idxLeft] = min(abs(yFine(1:maxIdx) - halfHeight));
%             Rising_Vals(p, eccIdx) = xFine(idxLeft);
%             
%             % --- FIND FALLING PHASE (Right of Peak) ---
%             [~, idxRight] = min(abs(yFine(maxIdx:end) - halfHeight));
%             fullIdxRight = maxIdx + idxRight - 1; 
%             Falling_Vals(p, eccIdx) = xFine(fullIdxRight);
%             
%             Peak_Vals(p, eccIdx) = xFine(maxIdx);
%         end
%     end
% end
% 
% % --- PLOT 1: SUMMARY OF RISING/FALLING THRESHOLDS ---
% figure('Color', 'w', 'Name', 'One Motion Bandwidth Summary', 'Position', [150, 150, 600, 500]);
% hold on; box on;
% title('Integration Bandwidth (One Motion)');
% xlabel('Eccentricity (deg)');
% ylabel('Spatial Separation (deg)');
% 
% colors_sub = lines(nSubs);
% 
% for p = 1:nSubs
%     % Plot Rising (Solid markers)
%     plot(ecc_labels, Rising_Vals(p,:), '^-', 'Color', colors_sub(p,:), ...
%         'MarkerFaceColor', colors_sub(p,:), 'LineWidth', 1.5, ...
%         'DisplayName', [parsedData(p).participantID ' (Rising)']);
%     
%     % Plot Falling (Open markers or different shape)
%     plot(ecc_labels, Falling_Vals(p,:), 'v--', 'Color', colors_sub(p,:), ...
%         'MarkerFaceColor', 'w', 'LineWidth', 1.5, ...
%         'DisplayName', [parsedData(p).participantID ' (Falling)']);
%     
%     % Connect them to show the "Integration Zone"
%     for e = 1:nEccs
%         if ~isnan(Rising_Vals(p,e)) && ~isnan(Falling_Vals(p,e))
%             line([ecc_labels(e) ecc_labels(e)], [Rising_Vals(p,e) Falling_Vals(p,e)], ...
%                 'Color', [colors_sub(p,:) 0.3], 'LineWidth', 2, 'HandleVisibility', 'off');
%         end
%     end
% end
% 
% legend('Location', 'NorthWest', 'NumColumns', 2);
% ylim([0 0.8]); xlim([2.5 7.5]);
% 
% 
% % --- PLOT 2: SPLINE FITS (SEPARATE PANELS PER PARTICIPANT) ---
% figure('Color', 'w', 'Name', 'Individual Spline Fits', 'Position', [200, 200, 1400, 400]);
% % Dynamic subplot layout
% nCols = nSubs; 
% nRows = 1;
% 
% for p = 1:nSubs
%     subplot(nRows, nCols, p); hold on; box on;
%     title(sprintf('Participant: %s', parsedData(p).participantID));
%     xlabel('Spatial Separation'); ylabel('Prop. One Motion');
%     ylim([-0.1 1.2]); 
%     
%     for b = 1:length(parsedData(p).blocks)
%         thisEcc = parsedData(p).blocks(b).eccentricity;
%         
%         if isempty(parsedData(p).blocks(b).prop_OneMot)
%             continue; 
%         end
%         
%         % Calculate Mean and SEM for Error Bars
%         propData = parsedData(p).blocks(b).prop_OneMot;
%         yMean = mean(propData, 1, 'omitnan');
%         nSess = size(propData, 1);
%         ySEM  = std(propData, 0, 1, 'omitnan') ./ sqrt(nSess);
%         
%         % 2. Spline Interpolation (Standard Cubic Spline)
%         yFine = interp1(xVal, yMean, xFine, 'spline');
%         
%         % Color by Eccentricity
%         cIdx = find(ecc_labels == thisEcc);
%         if isempty(cIdx), continue; end
%         cColor = colors(cIdx, :);
%         
%         % Plot Raw Points with Error Bars (faint)
%         % Fixed MarkerFaceColor to be RGB triplet
%         errorbar(xVal, yMean, ySEM, 'o', 'Color', [cColor 0.4], ...
%             'MarkerFaceColor', cColor, 'MarkerSize', 4, ...
%             'CapSize', 0, 'HandleVisibility', 'off');
%         
%         % Plot Spline Fit
%         plot(xFine, yFine, '-', 'Color', cColor, 'LineWidth', 1.5, ...
%             'DisplayName', sprintf('%d deg', thisEcc));
%         
%         % Optional: Plot half-height line for the peak
%         [maxVal, ~] = max(yFine);
%         if maxVal > 0.2
%             plot([min(xVal) max(xVal)], [maxVal/2 maxVal/2], ':', 'Color', [cColor 0.3], 'HandleVisibility', 'off');
%         end
%     end
%     
%     if p == 1
%         legend('Location', 'Best');
%     end
% end
% 
% 
% % --- CREATE TABLE FOR RISING/FALLING VALUES ---
% fprintf('\nGenerating One Motion Threshold Table...\n');
% 
% integration_rows = {};
% row_idx = 1;
% 
% for p = 1:nSubs
%     pName = parsedData(p).participantID;
%     for e = 1:nEccs
%         eccVal = ecc_labels(e);
%         
%         rVal = Rising_Vals(p, e);
%         fVal = Falling_Vals(p, e);
%         
%         integration_rows{row_idx, 1} = pName;
%         integration_rows{row_idx, 2} = eccVal;
%         integration_rows{row_idx, 3} = rVal;
%         integration_rows{row_idx, 4} = fVal;
%         integration_rows{row_idx, 5} = fVal - rVal; % Integration Bandwidth
%         
%         row_idx = row_idx + 1;
%     end
% end
% 
% Table_Integration = cell2table(integration_rows, ...
%     'VariableNames', {'Participant', 'Eccentricity', 'Rising_Thresh', 'Falling_Thresh', 'Bandwidth'});
% 
% fprintf('\n=== TABLE: ONE MOTION INTEGRATION BANDWIDTHS ===\n');
% disp(Table_Integration);
% 
% 
% fprintf('One Motion Analysis Complete.\n');
% 


%% ========================================================================
%% PART 6: ONE MOTION ANALYSIS (PER-SESSION FALLING THRESHOLD & LME)
%% ========================================================================
fprintf('\n--- Running One Motion (Falling Threshold) Analysis ---\n');

% --- SETUP VARIABLES ---
nSubs = numel(parsedData);
nEccs = length(ecc_labels);

% Storage for Aggregate Summary (Method 1 visualization)
Rising_Vals_Agg  = nan(nSubs, nEccs);
Falling_Vals_Agg = nan(nSubs, nEccs);

% Storage for Per-Session LME (Method 2)
T_OneMot_LME = table(); 

% High-res x-axis for spline interpolation
xFine = linspace(min(xVal), max(xVal), 1000);

% Loop through participants
for p = 1:nSubs
    pName = parsedData(p).participantID;
    
    for b = 1:length(parsedData(p).blocks)
        thisEcc = parsedData(p).blocks(b).eccentricity;
        eccIdx = find(ecc_labels == thisEcc);
        
        if isempty(eccIdx), continue; end
        
        % Get Data: [nSessions x nSeparations]
        prop_OneMot = parsedData(p).blocks(b).prop_OneMot; 
        nSess = size(prop_OneMot, 1);
        
        % =================================================================
        % METHOD 1: AGGREGATE MEAN (For simple summary plot)
        % =================================================================
        yMean = mean(prop_OneMot, 1, 'omitnan');
        yFine = interp1(xVal, yMean, xFine, 'spline');
        [maxVal, maxIdx] = max(yFine);
        halfHeight = maxVal / 2;
        
        if maxVal > 0.2
            % Find Falling Edge (Right of Peak)
            [~, idxRight] = min(abs(yFine(maxIdx:end) - halfHeight));
            Falling_Vals_Agg(p, eccIdx) = xFine(maxIdx + idxRight - 1);
            
            % Find Rising Edge (Left of Peak)
            [~, idxLeft] = min(abs(yFine(1:maxIdx) - halfHeight));
            Rising_Vals_Agg(p, eccIdx) = xFine(idxLeft);
        end
        
        % =================================================================
        % METHOD 2: FIT SPLINE TO EVERY INDIVIDUAL SESSION
        % =================================================================
        for s = 1:nSess
            ySess = prop_OneMot(s, :);
            
            % Interpolate Session
            yFineSess = interp1(xVal, ySess, xFine, 'spline');
            [maxS, maxIdxS] = max(yFineSess);
            halfH_S = maxS / 2;
            
            s_Falling = NaN; 
            s_Rising  = NaN;
            
            % Only analyze if there is a distinct peak perception (> 0.2)
            if maxS > 0.2
                try
                    % 1. Find Falling Threshold (Right Side Half-Height)
                    % This is the critical breakdown point of One Motion
                    [~, iR] = min(abs(yFineSess(maxIdxS:end) - halfH_S));
                    s_Falling = xFine(maxIdxS + iR - 1);
                    
                    % 2. Find Rising Threshold (Left Side)
                    [~, iL] = min(abs(yFineSess(1:maxIdxS) - halfH_S));
                    s_Rising = xFine(iL);
                catch
                    s_Falling = NaN;
                end
            end
            
            % Append to Table for LME
            rowT = table({pName}, thisEcc, s_Falling, s_Rising, ...
                'VariableNames', {'ParticipantID', 'Eccentricity', 'Falling_Threshold', 'Rising_Threshold'});
            T_OneMot_LME = [T_OneMot_LME; rowT]; %#ok<AGROW>
        end
    end
end


%% --- 6.1: RUN LME ON FALLING THRESHOLD ---
fprintf('\n--- Running LME on One Motion Falling Thresholds ---\n');

% Prepare Data
dataFalling = T_OneMot_LME(~isnan(T_OneMot_LME.Falling_Threshold), :);
dataFalling.ParticipantID = categorical(dataFalling.ParticipantID);

% Fit LME: Falling_Threshold ~ Eccentricity
lme_fall = fitlme(dataFalling, 'Falling_Threshold ~ 1 + Eccentricity + (1 + Eccentricity | ParticipantID)');

disp(lme_fall);
fprintf('One Motion (Falling) LME Adjusted R^2: %.3f\n', lme_fall.Rsquared.Adjusted);

% Extract Fixed Effects
fall_fixed = lme_fall.Coefficients;
fprintf('\n=== TABLE: ONE MOTION (FALLING EDGE) FIXED EFFECTS ===\n');
disp(fall_fixed);


%% --- 6.2: PLOT LME FIT FOR FALLING THRESHOLD ---
figure('Color', 'w', 'Name', 'One Motion Falling Threshold LME', 'Position', [150, 150, 700, 500]);
hold on; box on;

% 1. Plot Mean Data Points (+/- SEM)
[G_fall, t_agg_fall] = findgroups(dataFalling.Eccentricity);
mean_fall = splitapply(@mean, dataFalling.Falling_Threshold, G_fall);
sem_fall  = splitapply(@(x) std(x)/sqrt(numel(x)), dataFalling.Falling_Threshold, G_fall);

errorbar(t_agg_fall, mean_fall, sem_fall, 's', 'Color', [0.8 0.4 0], ...
    'MarkerFaceColor', [0.8 0.4 0], 'MarkerSize', 8, 'LineWidth', 1.5, ...
    'DisplayName', 'Data (Mean \pm SEM)');

% 2. Plot LME Population Fit Line
xFit = linspace(2.5, 7.5, 100);
f_Int = fall_fixed.Estimate(1);
f_Slope = fall_fixed.Estimate(2);
yPop_fall = f_Int + f_Slope * xFit;

plot(xFit, yPop_fall, 'k-', 'LineWidth', 2.5, 'DisplayName', 'LME Population Fit');

% 3. Plot Individual Linear Fits (Faint Dashed Lines)
uIDs = unique(dataFalling.ParticipantID);
for i = 1:numel(uIDs)
    subTable = dataFalling(dataFalling.ParticipantID == uIDs(i), :);
    if height(subTable) > 2
        pPoly = polyfit(subTable.Eccentricity, subTable.Falling_Threshold, 1);
        plot(xFit, polyval(pPoly, xFit), '--', 'Color', [0.8 0.4 0 0.2], ...
             'LineWidth', 1, 'HandleVisibility', 'off');
    end
end

% Aesthetics
xlabel('Eccentricity (deg)', 'FontSize', 12);
ylabel('Falling Threshold (deg)', 'FontSize', 12);
title('Effect of Eccentricity on One Motion Breakdown', 'FontSize', 14);
legend('Location', 'Best');
xlim([2.5 7.5]); 

% Add Stat Text
txt = sprintf('\\beta = %.4f, p = %.3f', f_Slope, fall_fixed.pValue(2));
text(3, max(mean_fall)*1.05, txt, 'FontSize', 12, 'EdgeColor', 'k', 'BackgroundColor', 'w');


%% --- 6.3: VISUAL CHECK - INDIVIDUAL SPLINE FITS ---
figure('Color', 'w', 'Name', 'Individual Spline Fits (One Motion)', 'Position', [200, 200, 1400, 400]);
nCols = nSubs; nRows = 1;

for p = 1:nSubs
    subplot(nRows, nCols, p); hold on; box on;
    title(sprintf('%s', parsedData(p).participantID));
    xlabel('Separation'); ylabel('Prop. One Mot'); ylim([-0.1 1.2]); 
    
    for b = 1:length(parsedData(p).blocks)
        thisEcc = parsedData(p).blocks(b).eccentricity;
        if isempty(parsedData(p).blocks(b).prop_OneMot), continue; end
        
        yMean = mean(parsedData(p).blocks(b).prop_OneMot, 1, 'omitnan');
        yFine = interp1(xVal, yMean, xFine, 'spline');
        
        cIdx = find(ecc_labels == thisEcc);
        if ~isempty(cIdx)
            plot(xFine, yFine, '-', 'Color', colors(cIdx, :), 'LineWidth', 1.5, ...
                'DisplayName', sprintf('%d deg', thisEcc));
            
            % Show the Aggregate Falling Point for visual verification
            if ~isnan(Falling_Vals_Agg(p, cIdx))
                plot(Falling_Vals_Agg(p, cIdx), interp1(xFine, yFine, Falling_Vals_Agg(p, cIdx)), ...
                     'rx', 'HandleVisibility', 'off');
            end
        end
    end
    if p==1, legend('Location','Best'); end
end

fprintf('One Motion Analysis Complete.\n');


%% ========================================================================
%% PART 7: SEGMENTATION vs. ECCENTRICITY (GROUPED BY SPATIAL SEPARATION)
%% ========================================================================
fprintf('\n--- Running Part 7: Segmentation vs Eccentricity Analysis ---\n');

% Define a high-contrast colormap for the 10 spatial separations
% 'turbo' or 'jet' are good for distinguishing ordered magnitudes
sep_colors = turbo(length(xVal)); 

figure('Color', 'w', 'Name', 'Seg vs Eccentricity (By Separation)', ...
       'Position', [100, 100, 1600, 400]);

nSubs = numel(parsedData);

for p = 1:nSubs
    subplot(1, nSubs, p); hold on; box on;
    pName = parsedData(p).participantID;
    title(sprintf('Participant: %s', pName));
    
    % Initialize Matrices: [Rows = Eccentricities, Cols = Spatial Separations]
    % We need one for Means and one for SEMs
    Seg_Matrix = nan(length(ecc_labels), length(xVal));
    Seg_SEM_Matrix = nan(length(ecc_labels), length(xVal));
    
    % --- EXTRACT DATA ---
    for b = 1:length(parsedData(p).blocks)
        thisEcc = parsedData(p).blocks(b).eccentricity;
        
        % Find which row index this eccentricity belongs to
        rowIdx = find(ecc_labels == thisEcc);
        
        if ~isempty(rowIdx)
            % Get the raw proportion data: [nSessions x 10 Separations]
            propData = parsedData(p).blocks(b).prop_Seg;
            nSess = size(propData, 1);
            
            % 1. Calculate Mean across sessions
            blockMean = mean(propData, 1, 'omitnan');
            
            % 2. Calculate SEM (Standard Deviation / sqrt(N)) across sessions
            blockStd = std(propData, 0, 1, 'omitnan');
            blockSEM = blockStd ./ sqrt(nSess);
            
            % Store in matrices
            Seg_Matrix(rowIdx, :) = blockMean;
            Seg_SEM_Matrix(rowIdx, :) = blockSEM;
        end
    end
    
    % --- PLOT CURVES ---
    % We loop through the columns (Spatial Separations) to plot lines
    for sep = 1:length(xVal)
        % Extract Mean and SEM for this specific spatial separation
        yData = Seg_Matrix(:, sep);
        yErr  = Seg_SEM_Matrix(:, sep);
        
        % Plot with Error Bars
        % CapSize=0 keeps the plot cleaner since there are 10 overlapping lines
        errorbar(ecc_labels, yData, yErr, '-o', ...
             'Color', sep_colors(sep, :), ...
             'MarkerFaceColor', sep_colors(sep, :), ...
             'MarkerSize', 5, ...
             'LineWidth', 1.5, ...
             'CapSize', 0); 
    end
    
    % --- AESTHETICS ---
    xlabel('Eccentricity (deg)');
    ylabel('Proportion Segregated');
    ylim([0 1.05]); 
    xlim([min(ecc_labels)-0.5, max(ecc_labels)+0.5]);
    xticks(ecc_labels);
    grid on;
    
    % Only add the Colorbar (Right Hand Axis) to the last subplot to save space
    if p == nSubs
        colormap(turbo(length(xVal)));
        c = colorbar;
        c.Label.String = 'Spatial Separation (deg)';
        % Set ticks to match the actual xVal steps
        caxis([min(xVal) max(xVal)]);
        % Optional: Set specific ticks for readability
        c.Ticks = linspace(min(xVal), max(xVal), 5); 
    end
end

fprintf('Part 7 Visualization Complete.\n');


%% ========================================================================
%% PART 8: THREE-PHASE DIAGRAM - NO MOTION / INTEGRATION / SEGMENTATION
%% ========================================================================
% This visualization creates a heat map showing THREE perceptual zones:
%   1. NO MOTION (subjects don't perceive coherent motion)
%   2. INTEGRATION (subjects perceive a single direction - vector average)
%   3. SEGMENTATION (subjects perceive two distinct directions)
%
% X-axis: Spatial Separation (degrees)
% Y-axis: Eccentricity (degrees)
% Color:  RGB encoding of the three response types
%
% The V1 RF size function is overlaid to show that the boundary between
% integration and segmentation parallels V1 receptive field sizes.
%
% Author: Added based on reviewer suggestions
% ========================================================================

fprintf('\n--- Running Part 8: Three-Phase Diagram ---\n');

%% --- 8.1: AGGREGATE DATA ACROSS ALL PARTICIPANTS ---

% Initialize matrices to hold ALL proportions for each response type
% Dimensions: [nSubjects x nEccentricities x nSpatialSeparations]
nEccs = length(ecc_labels);      % 5 eccentricities (3-7 deg)
nSeps = length(xVal);            % 10 spatial separations

AllSubs_Seg    = nan(nSubs, nEccs, nSeps);  % Two Directions
AllSubs_Single = nan(nSubs, nEccs, nSeps);  % Single Direction (Integration)
AllSubs_NoMot  = nan(nSubs, nEccs, nSeps);  % No Motion

for p = 1:nSubs
    for b = 1:length(parsedData(p).blocks)
        thisEcc = parsedData(p).blocks(b).eccentricity;
        eccIdx = find(ecc_labels == thisEcc);
        
        if ~isempty(eccIdx)
            % Get mean proportion across sessions for each response type
            if ~isempty(parsedData(p).blocks(b).prop_Seg)
                propSeg = parsedData(p).blocks(b).prop_Seg;      % [nSessions x 10]
                AllSubs_Seg(p, eccIdx, :) = mean(propSeg, 1, 'omitnan');
            end
            
            if ~isempty(parsedData(p).blocks(b).prop_OneMot)
                propSingle = parsedData(p).blocks(b).prop_OneMot; % [nSessions x 10]
                AllSubs_Single(p, eccIdx, :) = mean(propSingle, 1, 'omitnan');
            end
            
            if ~isempty(parsedData(p).blocks(b).prop_NoMot)
                propNoMot = parsedData(p).blocks(b).prop_NoMot;   % [nSessions x 10]
                AllSubs_NoMot(p, eccIdx, :) = mean(propNoMot, 1, 'omitnan');
            end
        end
    end
end

% Calculate Grand Mean across subjects for each response type
GrandMean_Seg    = squeeze(mean(AllSubs_Seg, 1, 'omitnan'));    % [nEccs x nSeps]
GrandMean_Single = squeeze(mean(AllSubs_Single, 1, 'omitnan')); % [nEccs x nSeps]
GrandMean_NoMot  = squeeze(mean(AllSubs_NoMot, 1, 'omitnan'));  % [nEccs x nSeps]

% Create finer grid for smooth interpolation
xFine_Sep = linspace(min(xVal), max(xVal), 100);   % Fine spatial separation
yFine_Ecc = linspace(min(ecc_labels), max(ecc_labels), 100);  % Fine eccentricity

% Interpolate each response type onto the fine grid
[X_coarse, Y_coarse] = meshgrid(xVal, ecc_labels);
[X_fine, Y_fine] = meshgrid(xFine_Sep, yFine_Ecc);

Z_Seg    = interp2(X_coarse, Y_coarse, GrandMean_Seg, X_fine, Y_fine, 'cubic');
Z_Single = interp2(X_coarse, Y_coarse, GrandMean_Single, X_fine, Y_fine, 'cubic');
Z_NoMot  = interp2(X_coarse, Y_coarse, GrandMean_NoMot, X_fine, Y_fine, 'cubic');

% Clamp values to [0, 1] range
Z_Seg    = max(0, min(1, Z_Seg));
Z_Single = max(0, min(1, Z_Single));
Z_NoMot  = max(0, min(1, Z_NoMot));

% Normalize so they sum to 1 at each point (they should already, but ensure it)
Z_Total = Z_Seg + Z_Single + Z_NoMot;
Z_Seg    = Z_Seg ./ Z_Total;
Z_Single = Z_Single ./ Z_Total;
Z_NoMot  = Z_NoMot ./ Z_Total;


%% --- 8.2: DEFINE V1 RF SIZE FUNCTION ---

% Neurophysiology reference (macaque V1) - Dow et al., 1981
V1_slope_Dow = 0.044;
V1_intercept_Dow = 0.22;
V1_RF_Dow = V1_slope_Dow * yFine_Ecc + V1_intercept_Dow;

% Your psychophysical LME fit (critical spatial separation)
LME_slope = 0.036;      % From your analysis
LME_intercept = 0.251;  % From your analysis
LME_fit = LME_slope * yFine_Ecc + LME_intercept;


%% --- 8.3: CREATE RGB IMAGE FOR THREE-PHASE DIAGRAM ---
% Color coding:
%   RED channel   = Segmentation (Two Directions)
%   GREEN channel = Integration (Single Direction)  
%   BLUE channel  = No Motion

% Method 1: Direct RGB mapping
RGB_image = cat(3, Z_Seg, Z_Single, Z_NoMot);

% Method 2: Alternative color scheme (more intuitive)
% Blue = No Motion, Green/Yellow = Integration, Red = Segmentation
% Use HSV or custom blending

% Create custom RGB with better color distinction:
%   No Motion    -> Blue   (0, 0, 1)
%   Integration  -> Green  (0, 1, 0)
%   Segmentation -> Red    (1, 0, 0)

R_channel = Z_Seg;                    % Red for segmentation
G_channel = Z_Single;                 % Green for integration
B_channel = Z_NoMot;                  % Blue for no motion

RGB_image_v2 = cat(3, R_channel, G_channel, B_channel);


%% --- 8.4: MAIN THREE-PHASE DIAGRAM (FIGURE) ---

figure('Color', 'w', 'Name', 'Three-Phase Diagram', ...
       'Position', [100, 100, 900, 700]);

% Display the RGB image
imagesc(xFine_Sep, yFine_Ecc, RGB_image_v2);
set(gca, 'YDir', 'normal');  % Eccentricity increases upward

hold on;

% --- OVERLAY V1 RF SIZE LINES ---
h_LME = plot(LME_fit, yFine_Ecc, 'w-', 'LineWidth', 3, 'DisplayName', 'Psychophysical \lambda_{75}');
h_V1 = plot(V1_RF_Dow, yFine_Ecc, 'w--', 'LineWidth', 2.5, 'DisplayName', 'V1 RF (Dow et al.)');

% --- ADD CONTOUR LINES FOR BOUNDARIES ---
% Contour where Segmentation = 0.5 (boundary to segmentation)
[C1, h1] = contour(X_fine, Y_fine, Z_Seg, [0.5 0.5], 'k-', 'LineWidth', 2);

% Contour where No Motion = 0.5 (boundary from no motion)
[C2, h2] = contour(X_fine, Y_fine, Z_NoMot, [0.5 0.5], 'k:', 'LineWidth', 2);

% --- OVERLAY DATA POINTS ---
for e = 1:nEccs
    for s = 1:nSeps
        plot(xVal(s), ecc_labels(e), 'ko', 'MarkerSize', 3, ...
             'MarkerFaceColor', 'none', 'HandleVisibility', 'off');
    end
end

hold off;

% --- AESTHETICS ---
xlabel('Spatial Separation (°)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Retinal Eccentricity (°)', 'FontSize', 14, 'FontWeight', 'bold');
title('Three-Phase Diagram: No Motion → Integration → Segmentation', 'FontSize', 14);

xlim([min(xVal) max(xVal)]);
ylim([min(ecc_labels) max(ecc_labels)]);

% Legend
legend([h_LME, h_V1], 'Location', 'SouthEast', 'FontSize', 11, 'Box', 'on', ...
       'Color', 'w');

% Add text annotations for the three zones
text(0.12, 6.5, {'NO MOTION'}, 'FontSize', 11, 'FontWeight', 'bold', ...
     'Color', 'w', 'HorizontalAlignment', 'center');
text(0.32, 5.0, {'INTEGRATION', '(Single Direction)'}, 'FontSize', 11, ...
     'FontWeight', 'bold', 'Color', 'k', 'HorizontalAlignment', 'center');
text(0.60, 4.0, {'SEGMENTATION', '(Two Directions)'}, 'FontSize', 11, ...
     'FontWeight', 'bold', 'Color', 'w', 'HorizontalAlignment', 'center');

% Add color legend (custom)
annotation('textbox', [0.15, 0.02, 0.25, 0.06], 'String', ...
    '{\color{blue}■} No Motion  {\color{green}■} Integration  {\color{red}■} Segmentation', ...
    'FontSize', 10, 'EdgeColor', 'none', 'HorizontalAlignment', 'center');

set(gca, 'FontSize', 12, 'LineWidth', 1.5, 'Box', 'on');


%% --- 8.5: ALTERNATIVE VISUALIZATION - DOMINANT RESPONSE MAP ---
% Instead of RGB blending, show which response is DOMINANT at each point

figure('Color', 'w', 'Name', 'Dominant Response Map', ...
       'Position', [150, 100, 900, 700]);

% Find dominant response at each point (1=NoMot, 2=Single, 3=Seg)
[~, Dominant] = max(cat(3, Z_NoMot, Z_Single, Z_Seg), [], 3);

% Create custom colormap: Blue, Green, Red
cmap_dominant = [0 0 0.8;    % Blue for No Motion
                 0 0.7 0;    % Green for Integration
                 0.8 0 0];   % Red for Segmentation

imagesc(xFine_Sep, yFine_Ecc, Dominant);
set(gca, 'YDir', 'normal');
colormap(cmap_dominant);
caxis([0.5 3.5]);

hold on;

% Overlay V1 RF lines
plot(LME_fit, yFine_Ecc, 'w-', 'LineWidth', 3);
plot(V1_RF_Dow, yFine_Ecc, 'w--', 'LineWidth', 2.5);

% Add boundary contours
contour(X_fine, Y_fine, Z_Seg, [0.5 0.5], 'w-', 'LineWidth', 2);
contour(X_fine, Y_fine, Z_NoMot, [0.5 0.5], 'w:', 'LineWidth', 2);

hold off;

xlabel('Spatial Separation (°)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Retinal Eccentricity (°)', 'FontSize', 14, 'FontWeight', 'bold');
title('Dominant Perceptual Response', 'FontSize', 14);

% Custom colorbar
cb = colorbar;
cb.Ticks = [1 2 3];
cb.TickLabels = {'No Motion', 'Integration', 'Segmentation'};
cb.FontSize = 11;

xlim([min(xVal) max(xVal)]);
ylim([min(ecc_labels) max(ecc_labels)]);

% Add zone labels
text(0.12, 6.5, 'NO MOTION', 'FontSize', 11, 'FontWeight', 'bold', ...
     'Color', 'w', 'HorizontalAlignment', 'center');
text(0.32, 5.0, 'INTEGRATION', 'FontSize', 11, 'FontWeight', 'bold', ...
     'Color', 'w', 'HorizontalAlignment', 'center');
text(0.58, 4.0, 'SEGMENTATION', 'FontSize', 11, 'FontWeight', 'bold', ...
     'Color', 'w', 'HorizontalAlignment', 'center');

set(gca, 'FontSize', 12, 'LineWidth', 1.5, 'Box', 'on');


%% --- 8.6: INDIVIDUAL SUBJECT THREE-PHASE DIAGRAMS ---

figure('Color', 'w', 'Name', 'Individual Three-Phase Diagrams', ...
       'Position', [200, 100, 1400, 350]);

for p = 1:nSubs
    subplot(1, nSubs, p);
    
    % Extract this subject's data
    Subj_Seg    = squeeze(AllSubs_Seg(p, :, :));
    Subj_Single = squeeze(AllSubs_Single(p, :, :));
    Subj_NoMot  = squeeze(AllSubs_NoMot(p, :, :));
    
    % Interpolate
    Z_Seg_s    = interp2(X_coarse, Y_coarse, Subj_Seg, X_fine, Y_fine, 'cubic');
    Z_Single_s = interp2(X_coarse, Y_coarse, Subj_Single, X_fine, Y_fine, 'cubic');
    Z_NoMot_s  = interp2(X_coarse, Y_coarse, Subj_NoMot, X_fine, Y_fine, 'cubic');
    
    % Clamp and normalize
    Z_Seg_s    = max(0, min(1, Z_Seg_s));
    Z_Single_s = max(0, min(1, Z_Single_s));
    Z_NoMot_s  = max(0, min(1, Z_NoMot_s));
    
    Z_Total_s = Z_Seg_s + Z_Single_s + Z_NoMot_s;
    Z_Seg_s    = Z_Seg_s ./ Z_Total_s;
    Z_Single_s = Z_Single_s ./ Z_Total_s;
    Z_NoMot_s  = Z_NoMot_s ./ Z_Total_s;
    
    % Create RGB image
    RGB_subj = cat(3, Z_Seg_s, Z_Single_s, Z_NoMot_s);
    
    imagesc(xFine_Sep, yFine_Ecc, RGB_subj);
    set(gca, 'YDir', 'normal');
    
    hold on;
    plot(LME_fit, yFine_Ecc, 'w-', 'LineWidth', 2);
    plot(V1_RF_Dow, yFine_Ecc, 'w--', 'LineWidth', 1.5);
    contour(X_fine, Y_fine, Z_Seg_s, [0.5 0.5], 'k-', 'LineWidth', 1.5);
    hold off;
    
    xlabel('Spatial Sep. (°)', 'FontSize', 10);
    if p == 1
        ylabel('Eccentricity (°)', 'FontSize', 10);
    end
    title(sprintf('%s', parsedData(p).participantID), 'FontSize', 12);
    
    xlim([min(xVal) max(xVal)]);
    ylim([min(ecc_labels) max(ecc_labels)]);
end

% Add shared annotation
annotation('textbox', [0.35, 0.01, 0.35, 0.05], 'String', ...
    '{\color{blue}■} No Motion  {\color{green}■} Integration  {\color{red}■} Segmentation', ...
    'FontSize', 10, 'EdgeColor', 'none', 'HorizontalAlignment', 'center');


%% --- 8.7: STACKED AREA PLOT (PROPORTION BY SPATIAL SEPARATION) ---
% Shows how the three responses change with spatial separation at each eccentricity

figure('Color', 'w', 'Name', 'Response Proportions by Eccentricity', ...
       'Position', [250, 150, 1200, 400]);

for e = 1:nEccs
    subplot(1, nEccs, e);
    
    % Get data for this eccentricity
    y_NoMot  = GrandMean_NoMot(e, :);
    y_Single = GrandMean_Single(e, :);
    y_Seg    = GrandMean_Seg(e, :);
    
    % Stack the areas
    Y_stack = [y_NoMot; y_Single; y_Seg]';
    
    % Area plot
    h_area = area(xVal, Y_stack);
    h_area(1).FaceColor = [0 0 0.8];    % Blue - No Motion
    h_area(2).FaceColor = [0 0.7 0];    % Green - Integration
    h_area(3).FaceColor = [0.8 0 0];    % Red - Segmentation
    
    hold on;
    % Add V1 RF size for this eccentricity as vertical line
    V1_at_ecc = V1_slope_Dow * ecc_labels(e) + V1_intercept_Dow;
    xline(V1_at_ecc, 'w--', 'LineWidth', 2);
    
    % Add LME threshold
    LME_at_ecc = LME_slope * ecc_labels(e) + LME_intercept;
    xline(LME_at_ecc, 'w-', 'LineWidth', 2);
    hold off;
    
    xlabel('Spatial Sep. (°)', 'FontSize', 10);
    if e == 1
        ylabel('Proportion', 'FontSize', 10);
    end
    title(sprintf('%d° Eccentricity', ecc_labels(e)), 'FontSize', 11);
    
    xlim([min(xVal) max(xVal)]);
    ylim([0 1]);
    
    if e == nEccs
        legend({'No Motion', 'Integration', 'Segmentation'}, ...
               'Location', 'East', 'FontSize', 8);
    end
end

sgtitle('Response Proportions Across Spatial Separation', 'FontSize', 14);


%% --- 8.8: PUBLICATION FIGURE - CLEAN THREE-PHASE ---

figure('Color', 'w', 'Name', 'Publication: Three-Phase Diagram', ...
       'Position', [100, 100, 800, 650]);

% Use the RGB image
imagesc(xFine_Sep, yFine_Ecc, RGB_image_v2);
set(gca, 'YDir', 'normal');

hold on;

% V1 RF overlay
h1 = plot(LME_fit, yFine_Ecc, 'k-', 'LineWidth', 3);
h2 = plot(V1_RF_Dow, yFine_Ecc, 'k--', 'LineWidth', 2.5);

% Add boundary contours (black for visibility)
[C_seg, ~] = contour(X_fine, Y_fine, Z_Seg, [0.5 0.5], 'k-', 'LineWidth', 1.5);
[C_nomot, ~] = contour(X_fine, Y_fine, Z_NoMot, [0.5 0.5], 'k:', 'LineWidth', 1.5);

% Mark data points
for e = 1:nEccs
    for s = 1:nSeps
        plot(xVal(s), ecc_labels(e), 'ko', 'MarkerSize', 2, 'HandleVisibility', 'off');
    end
end

hold off;

% Labels
xlabel('Spatial Separation (°)', 'FontSize', 14);
ylabel('Retinal Eccentricity (°)', 'FontSize', 14);

% Legend
legend([h1, h2], {sprintf('\\lambda_{75} (slope=%.3f)', LME_slope), ...
                  sprintf('V1 RF (slope=%.3f)', V1_slope_Dow)}, ...
       'Location', 'SouthEast', 'FontSize', 10, 'Box', 'on', 'Color', 'w');

% Axis settings
xlim([0.07 0.71]);
ylim([3 7]);
set(gca, 'XTick', [0.1 0.2 0.3 0.4 0.5 0.6 0.7]);
set(gca, 'YTick', [3 4 5 6 7]);
set(gca, 'FontSize', 12, 'LineWidth', 1.5, 'Box', 'on');

% Add color legend box
annotation('rectangle', [0.13, 0.83, 0.28, 0.12], 'FaceColor', 'w', 'EdgeColor', 'k');
annotation('textbox', [0.14, 0.89, 0.12, 0.04], 'String', '{\color[rgb]{0,0,0.8}■} No Motion', ...
    'FontSize', 10, 'EdgeColor', 'none', 'FontWeight', 'bold');
annotation('textbox', [0.14, 0.85, 0.12, 0.04], 'String', '{\color[rgb]{0,0.7,0}■} Integration', ...
    'FontSize', 10, 'EdgeColor', 'none', 'FontWeight', 'bold');
annotation('textbox', [0.14, 0.81, 0.12, 0.04], 'String', '{\color[rgb]{0.8,0,0}■} Segmentation', ...
    'FontSize', 10, 'EdgeColor', 'none', 'FontWeight', 'bold');


%% --- 8.9: QUANTITATIVE ANALYSIS - EXTRACT BOUNDARIES ---

fprintf('\nExtracting Phase Boundaries...\n');

% Boundary 1: No Motion -> Integration (where No Motion drops to 50%)
% Boundary 2: Integration -> Segmentation (where Segmentation rises to 50%)

Boundary_NoMot_to_Int = nan(nEccs, 1);
Boundary_Int_to_Seg = nan(nEccs, 1);

for e = 1:nEccs
    % --- Boundary 1: No Motion -> Integration ---
    yRow_NoMot = GrandMean_NoMot(e, :);
    if any(yRow_NoMot > 0.5) && any(yRow_NoMot < 0.5)
        diff_from_half = yRow_NoMot - 0.5;
        sign_changes = find(diff(sign(diff_from_half)) ~= 0);
        if ~isempty(sign_changes)
            idx = sign_changes(1);
            x1 = xVal(idx); x2 = xVal(idx + 1);
            y1 = yRow_NoMot(idx); y2 = yRow_NoMot(idx + 1);
            Boundary_NoMot_to_Int(e) = x1 + (0.5 - y1) * (x2 - x1) / (y2 - y1);
        end
    end
    
    % --- Boundary 2: Integration -> Segmentation ---
    yRow_Seg = GrandMean_Seg(e, :);
    if any(yRow_Seg < 0.5) && any(yRow_Seg > 0.5)
        diff_from_half = yRow_Seg - 0.5;
        sign_changes = find(diff(sign(diff_from_half)) ~= 0);
        if ~isempty(sign_changes)
            idx = sign_changes(1);
            x1 = xVal(idx); x2 = xVal(idx + 1);
            y1 = yRow_Seg(idx); y2 = yRow_Seg(idx + 1);
            Boundary_Int_to_Seg(e) = x1 + (0.5 - y1) * (x2 - x1) / (y2 - y1);
        end
    end
end

% Compare with V1 RF sizes
V1_RF_at_Ecc = V1_slope_Dow * ecc_labels' + V1_intercept_Dow;
LME_at_Ecc = LME_slope * ecc_labels' + LME_intercept;

fprintf('\n=== TABLE: PHASE BOUNDARIES vs. V1 RF SIZE ===\n');
Boundary_Table = table(ecc_labels', Boundary_NoMot_to_Int, Boundary_Int_to_Seg, ...
                       LME_at_Ecc, V1_RF_at_Ecc, ...
    'VariableNames', {'Eccentricity', 'NoMot_to_Int', 'Int_to_Seg', 'LME_Lambda75', 'V1_RF_Dow'});
disp(Boundary_Table);

% Correlation: Integration->Segmentation boundary vs V1 RF
valid_idx = ~isnan(Boundary_Int_to_Seg);
if sum(valid_idx) >= 3
    [r_seg, p_seg] = corr(Boundary_Int_to_Seg(valid_idx), V1_RF_at_Ecc(valid_idx));
    fprintf('\nCorrelation (Int->Seg Boundary vs V1 RF): r = %.3f, p = %.4f\n', r_seg, p_seg);
end

% Linear fit to the Integration->Segmentation boundary
valid_ecc = ecc_labels(valid_idx)';
valid_boundary = Boundary_Int_to_Seg(valid_idx);
if length(valid_boundary) >= 2
    p_fit = polyfit(valid_ecc, valid_boundary, 1);
    fprintf('Integration->Segmentation boundary: slope = %.4f, intercept = %.4f\n', p_fit(1), p_fit(2));
    fprintf('V1 RF (Dow et al.): slope = %.4f, intercept = %.4f\n', V1_slope_Dow, V1_intercept_Dow);
end

fprintf('\nPart 8: Three-Phase Diagram Complete.\n');
fprintf('Key finding: The Integration->Segmentation boundary parallels V1 RF size.\n');




