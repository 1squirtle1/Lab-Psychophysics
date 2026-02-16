% Load the data from an Excel file
filename = 'pairedunpaired_psychophysicsdatav2.xlsx';
sheets = {'RB_analysis', 'AP_analysis', 'BG_analysis'};

% Define Weibull function
weibull_func = @(p, x) 1 - exp(-((x./p(1)).^p(2)));

% Initial guess for parameters [lambda, k]
p0 = [0.6, 3];

% Define a color map for different degrees
colors = {'b', 'r', 'g', 'm', 'c'};

thresholds = zeros(5, length(sheets));


% Loop through each sheet
for sheet_idx = 1:length(sheets)
    sheet = sheets{sheet_idx};
    data = xlsread(filename, sheet);
%     x_values = data(1:10, 1:10);
    x_values = data(1:10, 1);% Assuming the first column is the x-values (stimulus pairing range)
    figure;
    % Loop through each eccentricity (3, 4, 5, 6, 7 deg)
    j = 1;
    for i = 2:2:11
        y_values = (data(1:10, i)*100) / 100; % Normalize the y-values
        
        % Fit the Weibull function
        [p_fit, ~] = lsqcurvefit(weibull_func, p0, x_values, y_values, [0, 0], [2, 50]);
        
        % Generate the fitted curve
        x_fit = linspace(min(x_values), max(x_values), 100);
        y_fit = weibull_func(p_fit, x_fit);

        % Plot with unique color
        plot(x_values, (data(1:10, i)*100), 'o', 'MarkerFaceColor', colors{j}, 'Color', colors{j}); % Original data points
        hold on;
        plot(x_fit, y_fit * 100, '--', 'Color', colors{j}); % Fitted curve (scaled back to percentages)
        
        % Display the fitted parameters
        lambda = p_fit(1);
        k = p_fit(2);
        fprintf('%s - %d deg:\n', sheet, i+2);
        fprintf('Lambda: %f\n', lambda);
        fprintf('k: %f\n', k);

        % Calculate the 50% threshold
%         threshold_50 = lambda * (log(2))^(1/k);
%         fprintf('50%% Threshold: %f\n\n', threshold_50);
%         thresholds(i-1, sheet_idx) = lambda * (log(2))^(1/k);
        %calculate the 75% threshold
        thresholds(j, sheet_idx) = lambda * (log(1/(1-0.75)))^(1/k);
        j=j+1;

    end
    
    disp('50% Thresholds for each sheet:');
    disp(array2table(thresholds, 'RowNames', {'3 deg', '4 deg', '5 deg', '6 deg', '7 deg'}, 'VariableNames', sheets));


    xlabel('Stimuli Pairing Range (deg)');
    ylabel('Percentage of Time Seen as Segmented (%)');
    title([sheet]);
    legend('3 deg','3 fit', '4 deg','4 fit', '5 deg','5 fit', '6 deg','6fit', '7 deg','7 fit',Location='best');
    grid on;
end
%convert threshold val to maximal spatial seperation
thresholds_max_spatial = thresholds.*sin(pi/4);

% Assuming 'thresholds' contains the 50% threshold values
% where each row corresponds to an eccentricity (3, 4, 5, 6, 7 deg)
avg_thresholds = mean(thresholds(:,[1 3:4]), 2); % Compute the average threshold for each eccentricity

% Define x values corresponding to the eccentricities
x = [3, 4, 5, 6, 7]';

% Fit a linear regression model
p = polyfit(x, avg_thresholds, 1);

% Extract the slope and intercept
slope = p(1);
intercept = p(2);

% Display the results
fprintf('Slope: %f\n', slope);
fprintf('Intercept: %f\n', intercept);

