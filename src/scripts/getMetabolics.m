clear; clc; close all;

%% Settings
dataFolder = 'results/best2';
save_location = dataFolder;

% Speed labels used in filenames
speedLabels = { '1','12', '14', '16', '18'};

% Metrics to extract
metrics = {'MetCost', 'MetRate', 'CoT'};

% Initialize result storage
resVar = cell(length(speedLabels) * 2, 6);
columnNames = {'Speed', 'Type', 'MetCost', 'MetCostPerMus', 'MetRate', 'CoT'};

row = 1;
for i = 1:length(speedLabels)
    speedLabel = speedLabels{i};

    for condition = ["ideal", "realistic"]
        fileName = fullfile(dataFolder, sprintf('%s%s.mat', condition, speedLabel));
        if ~isfile(fileName)
            warning('File not found: %s', fileName);
            continue;
        end

        res = load(fileName);
        result = res.result;

        % Get metabolic cost metrics
        [metCost, ~, metCostPerMus, metRate, CoT] = result.problem.getMetabolicCost(result.X);

        % Store results
        resVar{row, 1} = str2double(strrep(speedLabel, 'p', '.'));  % convert '12' to 1.2 if needed
        resVar{row, 2} = char(condition);
        resVar{row, 3} = metCost;
        resVar{row, 4} = metCostPerMus;
        resVar{row, 5} = metRate;
        resVar{row, 6} = CoT;

        row = row + 1;
    end
end

%% Convert to Table
resTable = cell2table(resVar, 'VariableNames', columnNames);

%% Summary: Mean values per speed and condition
results = struct();
for i = 1:length(metrics)
    metric = metrics{i};

    fprintf('\n%s:\n', metric);
    fprintf('%6s\t%10s\t%10s\n', 'Speed', 'Ideal', 'Realistic');

    for s = 1:length(speedLabels)
        speed = str2double(strrep(speedLabels{s}, 'p', '.'));
        idealRow = strcmp(resTable.Type, 'ideal') & resTable.Speed == speed;
        realRow  = strcmp(resTable.Type, 'realistic') & resTable.Speed == speed;

        if any(idealRow) && any(realRow)
            idealVal = resTable{idealRow, metric};
            realVal  = resTable{realRow, metric};

            results.(metric).speed(s) = speed;
            results.(metric).ideal(s) = idealVal;
            results.(metric).realistic(s) = realVal;

            fprintf('%6.2f\t%10.4f\t%10.4f\n', speed, idealVal, realVal);
        end
    end
end

%% Per-Muscle Mean Across All Speeds
numMuscles = size(resTable.MetCostPerMus{find(~cellfun(@isempty, resTable.MetCostPerMus), 1)}, 1);
muscleNames = arrayfun(@(x) sprintf('Muscle_%d', x), 1:numMuscles, 'UniformOutput', false);

% Preallocate for means
meanIdealPerMus = zeros(numMuscles, 1);
meanRealPerMus = zeros(numMuscles, 1);

% Stack muscle cost arrays
idealStacks = cat(3, resTable.MetCostPerMus{strcmp(resTable.Type, 'ideal')});
realStacks  = cat(3, resTable.MetCostPerMus{strcmp(resTable.Type, 'realistic')});

% Compute means
meanIdealPerMus = mean(idealStacks, 3);
meanRealPerMus = mean(realStacks, 3);

% Create summary table
muscleSummary = table(muscleNames', meanIdealPerMus, meanRealPerMus, ...
                      'VariableNames', {'Muscle', 'Mean_Ideal', 'Mean_Realistic'});

disp('Mean Metabolic Cost per Muscle across all speeds:');
disp(muscleSummary);

%% Plot % of Realistic MetRate relative to Ideal
% Your custom label order
customOrder = [1, 12, 14, 16, 18];
xLabels = { '1','1.2', '1.4', '1.6', '1.8'};

% Prepare y data in this same order
percentRealOverIdeal = zeros(length(customOrder), 1);
for i = 1:length(customOrder)
    spd = customOrder(i);
    idealRate = resTable.MetRate(strcmp(resTable.Type, 'ideal') & resTable.Speed == spd);
    realRate  = resTable.MetRate(strcmp(resTable.Type, 'realistic') & resTable.Speed == spd);

    if isempty(idealRate) || isempty(realRate)
        percentRealOverIdeal(i) = NaN;
    else
        percentRealOverIdeal(i) = 100 * realRate / idealRate;
    end
end

% Plot with categorical x-axis
figure; hold on; box off;
x = 1:length(customOrder);
plot(x, percentRealOverIdeal, '-o', 'LineWidth', 2, 'Color', [0.2 0.4 0.7]);
yline(100, '--k', 'LineWidth', 1.5);

xticks(x);
xticklabels(xLabels);
xlabel('Speed (m/s)');
ylabel('Realistic MetRate (% of Ideal)');
title('Relative Metabolic Rate: Realistic vs Ideal');
ylim([min(percentRealOverIdeal)-5, max(percentRealOverIdeal)+5]);
grid on;
set(gca, 'FontSize', 12);


% Setup
dataFolder = 'results/TCSG_improve/best/';
speedLabels = { '1','12', '14', '16', '18', };
xLabels = { '1','1.2', '1.4', '1.6', '1.8', };
nSpeeds = length(speedLabels);

% Joints of interest
joints = {'hip_flexion_r', 'knee_angle_r', 'ankle_angle_r'};
nJoints = length(joints);

% Initialize storage: [nJoints x nSpeeds]
maxMoment.ideal = NaN(nJoints, nSpeeds);
minMoment.ideal = NaN(nJoints, nSpeeds);
maxMoment.realistic = NaN(nJoints, nSpeeds);
minMoment.realistic = NaN(nJoints, nSpeeds);

% Settings for extraction
settings.moment = joints;

% Extract data
for i = 1:nSpeeds
    for cond = ["ideal", "realistic"]
        fileName = fullfile(dataFolder, sprintf('%s%s.mat', cond, speedLabels{i}));
        if ~isfile(fileName)
            warning('File not found: %s', fileName);
            continue;
        end

        % Load and extract
        res = load(fileName);
        result = res.result;
        extracted = result.problem.extractData(result.X, settings, [], 1);

        % Loop through joints and extract min/max
        for j = 1:nJoints
            idx = strcmp(extracted.name, joints{j}) & strcmp(extracted.type, 'moment');
            if any(idx)
                data = extracted.sim{idx};
                minVal = min(data);
                maxVal = max(data);

                if cond == "ideal"
                    minMoment.ideal(j, i) = minVal;
                    maxMoment.ideal(j, i) = maxVal;
                else
                    minMoment.realistic(j, i) = minVal;
                    maxMoment.realistic(j, i) = maxVal;
                end
            end
        end
    end
end



% Define speed order (corresponding to moment matrix columns)
customOrder = [1, 12, 14, 16, 18, 2]; % Match moment matrix ordering
nSpeeds = length(customOrder);
nJoints = length(joints);

% Compute metabolic rate ratio per speed
metRateRatio = NaN(nSpeeds, 1);
for i = 1:nSpeeds
    spd = customOrder(i);
    realIdx = strcmp(resTable.Type, 'realistic') & resTable.Speed == spd;
    idealIdx = strcmp(resTable.Type, 'ideal') & resTable.Speed == spd;

    if any(realIdx) && any(idealIdx)
        realMet = resTable.MetRate(realIdx);
        idealMet = resTable.MetRate(idealIdx);
        metRateRatio(i) = realMet / idealMet;
    end
end

% Compute difference in max and min moments per joint per speed
deltaMax = (maxMoment.realistic ./ maxMoment.ideal);  % [joint x speed]
deltaMin = (minMoment.realistic ./ minMoment.ideal);

% Now correlate for each joint
fprintf('\nCorrelation between Δ moment and metabolic rate ratio:\n');
for j = 1:nJoints
    % Extract differences across speeds for this joint
    dx_max = deltaMax(j, :)';
    dx_min = deltaMin(j, :)';

    % Pearson correlation (ignores NaNs automatically)
    [r_max, p_max] = corr(dx_max, metRateRatio, 'rows', 'complete');
    [r_min, p_min] = corr(dx_min, metRateRatio, 'rows', 'complete');

    fprintf('Joint: %-15s | Max: r = %.3f (p = %.3f), Min: r = %.3f (p = %.3f)\n', ...
        joints{j}, r_max, p_max, r_min, p_min);
end


jointIdx = find(strcmp(joints, 'ankle_angle_r'));

figure;
subplot(1,2,1)
scatter(deltaMax(jointIdx, :)', metRateRatio, 80, 'filled');
xlabel('Δ Max Moment (Nm/kg)');
ylabel('MetRate Ratio (Real/Ideal)');
title('Ankle Max Moment vs MetRate');

subplot(1,2,2)
scatter(deltaMin(jointIdx, :)', metRateRatio, 80, 'filled');
xlabel('Δ Min Moment (Nm/kg)');
ylabel('MetRate Ratio (Real/Ideal)');
title('Ankle Min Moment vs MetRate');
