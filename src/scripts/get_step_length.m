clear; clc; close all;

%% Settings
dataFolder = 'results/TCSG_improve/best';
save_location = dataFolder;

% Speed labels used in filenames
speedLabels = {'1', '12', '14', '16', '18', '2'};

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
        SL = result.problem.getStepLength(result.X, 'left');

        % Store results
        resVar{row, 1} = str2double(strrep(speedLabel, 'p', '.'));  % convert '12' to 1.2 if needed
        resVar{row, 2} = char(condition);
        resVar{row, 3} = SL;
        
        row = row + 1;
    end
end

%% Convert to Table
resTable = cell2table(resVar, 'VariableNames', columnNames);

