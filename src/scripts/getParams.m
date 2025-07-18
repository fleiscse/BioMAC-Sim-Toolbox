clear; clc; close all;

% Setup
dataFolder = 'results/TCSG_improve/best';
save_location = dataFolder;

speedLabels = {'1', '12', '14', '16', '18', '2'};
xLabels = {'1.0', '1.2', '1.4', '1.6', '1.8', '2.0'};
nSpeeds = length(speedLabels);

% Predefine settings
settings.angle  = {'hip_flexion_r', 'knee_angle_r', 'ankle_angle_r', 'pelvis_ty', 'pelvis_tilt'};
settings.moment = {'hip_flexion_r', 'knee_angle_r', 'ankle_angle_r', 'pelvis_ty', 'pelvis_tilt'};
settings.u      = {'hamstrings_r', 'bifemsh_r', 'glut_max_r', 'iliopsoas_r', 'rect_fem_r', 'vasti_r', 'gastroc_r', 'soleus_r', 'tib_ant_r'};
settings.a      = settings.u;

% Create empty containers for extracted data
extractedData.ideal = cell(1, nSpeeds);
extractedData.realistic = cell(1, nSpeeds);

% Load and extract for each speed and condition
for i = 1:nSpeeds
    for cond = ["ideal", "realistic"]
        fileName = fullfile(dataFolder, sprintf('%s%s.mat', cond, speedLabels{i}));
        if isfile(fileName)
            res = load(fileName);
            result = res.result;
            extracted = result.problem.extractData(result.X, settings, [], 1);
            extractedData.(cond){i} = extracted;  % ✅ now safe
        else
            warning('File not found: %s', fileName);
        end
    end
end


% Get all data types
dataTypes = unique(extractedData.ideal{1}.type);

%% Plot loop for each data type
for d = 1:length(dataTypes)
    dataType = dataTypes{d};
    
    % Find all variables of this type (assume same across speeds)
    allVars = extractedData.ideal{1}.name(strcmp(extractedData.ideal{1}.type, dataType));
    
    % Create figure
    figure;
    numVars = length(allVars);
    rows = ceil(sqrt(numVars));
    cols = ceil(numVars / rows);
    
    for v = 1:numVars
        varName = allVars{v};
        if ~contains(varName, '_r'), continue; end

        % Initialize for plotting
        idealTraces = NaN(nSpeeds, 100);
        realTraces  = NaN(nSpeeds, 100);

        for i = 1:nSpeeds
            idxIdeal = strcmp(extractedData.ideal{i}.name, varName) & strcmp(extractedData.ideal{i}.type, dataType);
            idxReal  = strcmp(extractedData.realistic{i}.name, varName) & strcmp(extractedData.realistic{i}.type, dataType);

            if any(idxIdeal)
                idealTraces(i, :) = extractedData.ideal{i}.sim{idxIdeal}(:)';
            end
            if any(idxReal)
                realTraces(i, :) = extractedData.realistic{i}.sim{idxReal}(:)';
            end
        end

        % Plot
        subplot(rows, cols, v);
        hold on;
        plot(idealTraces', '--', 'Color', [0.2 0.4 1], 'LineWidth', 1);  % Ideal
        plot(realTraces', '-',  'Color', [1 0.2 0.2], 'LineWidth', 1);  % Realistic
        title(varName, 'Interpreter', 'none');
        xlabel('% Stride'); ylabel('Value');
        xticks(linspace(1, 100, 5));
        box off; grid on;
        legend('Ideal', 'Realistic');
    end

    sgtitle(['Across Speeds – ', dataType], 'Interpreter', 'none');
    saveas(gcf, fullfile(save_location, ['AcrossSpeeds_', dataType, '.png']));
end
