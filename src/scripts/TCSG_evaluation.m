clear all
close all
clc

% Settings
% Path to the results
folderOverground = 'C:\Users\Sophie Fleischmann\Documents\ResearchProjects\BeReal\ISB\BioMAC-Sim-Toolbox\results\TCSG\overground';

folderReal = 'C:\Users\Sophie Fleischmann\Documents\ResearchProjects\BeReal\ISB\BioMAC-Sim-Toolbox\results\TCSG\real';
folderIdeal = 'C:\Users\Sophie Fleischmann\Documents\ResearchProjects\BeReal\ISB\BioMAC-Sim-Toolbox\results\TCSG\ideal';

folderIdeal = 'C:\Users\Sophie Fleischmann\Documents\ResearchProjects\BeReal\ISB\BioMAC-Sim-Toolbox\results\TCSG\18\ideal_from_overground6';
folderReal= 'C:\Users\Sophie Fleischmann\Documents\ResearchProjects\BeReal\ISB\BioMAC-Sim-Toolbox\results\TCSG\18\real_from_overground6';


save_location = 'C:\Users\Sophie Fleischmann\Documents\ResearchProjects\BeReal\ISB\BioMAC-Sim-Toolbox\results\TCSG';

% Initialize cell array with column names
resVar = cell(15, 2); % Preallocate for efficiency

% Define column headers
columnNames = {'ExtractedData', 'Type'};

% Process ideal data
for i = 1:5
     filePattern = fullfile(folderIdeal, ['*', num2str(i), '.mat']);
    fileList = dir(filePattern);
    resultFile = fullfile(folderIdeal, fileList(1).name);
    res = load(resultFile);
    result = res.result;
    
    settings.angle  = {'hip_flexion_r', 'knee_angle_r', 'ankle_angle_r', 'pelvis_ty', 'pelvis_tilt'};
    settings.moment = {'hip_flexion_r', 'knee_angle_r', 'ankle_angle_r', 'pelvis_ty', 'pelvis_tilt'};
    settings.u      = {'hamstrings_r', 'bifemsh_r', 'glut_max_r', 'iliopsoas_r', 'rect_fem_r', 'vasti_r', 'gastroc_r', 'soleus_r', 'tib_ant_r'};
    settings.a      = {'hamstrings_r', 'bifemsh_r', 'glut_max_r', 'iliopsoas_r', 'rect_fem_r', 'vasti_r', 'gastroc_r', 'soleus_r', 'tib_ant_r'};

    step_dur = result.X(result.problem.idx.dur); % Duration of movement
    
    % Extract data and store in resVar
    resVar{i, 1} = result.problem.extractData(result.X, settings, [], 1);
    resVar{i, 2} = 'ideal';

end
% Process real data
for i = 1:5
    filePattern = fullfile(folderReal, ['*', num2str(i), '.mat']);
    fileList = dir(filePattern);
    resultFile = fullfile(folderReal, fileList(1).name);
    res = load(resultFile);
    result = res.result;
    
    settings.angle  = {'hip_flexion_r', 'knee_angle_r', 'ankle_angle_r', 'pelvis_ty', 'pelvis_tilt'};
    settings.moment = {'hip_flexion_r', 'knee_angle_r', 'ankle_angle_r', 'pelvis_ty', 'pelvis_tilt'};
    settings.u      = {'hamstrings_r', 'bifemsh_r', 'glut_max_r', 'iliopsoas_r', 'rect_fem_r', 'vasti_r', 'gastroc_r', 'soleus_r', 'tib_ant_r'};
    settings.a      = {'hamstrings_r', 'bifemsh_r', 'glut_max_r', 'iliopsoas_r', 'rect_fem_r', 'vasti_r', 'gastroc_r', 'soleus_r', 'tib_ant_r'};

    step_dur = result.X(result.problem.idx.dur); % Duration of movement
    
    % Extract data and store in resVar
    resVar{i+5, 1} = result.problem.extractData(result.X, settings, [], 1);
    resVar{i+5, 2} = 'real';

  
end
for i = 1:5
    filePattern = fullfile(folderOverground, ['*', num2str(i), '.mat']);
    fileList = dir(filePattern);
    resultFile = fullfile(folderOverground, fileList(1).name);
    res = load(resultFile);
    result = res.result;
    
    settings.angle  = {'hip_flexion_r', 'knee_angle_r', 'ankle_angle_r', 'pelvis_ty', 'pelvis_tilt'};
    settings.moment = {'hip_flexion_r', 'knee_angle_r', 'ankle_angle_r', 'pelvis_ty', 'pelvis_tilt'};
    settings.u      = {'hamstrings_r', 'bifemsh_r', 'glut_max_r', 'iliopsoas_r', 'rect_fem_r', 'vasti_r', 'gastroc_r', 'soleus_r', 'tib_ant_r'};
    settings.a      = {'hamstrings_r', 'bifemsh_r', 'glut_max_r', 'iliopsoas_r', 'rect_fem_r', 'vasti_r', 'gastroc_r', 'soleus_r', 'tib_ant_r'};

    step_dur = result.X(result.problem.idx.dur); % Duration of movement
    
    % Extract data and store in resVar
    resVar{i+10, 1} = result.problem.extractData(result.X, settings, [], 1);
    resVar{i+10, 2} = 'overground';


  
end
% Convert resVar to a table for better readability
resTable = cell2table(resVar, 'VariableNames', columnNames);

meanDataStruct = struct();

% Extract unique data types
dataTypes = unique(resTable.ExtractedData{1}.type);

% Loop through each data type to create separate figures
for d = 1:length(dataTypes)
    dataType = dataTypes{d};
    
    % Extract all variable names under this data type
    extractedData = resTable.ExtractedData{1}; % Use first entry to get all variables
    varNames = extractedData.name(strcmp(extractedData.type, dataType));
    
    % Create a new figure for this data type
    figure;
    numVars = length(varNames);
    rows = ceil(sqrt(numVars)); % Adjust for grid layout
    cols = ceil(numVars / rows);
    

    meanDataStruct.(dataType) = struct();
    % Loop through each variable under this data type
    for v = 1:numVars
        varName = varNames{v};
         if ~contains(varName, '_r')
            continue;
        end
        
        % Initialize arrays for ideal and real
        idealData = [];
        realData = [];
        overgroundData = [];
        
        % Extract relevant data for each trial#
      
        for i = 1:15
            extractedData = resTable.ExtractedData{i};
       
            idx = strcmp(extractedData.name, varName) & strcmp(extractedData.type, dataType);

            simData = extractedData.sim(idx); % Cell array containing [Nx1 double]
            
            % Convert to matrix form
            dataMatrix = cell2mat(simData);
            
            if strcmp(resTable.Type{i}, 'ideal')
                idealData = cat(3, idealData, dataMatrix);
            elseif strcmp(resTable.Type{i}, 'real')
                realData = cat(3, realData, dataMatrix);
            else
                overgroundData = cat(3, overgroundData, dataMatrix);
            end
        end
        
        % Compute mean and std
        meanIdeal = mean(idealData, 3);
        stdIdeal = std(idealData, 0, 3);
        meanReal = mean(realData, 3);
        stdReal = std(realData, 0, 3);
        meanOverground = mean(overgroundData, 3);

        meanDataStruct.(dataType).(varName) = struct();
        meanDataStruct.(dataType).(varName).ideal = meanIdeal;
        meanDataStruct.(dataType).(varName).real = meanReal;
        meanDataStruct.(dataType).(varName).overground = meanOverground;
        
        % Time vector (assuming all have same length)
        timeVector = [1:1:100];
        
        % Create subplot
        subplot(rows, cols, v);
        hold on;
        
        % Plot mean and shaded std for ideal
%         fill([timeVector, fliplr(timeVector)], ...
%              [meanIdeal + stdIdeal; flipud(meanIdeal - stdIdeal)], ...
%              'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        plot(timeVector, meanIdeal, 'b');
        
%         % Plot mean and shaded std for real
%         fill([timeVector, fliplr(timeVector)], ...
%              [meanReal + stdReal; flipud(meanReal - stdReal)], ...
%              'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        plot(timeVector, meanReal, 'r');
      %  plot(timeVector, meanOverground, 'k');
        
        % Labels and legend
        title(varName, 'Interpreter', 'none');
        xlabel('Time (%)');
        ylabel('Value');
       % legend({'Ideal ± STD', 'Ideal Mean', 'Real ± STD', 'Real Mean'}, 'Location', 'best');
      
        hold off;

        





    end
    
    % Adjust figure layout
    sgtitle(['Mean for ', dataType], 'Interpreter', 'none');
    saveas(gcf, fullfile(save_location, ['18Mean_', dataType, '.png']));
end




% % Initialize arrays for stacked format (Method 1)
% Y = [];
% A = [];
% 
% % Initialize cell arrays for separated format (Method 2)
% idealTrials = {};
% realTrials = {};
% 
% % Loop through each trial and extract individual data
% for i = 1:15
%     extractedData = resTable.ExtractedData{i};
%     
%     % Extract all variable names under the first entry
%     if i == 1
%         dataTypes = unique(extractedData.type);
%     end
%     
%     for d = 1:length(dataTypes)
%         dataType = dataTypes{d};
%         
%         % Get relevant variable names
%         varNames = extractedData.name(strcmp(extractedData.type, dataType));
%         
%         for v = 1:length(varNames)
%             varName = varNames{v};
%             if ~contains(varName, '_r')
%                 continue;
%             end
%             
%             % Find the data for this variable
%             idx = strcmp(extractedData.name, varName) & strcmp(extractedData.type, dataType);
%             simData = extractedData.sim(idx);
%             dataMatrix = cell2mat(simData); % Convert cell to matrix
%             
%             % Append to stacked format
%             Y = [Y; dataMatrix(:)']; % Stack rows
%             if strcmp(resTable.Type{i}, 'ideal')
%                 A = [A; ones(size(dataMatrix, 1), 1)];  % Label "1" for ideal
%                 idealTrials{end+1} = dataMatrix;  % Store separately
%             elseif strcmp(resTable.Type{i}, 'real')
%                 A = [A; 2 * ones(size(dataMatrix, 1), 1)]; % Label "2" for real
%                 realTrials{end+1} = dataMatrix;  % Store separately
%             end
%         end
%     end
% end
% 
% % Convert separated trials to MATLAB struct
% separatedStruct.ideal = idealTrials;
% separatedStruct.real = realTrials;
% 
% % Save the variables for Python processing
% save(fullfile(save_location, 'spm1d_data.mat'), 'Y', 'A', 'separatedStruct');
% 
% 
% 
% 

% Compute the muscle activation differences (ideal - real)
muscles = {'hamstrings_r', 'bifemsh_r', 'glut_max_r', 'iliopsoas_r', 'rect_fem_r', ...
           'vasti_r', 'gastroc_r', 'soleus_r', 'tib_ant_r'};
%muscles = {'gastroc_r', 'soleus_r', 'tib_ant_r'};
%muscles = {'hamstrings_r', 'bifemsh_r',  'gastroc_r', 'rect_fem_r', ...
%           'vasti_r'};
activationDiffs = struct();
peaks = struct();
% Time vector (assuming all activations have the same length)
timeVector = 1:100;

for m = 1:length(muscles)
    muscle = muscles{m};
    
    % Extract mean activations for ideal and real conditions
    meanIdeal = meanDataStruct.a.(muscle).ideal;
    meanReal = meanDataStruct.a.(muscle).real;

    maxIdeal = max(meanIdeal)
    maxReal = max(meanReal)
    diff = maxReal / maxIdeal; %(maxReal - maxIdeal) / maxIdeal
    
    % Compute activation difference
    activationDiffs.(muscle) = meanIdeal - meanReal;
    peaks.(muscle) = diff;
end

% Plot all muscle activation differences on one plot
figure;
hold on;

colors = lines(length(muscles)); % Generate distinct colors

for m = 1:length(muscles)
    muscle = muscles{m};
    plot(timeVector, activationDiffs.(muscle), 'Color', colors(m, :), 'LineWidth', 2);
end

xlabel('Time (%)');
ylabel('Activation Difference (Ideal - Real)');
title('Muscle Activation Differences');
legend(muscles, 'Location', 'best');
grid on;
hold off;

% Save the plot
%saveas(gcf, fullfile(save_location, 'MuscleActivationDifferences.png'));


% Compute the muscle activation differences (ideal - real)
muscles_ankle = {'gastroc_r', 'soleus_r', 'tib_ant_r'};
muscles_knee = {'hamstrings_r', 'bifemsh_r', 'rect_fem_r', 'vasti_r'};
       
activationDiffs = struct();

% Load belt speed data


speeds_sim = [];

for i = 1:5
    filePattern = fullfile(folderReal, ['*', num2str(i), '.mat']);
    fileList = dir(filePattern);
    resultFile = fullfile(folderReal, fileList(1).name);
    res = load(resultFile);
    result = res.result;
    bR = result.problem.idx.belt_right;
    compSpeedR = result.X(bR);
    speeds_sim = cat(3, speeds_sim, compSpeedR);
end
meanSpeed_sim = mean(speeds_sim, 3);

% Define stance phase as ending at last index where simulated speed > 1.8
stance_end_idx = find(meanSpeed_sim > 1.805, 1, 'last');
stance_phase = 1:stance_end_idx;
meanSpeed_sim = meanSpeed_sim-1.8;

% Time vector (assuming all activations have the same length)
timeVector = 1:100;

% Compute activation differences for stance phase
activation_ankle = zeros(length(muscles_ankle), length(stance_phase));
activation_knee = zeros(length(muscles_knee), length(stance_phase));

for m = 1:length(muscles_ankle)
    muscle = muscles_ankle{m};
    meanIdeal = meanDataStruct.a.(muscle).ideal;
    meanReal = meanDataStruct.a.(muscle).real;
    activationDiff = meanIdeal - meanReal;
    activation_ankle(m, :) = activationDiff(stance_phase);
end

for m = 1:length(muscles_knee)
    muscle = muscles_knee{m};
    meanIdeal = meanDataStruct.a.(muscle).ideal;
    meanReal = meanDataStruct.a.(muscle).real;
    activationDiff = meanIdeal - meanReal;
    activation_knee(m, :) = activationDiff(stance_phase);
end

% Create subplots
fig = figure;
fig.Position(4) = 200;

% Subplot for ankle flexors/extensors
%subplot(2,2,1);

plot(stance_phase, meanSpeed_sim(stance_phase), 'Color', "#0072BD", 'LineWidth', 1);
%for m = 1:length(muscles_ankle)
hold on 
plot(stance_phase, activation_ankle(1, :), 'Color', "#0072BD", 'LineWidth', 1);
hold on 
plot(stance_phase, activation_ankle(2, :), 'Color', "#4DBEEE", 'LineWidth', 1);

plot(stance_phase, activation_ankle(3, :), 'Color', "#EDB120", 'LineWidth', 1);


xlabel('% of stride');
ylabel('activation difference');
%title('Ankle Flexors/Extensors');
legend([{'gastrocnemius', 'soleus' 'tibialis anterior'}], 'Location', 'southeast');
set(gca, 'box', 'off');
fontsize(gcf,scale=1.2);
%saveas(fig, 'muscleActDiff.png')
%grid on;

% % Subplot for knee flexors/extensors
% subplot(1,2,2);
% hold on;
% plot(stance_phase, meanSpeed_sim(stance_phase), 'Color', "#D95319", 'LineWidth', 1);
% for m = 1:length(muscles_knee)
%     plot(stance_phase, activation_knee(m, :), 'LineWidth', 1.5);
% end
% xlabel('% of Stance Phase');
% ylabel('Activation Difference & Belt Speed');
% title('Knee Flexors/Extensors');
% legend(['Measured Belt Speed', 'Simulated Belt Speed', muscles_knee], 'Location', 'southeast');
% set(gca, 'box', 'off');
% grid on;

% Save the plot

% Compute the mean muscle activation for all muscles (real and ideal)
% Compute the mean muscle activation for all muscles (real and ideal)
muscles = {'hamstrings_r', 'bifemsh_r', 'glut_max_r', 'iliopsoas_r', 'rect_fem_r', ...
           'vasti_r', 'gastroc_r', 'soleus_r', 'tib_ant_r'};
       
mean_activations_real = struct();
mean_activations_ideal = struct();
activation_differences = struct();

for m = 1:length(muscles)
    muscle = muscles{m};
    mean_ideal = mean(meanDataStruct.a.(muscle).ideal);
    mean_real = mean(meanDataStruct.a.(muscle).real);
    
    mean_activations_ideal.(muscle) = mean_ideal;
    mean_activations_real.(muscle) = mean_real;
    activation_differences.(muscle) = mean_ideal - mean_real;
end

% Print mean muscle activations for real and ideal conditions
disp('Mean Muscle Activations for Ideal Condition:');
disp(mean_activations_ideal);
disp('Mean Muscle Activations for Real Condition:');
disp(mean_activations_real);
disp('Difference in Mean Muscle Activations (Ideal - Real):');
disp(activation_differences);