clear all
close all
clc

% Settings
% Path to the results
folderReal = 'C:\Users\Sophie Fleischmann\Documents\ResearchProjects\BeReal\ISB\BioMAC-Sim-Toolbox\results\TCSG\real';
folderIdeal = 'C:\Users\Sophie Fleischmann\Documents\ResearchProjects\BeReal\ISB\BioMAC-Sim-Toolbox\results\TCSG\ideal';

save_location = 'C:\Users\Sophie Fleischmann\Documents\ResearchProjects\BeReal\ISB\BioMAC-Sim-Toolbox\results\TCSG';

% Initialize cell array with column names
resVar = cell(10, 6); % Preallocate for efficiency

% Define column headers
columnNames = {'ExtractedData', 'Type', 'MetCost', 'MetCostPerMus', 'MetRate', 'CoT'};

% Process ideal data
for i = 1:5
    resultFile = [folderIdeal, filesep, '2025_02_26_script2D12_from_standing', num2str(i)];
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

    % Compute metabolic cost
    [metCost, dmetCostdX, metCostPerMus, metRate, CoT] = result.problem.getMetabolicCost(result.X);
    
    resVar{i, 3} = metCost;
    resVar{i, 4} = metCostPerMus;
    resVar{i, 5} = metRate;
    resVar{i, 6} = CoT;
end

% Process real data
for i = 1:5
    resultFile = [folderReal, filesep, '2025_02_26_script2D_BeReal12_from_standing', num2str(i)];
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

    % Compute metabolic cost
    [metCost, dmetCostdX, metCostPerMus, metRate, CoT] = result.problem.getMetabolicCost(result.X);
    
    resVar{i+5, 3} = metCost;
    resVar{i+5, 4} = metCostPerMus;
    resVar{i+5, 5} = metRate;
    resVar{i+5, 6} = CoT;
end

% Convert resVar to a table for better readability
resTable = cell2table(resVar, 'VariableNames', columnNames);

% Extract indices for ideal and real trials
%idealIdx = strcmp(resTable.Type, 'ideal');
%realIdx = strcmp(resTable.Type, 'real');

% Initialize results structure
%results = struct();

% Compute Mean & Median for Overall Metabolic Measures
% metrics = {'MetCost', 'MetRate', 'CoT'};
% for i = 1:length(metrics)
%     metric = metrics{i};
%     
%     % Convert table column to array
%     idealValues = table2array(resTable(idealIdx, metric)); 
%     realValues = table2array(resTable(realIdx, metric));
% 
%     % Compute mean & median
%     results.(metric).mean.ideal = mean(idealValues, 'omitnan');
%     results.(metric).median.ideal = median(idealValues, 'omitnan');
%     results.(metric).mean.real = mean(realValues, 'omitnan');
%     results.(metric).median.real = median(realValues, 'omitnan');
%     
%     % Display results
%     fprintf('\n%s:\n', metric);
%     fprintf('  Ideal -> Mean: %.4f, Median: %.4f\n', results.(metric).mean.ideal, results.(metric).median.ideal);
%     fprintf('  Real  -> Mean: %.4f, Median: %.4f\n', results.(metric).mean.real, results.(metric).median.real);
% end
% 
% % Get number of muscles (assuming 18 muscles per trial)
% numMuscles = size(resTable.MetCostPerMus{1}, 1);
% 
% % Create muscle index labels (assuming ordered muscle names)
% muscleNames = arrayfun(@(x) sprintf('Muscle_%d', x), 1:numMuscles, 'UniformOutput', false);
% 
% % Initialize result tables
% metCostPerMus_mean = table('Size', [numMuscles, 3], ...
%                            'VariableTypes', {'string', 'double', 'double'}, ...
%                            'VariableNames', {'Muscle', 'Mean_Ideal', 'Mean_Real'});
% 
% metCostPerMus_median = table('Size', [numMuscles, 3], ...
%                              'VariableTypes', {'string', 'double', 'double'}, ...
%                              'VariableNames', {'Muscle', 'Median_Ideal', 'Median_Real'});
% 
% % Extract MetCostPerMus data from resTable
% idealValues = cat(3, resTable.MetCostPerMus{idealIdx}); % Stack along 3rd dimension
% realValues = cat(3, resTable.MetCostPerMus{realIdx});   % Stack along 3rd dimension
% 
% % Compute mean & median per muscle
% meanIdeal = mean(idealValues, 3, 'omitnan');
% medianIdeal = median(idealValues, 3, 'omitnan');
% meanReal = mean(realValues, 3, 'omitnan');
% medianReal = median(realValues, 3, 'omitnan');
% 
% % Fill result tables
% for i = 1:numMuscles
%     metCostPerMus_mean.Muscle(i) = muscleNames{i};
%     metCostPerMus_mean.Mean_Ideal(i) = meanIdeal(i);
%     metCostPerMus_mean.Mean_Real(i) = meanReal(i);
% 
%     metCostPerMus_median.Muscle(i) = muscleNames{i};
%     metCostPerMus_median.Median_Ideal(i) = medianIdeal(i);
%     metCostPerMus_median.Median_Real(i) = medianReal(i);
% end
% 
% % Display tables
% disp('Mean MetCost Per Muscle:');
% disp(metCostPerMus_mean);






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
    
    % Loop through each variable under this data type
    for v = 1:numVars
        varName = varNames{v};
         if ~contains(varName, '_r')
            continue;
        end
        
        % Initialize arrays for ideal and real
        idealData = [];
        realData = [];
        
        % Extract relevant data for each trial
        for i = 1:10
            extractedData = resTable.ExtractedData{i};
       
            idx = strcmp(extractedData.name, varName) & strcmp(extractedData.type, dataType);

            simData = extractedData.sim(idx); % Cell array containing [Nx1 double]
            
            % Convert to matrix form
            dataMatrix = cell2mat(simData);
            
            if strcmp(resTable.Type{i}, 'ideal')
                idealData = cat(3, idealData, dataMatrix);
            else
                realData = cat(3, realData, dataMatrix);
            end
        end
        
        % Compute mean and std
        meanIdeal = mean(idealData, 3);
        stdIdeal = std(idealData, 0, 3);
        meanReal = mean(realData, 3);
        stdReal = std(realData, 0, 3);
        
        % Time vector (assuming all have same length)
        timeVector = [1:1:100];
        
        % Create subplot
        subplot(rows, cols, v);
        hold on;
        
        % Plot mean and shaded std for ideal
        fill([timeVector, fliplr(timeVector)], ...
             [meanIdeal + stdIdeal; flipud(meanIdeal - stdIdeal)], ...
             'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        plot(timeVector, meanIdeal, 'b');
        
        % Plot mean and shaded std for real
        fill([timeVector, fliplr(timeVector)], ...
             [meanReal + stdReal; flipud(meanReal - stdReal)], ...
             'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        plot(timeVector, meanReal, 'r');
        
        % Labels and legend
        title(varName, 'Interpreter', 'none');
        xlabel('Time (%)');
        ylabel('Value');
       % legend({'Ideal ± STD', 'Ideal Mean', 'Real ± STD', 'Real Mean'}, 'Location', 'best');
        
        hold off;
    end
    
    % Adjust figure layout
    sgtitle(['Mean ± STD for ', dataType], 'Interpreter', 'none');
end



