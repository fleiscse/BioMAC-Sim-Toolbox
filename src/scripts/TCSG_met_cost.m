clear all
close all
clc

% Settings
% Path to the results
folderOverground = 'C:\Users\Sophie Fleischmann\Documents\ResearchProjects\BeReal\ISB\BioMAC-Sim-Toolbox\results\TCSG\18\overground6';
%folderReal = 'C:\Users\Sophie Fleischmann\Documents\ResearchProjects\BeReal\ISB\BioMAC-Sim-Toolbox\results\TCSG\real_from_overground';

folderReal = 'C:\Users\Sophie Fleischmann\Documents\ResearchProjects\BeReal\ISB\BioMAC-Sim-Toolbox\results\TCSG\18\real_from_overground6';
%folderIdeal = 'C:\Users\Sophie Fleischmann\Documents\ResearchProjects\BeReal\ISB\BioMAC-Sim-Toolbox\results\TCSG\ideal_from_overground';

folderIdeal = 'C:\Users\Sophie Fleischmann\Documents\ResearchProjects\BeReal\ISB\BioMAC-Sim-Toolbox\results\TCSG\18\ideal_from_overground6';

save_location = 'C:\Users\Sophie Fleischmann\Documents\ResearchProjects\BeReal\ISB\BioMAC-Sim-Toolbox\results\TCSG';

% Initialize cell array with column names
resVar = cell(15, 6); % Preallocate for efficiency

% Define column headers
columnNames = {'ExtractedData', 'Type', 'MetCost', 'MetCostPerMus', 'MetRate', 'CoT'};
total_dur = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOAD RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

    step_dur = result.X(result.problem.idx.dur) % Duration of movement
    total_dur = total_dur + step_dur;
    % Extract data and store in resVar
    resVar{i, 1} = 0;%result.problem.extractData(result.X, settings, [], 1);
    resVar{i, 2} = 'ideal';

    % Compute metabolic cost
    [metCost, dmetCostdX, metCostPerMus, metRate, CoT] = result.problem.getMetabolicCost(result.X);
    
    resVar{i, 3} = metCost;
    resVar{i, 4} = metCostPerMus;
    resVar{i, 5} = metRate;
    resVar{i, 6} = CoT;
    
end
step_dur = total_dur / 5
total_dur = 0;
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

    step_dur = result.X(result.problem.idx.dur) % Duration of movement
    total_dur = total_dur + step_dur;

    % Extract data and store in resVar
    resVar{i+5, 1} = 0;
    resVar{i+5, 2} = 'real';

    % Compute metabolic cost
    [metCost, dmetCostdX, metCostPerMus, metRate, CoT] = result.problem.getMetabolicCost(result.X);
    
    resVar{i+5, 3} = metCost;
    resVar{i+5, 4} = metCostPerMus;
    resVar{i+5, 5} = metRate;
    resVar{i+5, 6} = CoT;
end
step_dur = total_dur/5
% Process real data
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
    resVar{i+10, 1} = 0;
    resVar{i+10, 2} = 'overground';

    % Compute metabolic cost
    [metCost, dmetCostdX, metCostPerMus, metRate, CoT] = result.problem.getMetabolicCost(result.X);
    
    resVar{i+10, 3} = metCost;
    resVar{i+10, 4} = metCostPerMus;
    resVar{i+10, 5} = metRate;
    resVar{i+10, 6} = CoT;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GET MET COST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

resTable = cell2table(resVar, 'VariableNames', columnNames);

% Extract indices for ideal and real trials
idealIdx = strcmp(resTable.Type, 'ideal');
realIdx = strcmp(resTable.Type, 'real');
overgroundIdx = strcmp(resTable.Type, 'overground');

%Initialize results structure
results = struct();

%Compute Mean & Median for Overall Metabolic Measures
metrics = {'MetCost', 'MetRate', 'CoT'};
for i = 1:length(metrics)
    metric = metrics{i};
    
    % Convert table column to array
    idealValues = table2array(resTable(idealIdx, metric)); 
    realValues = table2array(resTable(realIdx, metric));
    overgroundValues = table2array(resTable(overgroundIdx, metric));

    % Compute mean & median
    results.(metric).mean.ideal = mean(idealValues);
    results.(metric).median.ideal = median(idealValues);
    results.(metric).mean.real = mean(realValues);
    results.(metric).median.real = median(realValues);
    results.(metric).median.overground = median(overgroundValues);
    results.(metric).mean.overground = mean(overgroundValues);
    % Display results
    fprintf('\n%s:\n', metric);
    fprintf('%s %s %s\n', 'Ideal', 'Real', 'Overground');

    fprintf('%f %f %f\n', results.(metric).mean.ideal, results.(metric).mean.real, results.(metric).mean.overground);

end

% Get number of muscles (assuming 18 muscles per trial)
numMuscles = size(resTable.MetCostPerMus{1}, 1);

% Create muscle index labels (assuming ordered muscle names)
muscleNames = arrayfun(@(x) sprintf('Muscle_%d', x), 1:numMuscles, 'UniformOutput', false);

% Initialize result tables
metCostPerMus_mean = table('Size', [numMuscles, 4], ...
                           'VariableTypes', {'string', 'double', 'double', 'double'}, ...
                           'VariableNames', {'Muscle', 'Mean_Ideal', 'Mean_Real', 'Mean_Overground'});



% Extract MetCostPerMus data from resTable
idealValues = cat(3, resTable.MetCostPerMus{idealIdx}); % Stack along 3rd dimension
realValues = cat(3, resTable.MetCostPerMus{realIdx});   % Stack along 3rd dimension
overgroundValues = cat(3, resTable.MetCostPerMus{overgroundIdx});   % Stack along 3rd dimension

% Compute mean & median per muscle
meanIdeal = mean(idealValues, 3);
medianIdeal = median(idealValues, 3);
meanReal = mean(realValues, 3);
medianReal = median(realValues, 3);
medianOverground = median(overgroundValues, 3);
meanOverground = mean(overgroundValues, 3);
% Fill result tables
for i = 1:numMuscles
    metCostPerMus_mean.Muscle(i) = muscleNames{i};
    metCostPerMus_mean.Mean_Ideal(i) = meanIdeal(i);
    metCostPerMus_mean.Mean_Real(i) = meanReal(i);
    metCostPerMus_mean.Mean_Overground(i) = meanOverground(i);

  %  metCostPerMus_median.Muscle(i) = muscleNames{i};
   % metCostPerMus_median.Median_Ideal(i) = medianIdeal(i);
   % metCostPerMus_median.Median_Real(i) = medianReal(i);
end

% Display tables
disp('Mean MetCost Per Muscle:');
disp(metCostPerMus_mean);






