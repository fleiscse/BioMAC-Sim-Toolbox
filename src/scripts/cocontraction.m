clear; clc;

% Settings
dataFolder = 'results/TCSG_valid';
speedLabels = {'1', '12', '14', '16', '18'};
xLabels = {'1.0', '1.2', '1.4', '1.6', '1.8'};
nSpeeds = length(speedLabels);

% Relevant ankle muscles
muscleNames = {'gastroc_r', 'soleus_r', 'tib_ant_r'};

% Storage: CCI across speeds
CCI.ideal = NaN(1, nSpeeds);
CCI.real  = NaN(1, nSpeeds);

for i = 1:nSpeeds
    speedFolder = fullfile(dataFolder, speedLabels{i});

    idealFiles = dir(fullfile(speedFolder, 'ideal*.mat'));
    realFiles  = dir(fullfile(speedFolder, 'realistic*.mat'));

    if isempty(idealFiles) || isempty(realFiles)
        warning('Missing trials at speed %s', speedLabels{i});
        continue;
    end

    nTrials = min(length(idealFiles), length(realFiles));
    settings.a = muscleNames;
% Initialize storage for CoA
CoA.ideal = NaN(1, nSpeeds);
CoA.real  = NaN(1, nSpeeds);

for i = 1:nSpeeds
    speedFolder = fullfile(dataFolder, speedLabels{i});
    idealFiles = dir(fullfile(speedFolder, 'ideal*.mat'));
    realFiles  = dir(fullfile(speedFolder, 'realistic*.mat'));

    if isempty(idealFiles) || isempty(realFiles)
        warning('Missing data for speed %s', speedLabels{i});
        continue;
    end

    nTrials = min(length(idealFiles), length(realFiles));
    settings.a = {'gastroc_r', 'tib_ant_r'};

    % Trial-wise storage
    CoA_trials_ideal = NaN(1, nTrials);
    CoA_trials_real  = NaN(1, nTrials);

    for k = 1:nTrials
        for condIdx = 1:2
            cond = ["ideal", "realistic"];
            fileList = {idealFiles, realFiles};

            try
                res = load(fullfile(speedFolder, fileList{condIdx}(k).name));
                result = res.result;
                extracted = result.problem.extractData(result.X, settings, [], 1);

                % Get activations
                a_SOL = extracted.sim{strcmp(extracted.name, 'gastroc_r') & strcmp(extracted.type, 'a')};
                a_TA  = extracted.sim{strcmp(extracted.name, 'tib_ant_r') & strcmp(extracted.type, 'a')};

                % Compute overlap area (approximate integral via sum)
                min_overlap = min(a_SOL, a_TA);
                overlap_area = sum(min_overlap);  % time-normalized → unitless area

                % Normalize by stride duration (assume 100 samples = 1 stride)
                phase_duration = length(a_SOL);  % could also use 100 if known

                coact_index = overlap_area / phase_duration;

                if cond(condIdx) == "ideal"
                    CoA_trials_ideal(k) = coact_index;
                else
                    CoA_trials_real(k)  = coact_index;
                end

            catch
                warning('Failed trial %d (%s) at speed %s', k, cond(condIdx), speedLabels{i});
            end
        end
    end

    % Average across trials
    CoA.ideal(i) = mean(CoA_trials_ideal, 'omitnan');
    CoA.real(i)  = mean(CoA_trials_real, 'omitnan');
end
end