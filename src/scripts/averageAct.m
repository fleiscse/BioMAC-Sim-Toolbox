clear; clc;

% Settings
dataFolder = 'results/TCSG_valid';
speedLabels = {'1', '12', '14', '16', '18'};
xLabels = {'1.0', '1.2', '1.4', '1.6', '1.8'};
nSpeeds = length(speedLabels);

% Muscle list (same as you previously used)
muscleNames = {'hamstrings_r', 'bifemsh_r', 'glut_max_r', 'iliopsoas_r', ...
               'rect_fem_r', 'vasti_r', 'gastroc_r', 'soleus_r', 'tib_ant_r'};
nMuscles = length(muscleNames);

% Storage: [muscle x speed]
peakIdeal = NaN(nMuscles, nSpeeds);
peakReal  = NaN(nMuscles, nSpeeds);

for i = 1:nSpeeds
    speedFolder = fullfile(dataFolder, speedLabels{i});

    idealFiles = dir(fullfile(speedFolder, 'ideal*.mat'));
    realFiles  = dir(fullfile(speedFolder, 'realistic*.mat'));

    if isempty(idealFiles) || isempty(realFiles)
        warning('Missing trials for speed %s', speedLabels{i});
        continue;
    end

    % TEMP storage for peaks across trials
    idealPeaks = NaN(nMuscles, length(idealFiles));
    realPeaks  = NaN(nMuscles, length(realFiles));

    % Ideal trials
    for k = 1:length(idealFiles)
        try
            res = load(fullfile(speedFolder, idealFiles(k).name));
            settings.a = muscleNames;
            extracted = res.result.problem.extractData(res.result.X, settings, [], 1);

            for m = 1:nMuscles
                idx = strcmp(extracted.name, muscleNames{m}) & strcmp(extracted.type, 'a');
                if any(idx)
                    idealPeaks(m, k) = max(extracted.sim{idx});
                end
            end
        catch
            warning('Ideal trial %d failed at speed %s', k, speedLabels{i});
        end
    end

    % Realistic trials
    for k = 1:length(realFiles)
        try
            res = load(fullfile(speedFolder, realFiles(k).name));
            settings.a = muscleNames;
            extracted = res.result.problem.extractData(res.result.X, settings, [], 1);

            for m = 1:nMuscles
                idx = strcmp(extracted.name, muscleNames{m}) & strcmp(extracted.type, 'a');
                if any(idx)
                    realPeaks(m, k) = max(extracted.sim{idx});
                end
            end
        catch
            warning('Realistic trial %d failed at speed %s', k, speedLabels{i});
        end
    end

    % Store average peak per muscle
    peakIdeal(:, i) = mean(idealPeaks, 2, 'omitnan');
    peakReal(:, i)  = mean(realPeaks, 2, 'omitnan');
end

% Compute difference: realistic - ideal
deltaPeak = (peakReal - peakIdeal);

% Output table
resultTable = array2table(deltaPeak, 'VariableNames', xLabels, 'RowNames', muscleNames);
disp('Δ Peak Muscle Activation (Realistic - Ideal):');
disp(resultTable);

% Optional: save
writetable(resultTable, 'peak_activation_difference.csv', 'WriteRowNames', true);
%% Plot full activation curves over time (mean of trials)
for i = 3:3
    speedFolder = fullfile(dataFolder, speedLabels{i});
    idealFiles = dir(fullfile(speedFolder, 'ideal*.mat'));
    realFiles  = dir(fullfile(speedFolder, 'realistic*.mat'));

    if isempty(idealFiles) || isempty(realFiles)
        continue;
    end

    % Initialize: [muscle x time x trials]
    nTrials = min(length(idealFiles), length(realFiles));
    settings.a = muscleNames;

    % Get time vector from first trial
    temp = load(fullfile(speedFolder, idealFiles(1).name));
    extracted = temp.result.problem.extractData(temp.result.X, settings, [], 1);
    firstIdx = find(strcmp(extracted.type, 'a'), 1);
    timeVec = linspace(0, 100, length(extracted.sim{firstIdx}));  % % gait cycle

    nTime = length(timeVec);
    A_ideal = NaN(nMuscles, nTime, nTrials);
    A_real  = NaN(nMuscles, nTime, nTrials);

    % Loop over trials
    for k = 1:nTrials
        % Ideal
        try
            res = load(fullfile(speedFolder, idealFiles(k).name));
            extracted = res.result.problem.extractData(res.result.X, settings, [], 1);
            for m = 1:nMuscles
                idx = strcmp(extracted.name, muscleNames{m}) & strcmp(extracted.type, 'a');
                if any(idx)
                    A_ideal(m, :, k) = extracted.sim{idx};
                end
            end
        catch
            warning('Failed ideal %d at speed %s', k, speedLabels{i});
        end

        % Realistic
        try
            res = load(fullfile(speedFolder, realFiles(k).name));
            extracted = res.result.problem.extractData(res.result.X, settings, [], 1);
            for m = 1:nMuscles
                idx = strcmp(extracted.name, muscleNames{m}) & strcmp(extracted.type, 'a');
                if any(idx)
                    A_real(m, :, k) = extracted.sim{idx};
                end
            end
        catch
            warning('Failed realistic %d at speed %s', k, speedLabels{i});
        end
    end

    % Plot
    fig = figure('Name', ['Activations at speed ' xLabels{i}], 'NumberTitle', 'off');
    for m = 9:9
         hold on;

        meanIdeal = mean(squeeze(A_ideal(m, :, :)), 2, 'omitnan');
        meanReal  = mean(squeeze(A_real(m, :, :)), 2, 'omitnan');

        plot(timeVec, meanIdeal, '-', 'LineWidth', 1.5, 'DisplayName', 'Ideal');
        plot(timeVec, meanReal, '-',  'LineWidth', 1.5, 'DisplayName', 'Realistic');

        title(muscleNames{m}, 'Interpreter', 'none');
        xlim([0 100]); ylim([0 1]);
        xlabel('% Gait'); ylabel('Activation');
        grid off;
        if m == 1
            legend('Location', 'best');
        end
    end
  %  sgtitle(['Muscle Activations at Speed ' xLabels{i} ' m/s']);
    set(gca, 'FontSize', 16);
    saveas(fig, 'tib14.svg')
end

% Define joint pairs for CCI
jointNames = {'hip', 'knee', 'ankle'};
pairA = {'iliopsoas_r', 'hamstrings_r', 'tib_ant_r'};
pairB = {'glut_max_r',  'vasti_r',      'gastroc_r'};

% Initialize: [joint x speed]
CCI.ideal = NaN(length(jointNames), nSpeeds);
CCI.real  = NaN(length(jointNames), nSpeeds);

for i = 1:nSpeeds
    speedFolder = fullfile(dataFolder, speedLabels{i});
    idealFiles = dir(fullfile(speedFolder, 'ideal*.mat'));
    realFiles  = dir(fullfile(speedFolder, 'realistic*.mat'));

    if isempty(idealFiles) || isempty(realFiles)
        continue;
    end

    nTrials = min(length(idealFiles), length(realFiles));
    settings.a = muscleNames;

    % Store trial-wise CCI
    CCI_trial.ideal = NaN(length(jointNames), nTrials);
    CCI_trial.real  = NaN(length(jointNames), nTrials);

    for k = 1:nTrials
        for condIdx = 1:2
            cond = ["ideal", "realistic"];
            fileList = {idealFiles, realFiles};
            res = load(fullfile(speedFolder, fileList{condIdx}(k).name));
            result = res.result;
            extracted = result.problem.extractData(result.X, settings, [], 1);

            % Initialize trial activations
            A = containers.Map;
            for m = 1:nMuscles
                idx = strcmp(extracted.name, muscleNames{m}) & strcmp(extracted.type, 'a');
                if any(idx)
                    A(muscleNames{m}) = extracted.sim{idx};
                end
            end

            CCIvals = NaN(length(jointNames), 1);
            for j = 1:length(jointNames)
                if isKey(A, pairA{j}) && isKey(A, pairB{j})
                    a1 = A(pairA{j});
                    a2 = A(pairB{j});
                    cci_t = 2 * min(a1, a2) ./ (a1 + a2 + eps);  % eps avoids division by 0
                    CCIvals(j) = mean(cci_t, 'omitnan');
                end
            end

            if cond(condIdx) == "ideal"
                CCI_trial.ideal(:, k) = CCIvals;
            else
                CCI_trial.real(:, k) = CCIvals;
            end
        end
    end

    % Average across trials
    CCI.ideal(:, i) = mean(CCI_trial.ideal, 2, 'omitnan');
    CCI.real(:, i)  = mean(CCI_trial.real, 2, 'omitnan');
end


% Plot CCI curves per joint
figure;
hold on;
for j = 1:length(jointNames)
    plot(1:nSpeeds, CCI.ideal(j,:), '--o', 'LineWidth', 1.5, 'DisplayName', [jointNames{j} ' (Ideal)']);
    plot(1:nSpeeds, CCI.real(j,:), '-s', 'LineWidth', 1.5, 'DisplayName', [jointNames{j} ' (Realistic)']);
end
xticks(1:nSpeeds); xticklabels(xLabels);
xlabel('Speed (m/s)');
ylabel('CCI');
title('Co-Contraction Index Across Speeds');
legend('Location', 'best');
grid on;
set(gca, 'FontSize', 12);
