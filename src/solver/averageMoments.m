clear; clc;
clear all;
% Setup
dataFolder = 'results/TCSG_valid/';
speedLabels = {'1', '12', '14', '16', '18'};
xLabels = {'1.0', '1.2', '1.4', '1.6', '1.8'};
nSpeeds = length(speedLabels);

joints = {'hip_flexion_r', 'knee_angle_r', 'ankle_angle_r'};
nJoints = length(joints);

% Initialize result storage
maxMoment.ideal = NaN(nJoints, nSpeeds);
minMoment.ideal = NaN(nJoints, nSpeeds);
maxMoment.realistic = NaN(nJoints, nSpeeds);
minMoment.realistic = NaN(nJoints, nSpeeds);

settings.moment = joints;

% Loop over speeds and conditions
for i = 1:nSpeeds
    speedFolder = fullfile(dataFolder, speedLabels{i});

    for cond = ["ideal", "realistic"]
        files = dir(fullfile(speedFolder, cond + "*.mat"));
        if isempty(files)
            warning('No %s files found for speed %s', cond, speedLabels{i});
            continue;
        end

        % Temporary storage for this condition & speed
        tempMin = NaN(nJoints, length(files));
        tempMax = NaN(nJoints, length(files));

        for f = 1:length(files)
            try
                res = load(fullfile(speedFolder, files(f).name));
                result = res.result;
                extracted = result.problem.extractData(result.X, settings, [], 1);

                for j = 1:nJoints
                    idx = strcmp(extracted.name, joints{j}) & strcmp(extracted.type, 'moment');
                    if any(idx)
                        data = extracted.sim{idx};
                        tempMin(j,f) = min(data);
                        tempMax(j,f) = max(data);
                    end
                end
            catch
                warning('Failed to process file: %s', files(f).name);
            end
        end

        % Store average over trials
        if cond == "ideal"
            minMoment.ideal(:,i) = mean(tempMin, 2, 'omitnan');
            maxMoment.ideal(:,i) = mean(tempMax, 2, 'omitnan');
        else
            minMoment.realistic(:,i) = mean(tempMin, 2, 'omitnan');
            maxMoment.realistic(:,i) = mean(tempMax, 2, 'omitnan');
        end
    end
end



% Plot max moments
fig = figure;
for j = 2:2
   
    hold on 
    bar(1:nSpeeds, maxMoment.ideal(j,:), 0.4, 'FaceColor', '#0072BD', 'DisplayName', 'Ideal','EdgeColor','none' );
    bar((1:nSpeeds)+0.4, maxMoment.realistic(j,:), 0.4, 'FaceColor', '#EDA120', 'DisplayName', 'Realistic', 'EdgeColor','none');
    xticks(1:nSpeeds + 0.2);
    xticklabels(xLabels);
  
    ylabel('Knee extension moment (Nm)'); xlabel('Speed (m/s)');
    box off; grid off;
    set(gca, 'FontSize', 14);
    xlim([0.5 6])
    if j == 2
        legend('Location', 'northwest');
    end
end
set(gca, 'FontSize', 16);
saveas(fig, 'knee_extension_moment.svg')

% Plot max moments
fig = figure;
for j = 3:3
   
    hold on 
    bar(1:nSpeeds, maxMoment.ideal(j,:), 0.4, 'FaceColor', '#0072BD', 'DisplayName', 'Ideal','EdgeColor','none' );
    bar((1:nSpeeds)+0.4, maxMoment.realistic(j,:), 0.4, 'FaceColor', '#EDA120', 'DisplayName', 'Realistic', 'EdgeColor','none');
    xticks(1:nSpeeds + 0.2);
    xticklabels(xLabels);
  
    ylabel('Dorsiflexion moment (Nm)'); xlabel('Speed (m/s)');
    box off; grid off;
    set(gca, 'FontSize', 16);
    xlim([0.5 6])
    if j == 2
        legend('Location', 'northwest');
    end
end
saveas(fig, 'dorsiflexion_moment.svg')
%sgtitle('Max Joint Moments: Realistic vs Ideal');

% Plot min moments (unchanged)
figure;
for j = 1:3
    subplot(2, ceil(nJoints/2), j);
    hold on;
    bar(1:nSpeeds, minMoment.ideal(j,:), 0.4, 'FaceColor', [0.2 0.4 1], 'DisplayName', 'Ideal');
    bar((1:nSpeeds)+0.4, minMoment.realistic(j,:), 0.4, 'FaceColor', [1 0.2 0.2], 'DisplayName', 'Realistic');
    xticks(1:nSpeeds + 0.2);
    xticklabels(xLabels);
    title(['Min ', joints{j}], 'Interpreter', 'none');
    ylabel('Nm/kg'); xlabel('Speed (m/s)');
    box off; grid on;
    if j == 1
        legend('Location', 'best');
    end
end
sgtitle('Min Joint Moments: Realistic vs Ideal');


% Initialize storage for RMSE and ROM
RMSE = NaN(nJoints, nSpeeds);
ROM.ideal = NaN(nJoints, nSpeeds);
ROM.realistic = NaN(nJoints, nSpeeds);

settings.angle = joints;  % now extract angles too

for i = 1:nSpeeds
    speedFolder = fullfile(dataFolder, speedLabels{i});

    % Load all ideal and realistic trials
    idealFiles = dir(fullfile(speedFolder, "ideal*.mat"));
    realFiles  = dir(fullfile(speedFolder, "realistic*.mat"));

    if isempty(idealFiles) || isempty(realFiles)
        warning('Missing files for speed %s', speedLabels{i});
        continue;
    end

    % Match files by count or use first available pairing
    nTrials = min(length(idealFiles), length(realFiles));

    for j = 1:nJoints
        angleIdealAll = [];
        angleRealAll = [];

        for k = 1:nTrials
            try
                % Load both trials
                resIdeal = load(fullfile(speedFolder, idealFiles(k).name));
                resReal  = load(fullfile(speedFolder, realFiles(k).name));

                X_ideal = resIdeal.result.X;
                X_real  = resReal.result.X;

                % Extract joint angles
                extractedIdeal = resIdeal.result.problem.extractData(X_ideal, settings, [], 1);
                extractedReal  = resReal.result.problem.extractData(X_real, settings, [], 1);

                idxIdeal = strcmp(extractedIdeal.name, joints{j}) & strcmp(extractedIdeal.type, 'angle');
                idxReal  = strcmp(extractedReal.name, joints{j})  & strcmp(extractedReal.type, 'angle');

                if any(idxIdeal) && any(idxReal)
                    angleIdeal = extractedIdeal.sim{idxIdeal};
                    angleReal  = extractedReal.sim{idxReal};

                    % Store for group RMSE and ROM
                    angleIdealAll = [angleIdealAll; angleIdeal(:)'];
                    angleRealAll  = [angleRealAll; angleReal(:)'];
                end
            catch
                warning('Failed to extract joint angle for speed %s, joint %s, trial %d', ...
                        speedLabels{i}, joints{j}, k);
            end
        end

        % Only compute if valid trials available
        if ~isempty(angleIdealAll) && ~isempty(angleRealAll)
            % Mean signals across trials
            meanIdeal = mean(angleIdealAll, 1, 'omitnan');
            meanReal  = mean(angleRealAll, 1, 'omitnan');

            % Compute RMSE over gait cycle
            RMSE(j, i) = sqrt(mean((meanReal - meanIdeal).^2, 'omitnan'));

            % ROM
            ROM.ideal(j, i) = max(meanIdeal) - min(meanIdeal);
            ROM.realistic(j, i) = max(meanReal) - min(meanReal);
        end
    end
end
fprintf('\nRMSE (degrees or radians) between Ideal and Realistic Joint Angles:\n');
for j = 1:nJoints
    fprintf('%-15s: ', joints{j});
    fprintf('%6.3f ', RMSE(j, :));
    fprintf('\n');
end


figure;
for j = 1:nJoints
    subplot(2, ceil(nJoints/2), j);
    hold on;
    bar(1:nSpeeds, ROM.ideal(j,:), 0.4, 'FaceColor', [0.2 0.4 1], 'DisplayName', 'Ideal');
    bar((1:nSpeeds)+0.4, ROM.realistic(j,:), 0.4, 'FaceColor', [1 0.2 0.2], 'DisplayName', 'Realistic');
    xticks(1:nSpeeds + 0.2);
    xticklabels(xLabels);
    title(['ROM ', joints{j}], 'Interpreter', 'none');
    ylabel('Range (deg or rad)'); xlabel('Speed (m/s)');
    box off; grid on;
    if j == 1
        legend('Location', 'best');
    end
end
sgtitle('Range of Motion: Realistic vs Ideal');