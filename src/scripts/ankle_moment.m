clear; clc;

% Setup
dataFolder = 'results/TCSG_improve/best/';
speedLabels = {'1', '12', '14', '16', '18', '2'};
xLabels = {'1.0', '1.2', '1.4', '1.6', '1.8', '2.0'};
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


% Plot max moments
figure;
for j = 1:nJoints
    subplot(2, ceil(nJoints/2), j);
    hold on;
    bar(1:nSpeeds, maxMoment.ideal(j,:), 0.4, 'FaceColor', [0.2 0.4 1], 'DisplayName', 'Ideal');
    bar((1:nSpeeds)+0.4, maxMoment.realistic(j,:), 0.4, 'FaceColor', [1 0.2 0.2], 'DisplayName', 'Realistic');
    xticks(1:nSpeeds + 0.2);
    xticklabels(xLabels);
    title(['Max ', joints{j}], 'Interpreter', 'none');
    ylabel('Nm/kg'); xlabel('Speed (m/s)');
    box off; grid on;
    if j == 1
        legend('Location', 'best');
    end
end
sgtitle('Max Joint Moments: Realistic vs Ideal');


% Plot min moments
figure;
for j = 1:nJoints
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