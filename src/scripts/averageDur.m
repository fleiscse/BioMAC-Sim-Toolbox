clear; clc; close all;

%% Settings
dataFolder = 'results/TCSG_valid';  % top-level folder with subfolders
speedLabels = {'1','12', '14', '16', '18'};
metrics = {'MetCost', 'MetRate', 'CoT'};

resVar = cell(length(speedLabels) * 2, 2);
columnNames = {'Speed', 'Type', 'dur'};

row = 1;
for i = 1:length(speedLabels)
    speedLabel = speedLabels{i};
    speedFolder = fullfile(dataFolder, speedLabel);

    for condition = ["ideal", "realistic"]
        % Find all relevant files
        fileList = dir(fullfile(speedFolder, condition + "*.mat"));

        if isempty(fileList)
            warning('No %s files found in %s', condition, speedFolder);
            continue;
        end

        durList = [];
       

        for f = 1:length(fileList)
            fileName = fullfile(speedFolder, fileList(f).name);
            res = load(fileName);
            result = res.result;

            try
                if condition =="realistic"
                    d = result.X(9089); %9089 for real
                else
                    d = result.X(8889);
                   end
               % if d
                durList(end+1) = d;
               
            catch
                warning('Failed to process file: %s', fileName);
            end
        end

        if isempty(durList)
            warning('No valid data found for %s at speed %s', condition, speedLabel);
            continue;
        end

        % Average metrics
        resVar{row, 1} = str2double(strrep(speedLabel, 'p', '.'));
        resVar{row, 2} = char(condition);
        resVar{row, 3} = mean(durList);
        row = row + 1;
    end
end

%% Continue with existing sections
resTable = cell2table(resVar, 'VariableNames', columnNames);
%% Plot % of Realistic MetRate relative to Ideal



%% Plot: MetRate for Ideal vs Realistic using custom speed order (no lines)
fig=figure; hold on; box off;

customOrder = [1, 12, 14, 16, 18];
xLabels = { '1 m/s', '1.2m/s', '1.4m/s', '1.6m/s', '1.8m/s' };
exp = [1.170,1.122,  1.048, 1.009, 0.940];
ideal = resTable.dur(1:2:10)
real = resTable.dur(2:2:10)
x = 1:length(customOrder);
scatter(x , ideal, 80, 'filled', 'DisplayName', 'Ideal', 'MarkerFaceColor','#0072BD')       % slight left shift
scatter(x, real, 80, 'filled', 'DisplayName', 'Realistic', 'MarkerFaceColor', '#EDA120');   % slight right shift
%scatter(x , exp, 80, 'filled', 'DisplayName', 'Experimental', 'MarkerFaceColor',[128, 128, 128]/255)       % slight left shift

xticks(x);
xticklabels(xLabels);
%xlabel('Speed (m/s)');
ylabel('stride duration (s)');
%title('Metabolic Rate vs Speed (Ideal vs Realistic)');
%legend('Location', 'southwest');
set(gcf, 'Position', [100, 100, 600, 250]);

grid off;
set(gca, 'FontSize', 14);
saveas(fig, 'duration.svg')


%% Plot: MetRate for Ideal vs Realistic using custom speed order (no lines)
fig=figure; hold on; box off;

customOrder = [1, 12, 14, 16, 18];
xLabels = { '1', '1.2', '1.4', '1.6', '1.8' };
exp = [1.170,1.122,  1.048, 1.009, 0.940];
ideal = resTable.dur(1:2:10)
real = resTable.dur(2:2:10)
x = 1:length(customOrder);
scatter(x , real - ideal, 80, 'filled', 'DisplayName', 'Ideal', 'MarkerFaceColor','#0072BD')       % slight left shift
%scatter(x, real, 80, 'filled', 'DisplayName', 'Realistic', 'MarkerFaceColor', '#EDA120');   % slight right shift

xticks(x);
xticklabels(xLabels);
xlabel('Speed (m/s)');
ylabel('Metabolic Rate (W/kg)');
%title('Metabolic Rate vs Speed (Ideal vs Realistic)');
legend('Location', 'northwest');
set(gcf, 'Position', [100, 100, 400, 400]);

grid off;
set(gca, 'FontSize', 14);
saveas(fig, 'metRate2.svg')