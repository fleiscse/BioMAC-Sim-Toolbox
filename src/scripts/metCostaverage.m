clear; clc; close all;

%% Settings
dataFolder = 'results/TCSG_valid';  % top-level folder with subfolders
speedLabels = {'1','12', '14', '16', '18'};
metrics = {'MetCost', 'MetRate', 'CoT'};

resVar = cell(length(speedLabels) * 2, 6);
columnNames = {'Speed', 'Type', 'MetCost', 'MetCostPerMus', 'MetRate', 'CoT'};

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

        metCostList = [];
        metRateList = [];
        CoTList = [];
        metCostPerMusList = [];

        for f = 1:length(fileList)
            fileName = fullfile(speedFolder, fileList(f).name);
            res = load(fileName);
            result = res.result;

            try
                [metCost, ~, metCostPerMus, metRate, CoT] = result.problem.getMetabolicCost(result.X);
                metCostList(end+1) = metCost;
                metRateList(end+1) = metRate;
                CoTList(end+1) = CoT;
                metCostPerMusList(:,end+1) = metCostPerMus;  % [muscle x sample]
            catch
                warning('Failed to process file: %s', fileName);
            end
        end

        if isempty(metCostList)
            warning('No valid data found for %s at speed %s', condition, speedLabel);
            continue;
        end

        % Average metrics
        resVar{row, 1} = str2double(strrep(speedLabel, 'p', '.'));
        resVar{row, 2} = char(condition);
        resVar{row, 3} = mean(metCostList);
        resVar{row, 4} = mean(metCostPerMusList, 2);  % mean per muscle
        resVar{row, 5} = mean(metRateList);
        resVar{row, 6} = mean(CoTList);
        row = row + 1;
    end
end

%% Continue with existing sections
resTable = cell2table(resVar, 'VariableNames', columnNames);
%% Plot % of Realistic MetRate relative to Ideal
% Your custom label order
customOrder = [1, 12, 14, 16, 18];
xLabels = { '1','1.2', '1.4', '1.6', '1.8'};
total = 0

% Prepare y data in this same order
percentRealOverIdeal = zeros(length(customOrder), 1);
for i = 1:length(customOrder)
    spd = customOrder(i);
    idealRate = resTable.MetRate(strcmp(resTable.Type, 'ideal') & resTable.Speed == spd);
    realRate  = resTable.MetRate(strcmp(resTable.Type, 'realistic') & resTable.Speed == spd);

    if isempty(idealRate) || isempty(realRate)
        percentRealOverIdeal(i) = NaN;
    else
        percentRealOverIdeal(i) = (realRate - idealRate);
        total=total + percentRealOverIdeal(i);
    end
end
mean = total/5

% Plot with categorical x-axis
figure; hold on; box off;
x = 1:length(customOrder);
plot(x, percentRealOverIdeal, '-o', 'LineWidth', 2, 'Color', [0.2 0.4 0.7]);
yline(1, '--k', 'LineWidth', 1.5);

xticks(x);
xticklabels(xLabels);
xlabel('Speed (m/s)');
ylabel('Realistic MetRate (% of Ideal)');
title('Relative Metabolic Rate: Realistic vs Ideal');
ylim([min(percentRealOverIdeal)-0.05, max(percentRealOverIdeal)+0.05]);
grid on;
set(gca, 'FontSize', 14);


%% Plot: MetRate for Ideal vs Realistic using custom speed order (no lines)
fig=figure; hold on; box off;

customOrder = [1, 12, 14, 16, 18];
xLabels = { '1', '1.2', '1.4', '1.6', '1.8' };

idealRates = NaN(size(customOrder));
realRates  = NaN(size(customOrder));

for i = 1:length(customOrder)
    spd = customOrder(i);

    idealIdx = strcmp(resTable.Type, 'ideal') & resTable.Speed == spd;
    realIdx  = strcmp(resTable.Type, 'realistic') & resTable.Speed == spd;

    if any(idealIdx)
        idealRates(i) = resTable.MetRate(idealIdx);
    end
    if any(realIdx)
        realRates(i) = resTable.MetRate(realIdx);
    end
end

% Plot as scatter points only (no connecting lines)
x = 1:length(customOrder);
scatter(x , idealRates, 80, 'filled', 'DisplayName', 'Ideal', 'MarkerFaceColor','#0072BD')       % slight left shift
scatter(x, realRates, 80, 'filled', 'DisplayName', 'Realistic', 'MarkerFaceColor', '#EDA120');   % slight right shift

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

% %% Plot: MetRate for Ideal vs Realistic with fitted curves (no connecting lines)
% figure; hold on; box off;
% 
% customOrder = [1, 12, 14, 16, 18];  % original speed labels
% xLabels = { '1', '1.2', '1.4', '1.6', '1.8' };
% realSpeeds = [1.0, 1.2, 1.4, 1.6, 1.8];  % actual m/s values for plotting
% 
% idealRates = NaN(size(customOrder));
% realRates  = NaN(size(customOrder));
% 
% for i = 1:length(customOrder)
%     spd = customOrder(i);
% 
%     idealIdx = strcmp(resTable.Type, 'ideal') & resTable.Speed == spd;
%     realIdx  = strcmp(resTable.Type, 'realistic') & resTable.Speed == spd;
% 
%     if any(idealIdx)
%         idealRates(i) = resTable.MetRate(idealIdx);
%     end
%     if any(realIdx)
%         realRates(i) = resTable.MetRate(realIdx);
%     end
% end
% 
% % Scatter points
% scatter(realSpeeds , idealRates, 80, 'filled', 'DisplayName', 'Ideal');       % shift left
% scatter(realSpeeds , realRates, 80, 'filled', 'DisplayName', 'Realistic');    % shift right
% 
% % Fit and plot 2nd-order polynomials
% fitSpeeds = linspace(min(realSpeeds), max(realSpeeds), 100);  % smooth x for fit line
% 
% % Ideal fit
% validIdeal = ~isnan(idealRates);
% p_ideal = polyfit(realSpeeds(validIdeal), idealRates(validIdeal), 2);
% y_fit_ideal = polyval(p_ideal, fitSpeeds);
% plot(fitSpeeds, y_fit_ideal, '-', 'LineWidth', 2, 'DisplayName', 'Ideal Fit');
% 
% % Realistic fit
% validReal = ~isnan(realRates);
% p_real = polyfit(realSpeeds(validReal), realRates(validReal), 2);
% y_fit_real = polyval(p_real, fitSpeeds);
% plot(fitSpeeds, y_fit_real, '-', 'LineWidth', 2, 'DisplayName', 'Realistic Fit');
% 
% % Aesthetics
% xlabel('Speed (m/s)');
% ylabel('Metabolic Rate (W/kg)');
% title('Metabolic Rate vs Speed (with Polynomial Fit)');
% xticks(realSpeeds);
% xticklabels(xLabels);
% legend('Location', 'northwest');
% grid on;
% set(gca, 'FontSize', 12);
