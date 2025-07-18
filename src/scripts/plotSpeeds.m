speeds = load("data/Walking/experimentalSpeeds.mat");

right12 = (speeds.speed_12).';% -1.2;
right18 = (speeds.speed_18).';
%xq = linspace(1, 65, 100);

%fy12 = fy_speeds.speed_12;
%idx = find(fy12 < 10, 1, 'first');
%speed12 = interp1((1:idx), right12(1:idx), xq, 'linear');

fig = figure;
fig.Position(4) = 300;
%fy18 = fy_speeds.speed_18;
%idx2 = find(fy18 < 10, 1, 'first');
%speed18 = interp1((1:idx2), right18(1:idx2), xq, 'linear');

subplot(1,2,1);
plot((1:100), right12, '--', 'Color', "#0072BD", 'LineWidth', 1 )

hold on 


folderReal = 'results/TCSG_valid/12';
speeds = [];


filePattern = fullfile(folderReal, ['*',  'realistic1', '.mat']);
fileList = dir(filePattern);
resultFile = fullfile(folderReal, fileList(1).name);

res = load(resultFile);
result = res.result;
bR = result.problem.idx.belt_right;
compSpeedR = result.X(bR);

plot(compSpeedR, 'Color', "#0072BD", 'LineWidth', 1 )
set(gca,'box','off')
legend('measured', 'simulated', 'Location', 'southeast')
xlabel('% of stride')
ylabel('belt speed (m/s)')
ylim([1.155 1.22])
folderReal = 'results/TCSG_valid/18';
subplot(1,2,2)
plot((1:100), right18, '--','Color',  "#D95319", 'LineWidth', 1 )
filePattern = fullfile(folderReal, ['*', 'realistic2', '.mat']);
fileList = dir(filePattern);
resultFile = fullfile(folderReal, fileList(1).name);

res = load(resultFile);
result = res.result;
bR = result.problem.idx.belt_right;
compSpeedR = result.X(bR);
hold on

plot(compSpeedR, 'Color', "#D95319", 'LineWidth', 1 )


xlabel('% of stride')
ylabel('belt speed (m/s)')

legend('measured', 'simulated', 'Location', 'southeast')
set(gca,'box','off')
fontsize(gcf,scale=1.2)
saveas(fig, 'beltspeeds.svg')