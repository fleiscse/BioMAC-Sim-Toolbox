fy_speeds = load("data\Walking\grfy.mat");
speeds = load("data\Walking\experimentalSpeeds.mat");

right12 = (speeds.speed_12).';% -1.2;
right18 = (speeds.speed_18).';
xq = linspace(1, 65, 100);

fy12 = fy_speeds.speed_12;
idx = find(fy12 < 10, 1, 'first');
speed12 = interp1((1:idx), right12(1:idx), xq, 'linear');

fig = figure;
fig.Position(4) = 300;
fy18 = fy_speeds.speed_18;
idx2 = find(fy18 < 10, 1, 'first');
speed18 = interp1((1:idx2), right18(1:idx2), xq, 'linear');

subplot(1,2,1);
plot((1:100), right12, '--', 'Color', "#0072BD", 'LineWidth', 1 )

hold on 
%plot((1:100), right18, '--','Color',  "#D95319", 'LineWidth', 1 )

folderReal = 'C:\Users\Sophie Fleischmann\Documents\ResearchProjects\BeReal\ISB\BioMAC-Sim-Toolbox\results\TCSG\real_from_overground';
speeds = [];

for i = 1:5
    filePattern = fullfile(folderReal, ['*', num2str(i), '.mat']);
    fileList = dir(filePattern);
    resultFile = fullfile(folderReal, fileList(1).name);
    
    res = load(resultFile);
    result = res.result;
    bR = result.problem.idx.belt_right;
    compSpeedR = result.X(bR);
    %%plot(compSpeedR)
    speeds = cat(3, speeds, compSpeedR);

%     settings.GRF   = {'GRF_x_r', 'GRF_y_r'};
%     res = result.problem.extractData(result.X, settings, [], 1);
%     idxFy = find(strcmp(res.name, 'GRF_y_r') ==1);
%     Fy = res.sim{idxFy} ;
end

meanSpeed = mean(speeds, 3) ;
plot(meanSpeed, 'Color', "#0072BD", 'LineWidth', 1 )
set(gca,'box','off')
legend('measured', 'simulated', 'Location', 'southeast')
xlabel('% of stride')
ylabel('belt speed (m/s)')
subplot(1,2,2)
folderReal18 = 'C:\Users\Sophie Fleischmann\Documents\ResearchProjects\BeReal\ISB\BioMAC-Sim-Toolbox\results\TCSG\18\real_from_overground6';
speeds18 = [];       
for i = 1:5
    filePattern = fullfile(folderReal18, ['*', num2str(i), '.mat']);
    fileList = dir(filePattern);
    resultFile = fullfile(folderReal18, fileList(1).name);
    
    res = load(resultFile);
    result = res.result;
    bR = result.problem.idx.belt_right;
    compSpeedR = result.X(bR);
   % plot(compSpeedR)
    %hold on 

    speeds18 = cat(3, speeds18, compSpeedR);
       

end 

idx = find(meanSpeed > 1.2, 1, 'last');

meanSpeed18 = mean(speeds18, 3);
idx18 = find(meanSpeed18 > 1.8, 1, 'last');
plot((1:100), right18, '--','Color',  "#D95319", 'LineWidth', 1 )
hold on
sim12 = interp1((1:idx), meanSpeed(1:idx), xq, 'linear') - 1.2;
plot(meanSpeed18,'Color',  "#D95319", 'LineWidth', 1 )
 


sim18 = interp1((1:idx18), meanSpeed18(1:idx18), xq, 'linear') - 1.8;

xlabel('% of stride')
ylabel('belt speed (m/s)')

legend('measured', 'simulated', 'Location', 'southeast')
set(gca,'box','off')
fontsize(gcf,scale=1.2)
saveas(fig, 'beltspeeds.png')