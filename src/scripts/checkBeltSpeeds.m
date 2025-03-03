folderReal = 'C:\Users\Sophie Fleischmann\Documents\ResearchProjects\BeReal\ISB\BioMAC-Sim-Toolbox\results\TCSG\real_from_overgroundHighEffort';
speeds = [];

for i = 1:5
    resultFile = [folderReal, filesep, '2025_02_27_script2D_BeReal12_', num2str(i)];
    res = load(resultFile);
    result = res.result;
    bR = result.problem.idx.belt_right;
    compSpeedR = result.X(bR);
    speeds = cat(3, speeds, compSpeedR);
        
   
       

end 

meanSpeed = mean(speeds, 3);
figure
plot(meanSpeed)