%======================================================================
%> @file Treadmill/script2D.m
%> @brief Script to show how to simulate a treadmill and with the 2D 
%> OpenSim model
%> @example
%> @details
%> This is a script that you can use to see how the treadmill works. 
%>
%> @author Anne Koelewijn
%> @date August, 2024
%======================================================================

clear all 
close all
clc

%% Settings
% Get path of this script
filePath = fileparts(mfilename('fullpath'));
% Path to your repository
path2repo = [filePath filesep '..' filesep '..' filesep];

% Fixed settings
dataFolder     = 'data';                % Relative from the path of the repository
dataFile       = 'experimentalData.mat';           % Running data from Fukuchi 2017 subject S001 with 3.5 m/s
modelFile      = 'gait2d.osim';                 % Name of the OpenSim model with lumbar joint locked modelFile      = 'gait10dof18musc.osim';      % Name of the base model from OpenSim 
resultFolder   = 'results/Treadmill'; 	        % Relative from the path of the repository

%% Initalization
% Get date
dateString = datestr(date, 'yyyy_mm_dd');
dataFile = [path2repo,filesep,dataFolder,  filesep,dataFile];
resultFile  = [path2repo,filesep,resultFolder,filesep,dateString,'_', mfilename,'_params_idealTreadmill'];

% Create resultfolder if it does not exist
if ~exist([path2repo,filesep,resultFolder], 'dir')
    mkdir([path2repo,filesep,resultFolder]);
end


%% Simulate walking on a single belt or split-belt treadmill
% Load tracking data struct and create a TrackingData object
experimentalData  = load(dataFile);
speed = experimentalData.speeds.speed_12;
grfx = experimentalData.grfx.speed_12;
grfy = experimentalData.grfy.speed_12;



targetSpeedTreadmill = 1.2;
delay = 5;


problemWalking = Treadmill.params_BeReal(speed, grfx, grfy, targetSpeedTreadmill, delay, resultFile);

% Create solver and change solver settings
solver = IPOPT();
solver.setOptionField('max_iter', 5000);
solver.setOptionField('tol', 0.0005);

% Solve the optimization problem and save the result. 
resultWalking = solver.solve(problemWalking);
resultWalking.save(resultFileWalking); 
%investigateBeltControl(resultWalking)
% If you want to create plots, take a look at one of the other examples.
resultWalking.problem.writeMovie(resultWalking.X, resultWalking.filename);


settings.plotInitialGuess = 1;
style.figureSize = [0 0 16 26];
% When also giving a filename as input, the function will automatically
% create a pdf summarizing all information and plots. Take a look at it!
% You might have to press Enter a couple of times in the Command Window
% for the saving to continue.
%resultWalking.report(settings, style, resultFileWalking);
%resultWalking.problem.getMetabolicCost(resultWalking.X)