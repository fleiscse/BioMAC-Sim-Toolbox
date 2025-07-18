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
from_standing = 0;
% Get path of this script
filePath = fileparts(mfilename('fullpath'));
% Path to your repository
path2repo = [filePath filesep '..' filesep '..' filesep];

% Fixed settings
dataFolder     = 'data/Walking';                % Relative from the path of the repository
dataFile       = 'Winter_normal_var_18.mat';           % Running data from Fukuchi 2017 subject S001 with 3.5 m/s
modelFile      = 'gait2d.osim';                 % Name of the OpenSim model with lumbar joint locked modelFile      = 'gait10dof18musc.osim';      % Name of the base model from OpenSim 
resultFolder   = 'results/TCSG_improve/18';%'results/TCSG/ideal'; 
resultFolder2   = 'results/TCSG_improve/18/real_from_overground6'; 	        % Relative from the path of the repository
% Relative from the path of the repository

%% Initalization
% Get date

dateString = datestr(date, 'yyyy_mm_dd');

% Get absolute file names
resultFileStanding = [resultFolder,filesep,'standing'];



%for 1.2:
%for 1.8:
trackingData = TrackingData.loadStruct(dataFile);
% Get absolute file names
%resultFileWalkingIdealTreadmill = '/home/rzlin/ys64ofuj/BeReal/BioMAC-Sim-Toolbox/results/overground/12/overground1';
for i = 1:5
    resultFileWalkingIdealTreadmill = sprintf('results/TCSG_improve_from_scratch/1/overground%d', i);



    resultFileWalking = sprintf('results/TCSG_improve_from_scratch/1/realistic2_%d', i);
dataFile           = [pwd,filesep,dataFolder,  filesep,dataFile];


% Create resultfolder if it does not exist
if ~exist([path2repo,filesep,resultFolder], 'dir')
    mkdir([path2repo,filesep,resultFolder]);
end



targetSpeed = 0; % m/s

% Create and initia  lize an instance of the OpenSim 2D model class.
model = Gait2d_osim_tread(modelFile, 1.9, 100);
%model = Gait2d_osim_tread(modelFile);
singlespeed = 1; %Change to 0 for split-belt treadmill simulation
if singlespeed
    model.setTreadmillSpeed(1.0);
else
    speed.left = 1.8;
    speed.right = 1.8;
    model.setTreadmillSpeed(speed);
end

targetSpeedTreadmill = 1.0;

isSymmetric = 0;
if from_standing
    initialGuess = resultFileStanding;
else
    initialGuess = resultFileWalkingIdealTreadmill;
end

problemWalking = Treadmill.walking2D_BeReal(model, resultFileWalking, trackingData, targetSpeed, targetSpeedTreadmill, isSymmetric, initialGuess);

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
resultWalking.problem.plotBeltSpeed(resultWalking.X)
%resultWalking.report(settings, style, resultFileWalking);
%resultWalking.problem.getMetabolicCost(resultWalking.X);
end