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
dataFile       = 'Winter_normal'%;_var_18.mat';           % Running data from Fukuchi 2017 subject S001 with 3.5 m/s
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
resultFileWalkingIdealTreadmill = '/home/rzlin/ys64ofuj/BeReal/BioMAC-Sim-Toolbox/results/overground/18/overground5';
resultFileStanding = [path2repo,filesep,resultFolder,filesep,'standing'];


resultFileWalking  = '/home/rzlin/ys64ofuj/BeReal/BioMAC-Sim-Toolbox/results/TCSG_improve/18/sliding_window_normalVar/realistic5';
dataFile           = [pwd,filesep,dataFolder,  filesep,dataFile];


% Create resultfolder if it does not exist
if ~exist([path2repo,filesep,resultFolder], 'dir')
    mkdir([path2repo,filesep,resultFolder]);
end

%% Standing: Simulate standing with minimal effort without tracking data for one point in time (static)
% Create an instane of the OpenSim 2D model class using the default settings
if from_standing

    model = Gait2d_osim(modelFile, 1.9, 100);
%model.bodymass=100;
% model = Gait2dc(modelFile);

% Call IntroductionExamples.standing2D() to specify the optimizaton problem
% We use the same function as in IntroductionExamples, since we are solving
% the same problem.
    problemStanding = IntroductionExamples.standing2D(model, resultFileStanding);

% Create an object of class solver. We use most of the time the IPOPT here.
    solver = IPOPT();

% Change settings of the solver
    solver.setOptionField('tol', 0.0000001);
%solver.setOptionField('constr_viol_tol', 0.000001);

% Solve the optimization problem
    resultStanding = solver.solve(problemStanding);

% Save the result
    resultStanding.save(resultFileStanding);

% To plot the result we have to extract the states x from the result vector X
    x = resultStanding.X(resultStanding.problem.idx.states);

% Now, we can plot the stick figure visualizing the result
    figure();
    resultStanding.problem.model.showStick(x);
    title('2D Standing');
end

% If the model is standing on the toes, the optimization ends in a local optimum and not the global one. Rerun this
% section and you should find a different solution, due to a different random
% initial guess. You can run it a couple of times until you find a good
% solution, standing on flat feet.

%% Simulate walking on a single belt or split-belt treadmill
% Load tracking data struct and create a TrackingData object
trackingData = TrackingData.loadStruct(dataFile);

%beltspeed = load('data/Walking/experimentalSpeeds.mat');
%rightSpeed = beltspeed.speed_12;
%rightSpeed(65:100) = 1.2;
%first_half = rightSpeed(1:50);   % First 50 samples
%second_half = rightSpeed(51:100); % Last 50 samples

% Rearrange: Second half first, then first half
% leftSpeed = [second_half, first_half];
% 
% newRows = table(...
%     {'belt_r'; 'belt_l'}, ...    % Column 1: name
%     {'beltspeed'; 'beltspeed'}, ...          % Column 2: type (adjust if needed)
%     {rightSpeed.';leftSpeed.'}, ... % Column 3: mean (replace with actual data)
%     {ones(100,1); ones(100,1)}, ... % Column 4: var (replace with actual variance)
%     {'BW'; 'BW'}, ...            % Column 5: unit
%     'VariableNames', trackingData.variables.Properties.VariableNames ... % Ensure same column names
% );
% 
% % Append the new rows to the existing table
% trackingData = trackingData.setVariables([trackingData.variables; newRows]);

% The global speed is going to be zeros, as the movement now comes from the treadmill
targetSpeed = 0; % m/s

% Create and initia  lize an instance of the OpenSim 2D model class.
model = Gait2d_osim_tread(modelFile, 1.9, 100);
%model = Gait2d_osim_tread(modelFile);
singlespeed = 1; %Change to 0 for split-belt treadmill simulation
if singlespeed
    model.setTreadmillSpeed(1.8);
else
    speed.left = 1.8;
    speed.right = 1.8;
    model.setTreadmillSpeed(speed);
end

targetSpeedTreadmill = 1.8;

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
resultWalking.report(settings, style, resultFileWalking);
resultWalking.problem.getMetabolicCost(resultWalking.X);