%======================================================================
%> @file TreadMillWalking/walking2D.m
%> @brief Function to specify the optimization problem for 2D running
%>
%> @author Anne Koelewijn
%> @date August 2024
%======================================================================

% ======================================================================
%> @brief Function to specify the optimization problem for 2D running
%>
%> @param   model          Gait2dc: Model which should be used for the simulation
%> @param   resultfile     String: Name of the resultfile including path
%> @param   trackingData   TrackingData: Tracking Data containing angles and GRFs data 
%> @param   targetSpeed    Double: Target speed of the movement in x direction in m/s. 
%>                         This target speed will be enforced.
%> @param   isSymmetric    Bool: Specifys weather we assume symmetry of the
%>                         movement. If we assume symmetry, we simulate only one half 
%>                         of gait cycle. This has to fit to the tracking data.
%> @param   initialGuess   String: Filename with path specifying the initial guess 
%> @retval  problem        Collocation: Optimization problem for 2D running
% ======================================================================
function problem = params_BeReal(speed, grfx, grfy, targetSpeedTreadmill, delay, resultfile)

%% Fixed settings
% We can choose the number of collocation nodes.
nNodes = length(speed);   
% Most of the time we use backard euler for discretization which is encoded with 'BE'.
Euler = 'BE';
% We usually use the name of the resultfile for the name of the logfile
logfile = resultfile;
% We want to plot intermediate results during solving of the problem.
plotLog = 0;


%% Create collocation problem
problem = CollocationParams(nNodes, Euler, logfile, plotLog);


%% Add states and controls including their bounds and initial values to the problem
% Get upper and lower bounds of the model and resize it


% Adapt the bounds for the first node to start the movement at X = 0


% Add states (initial values will be specified later)

% Add the treadmill speed as optimization variable (initial values will be specified later)
problem.addOptimVar('belt_left', repmat(0.9*targetSpeedTreadmill,1,nNodes), repmat(1.05*targetSpeedTreadmill,1,nNodes)); %one speed at every node, add extra node

%problem.addOptimVar('grf_delay', 2, 9);
% % problem.addOptimVar('Kfy', -0.005, -0.00004, -0.0002);
% % problem.addOptimVar('Kfx', 0.0001, 0.01, 0.0006695288);
% % problem.addOptimVar('Kp', 0.01, 1, 0.1558);
% % problem.addOptimVar('Kd', -0.1, -0.000001, -0.005853);

problem.addOptimVar('Kfy', -0.005, -0.00004, -0.19466999407124552);
problem.addOptimVar('Kfx', 0.00001, 0.1, 4.0933641847270306e-06);
problem.addOptimVar('Kp', 1, 6, 4.628766189441804);
problem.addOptimVar('Kd', -0.5, 0, -0.05614093069537616);
%problem.addOptimVar('c', 0.01, 0.01);


% Add speed in x direction of the movement. We choose here targetspeed for
% the lower and upper bound to ensure that we have exactly this speed. You
% could also choose other bounds.

% After adding all the components to X, we can use a previous solution as
% initial guess. In this case, we use the standing solution we produced
% before.



%% Add objective terms to the problem
W_track = 1;
problem.addObjective(@trackBeltSpeed, W_track, speed); 




%% Add constraints to the problem

lb = repmat(-0.1001,1,nNodes);
up = repmat(0.1001,1,nNodes);
%lb(80:100) = -0.5;
%up(80:100) = 0.5;
problem.addConstraint(@treadmillSpeedConstraintsParamsSigmoid,lb,up, grfx, grfy, delay)
problem.derivativetest()
fprintf('passed test')


end
