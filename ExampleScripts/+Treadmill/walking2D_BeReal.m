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
function problem = walking2D_BeReal(model, resultfile, trackingData, targetSpeed, targetSpeedTreadmill, isSymmetric, initialGuess)

%% Fixed settings
% We can coose the number of collocation nodes.
nNodes = 100;   
% Most of the time we use backard euler for discretization which is encoded with 'BE'.
Euler = 'BE';
% We usually use the name of the resultfile for the name of the logfile
logfile = resultfile;
% We want to plot intermediate results during solving of the problem.
plotLog = 1;


%% Create collocation problem
problem = Collocation(model, nNodes, Euler, logfile, plotLog);


%% Add states and controls including their bounds and initial values to the problem
% Get upper and lower bounds of the model and resize it
xmin = repmat(model.states.xmin, 1, nNodes+1); %We use one more than the number of nodes since the last node (the +1) one is then used in the periodicity constraint
xmax = repmat(model.states.xmax, 1, nNodes+1);

% Adapt the bounds for the first node to start the movement at X = 0
xmax(model.extractState('q', 'pelvis_tx'), 1) = 0;
xmin(model.extractState('q', 'pelvis_tx'), 1) = 0; 

% Add states (initial values will be specified later)
problem.addOptimVar('states', xmin, xmax);

% Add the treadmill speed as optimization variable (initial values will be specified later)
problem.addOptimVar('belt_left', repmat(0.9*targetSpeedTreadmill,1,nNodes), repmat(1.05*targetSpeedTreadmill,1,nNodes)); %one speed at every node,  add extra node
problem.addOptimVar('belt_right', repmat(0.9*targetSpeedTreadmill,1,nNodes), repmat(1.05*targetSpeedTreadmill,1,nNodes)); %one speed at every node, add extra node

problem.addOptimVar('belt_change_r', repmat(targetSpeedTreadmill - 1.05*targetSpeedTreadmill,1,nNodes), repmat(targetSpeedTreadmill - 0.9*targetSpeedTreadmill, 1, nNodes)); %one speed at every node,  add extra node
problem.addOptimVar('belt_change_l', repmat(targetSpeedTreadmill - 1.05*targetSpeedTreadmill,1,nNodes), repmat(targetSpeedTreadmill - 0.9*targetSpeedTreadmill, 1, nNodes));%speed at every node, add extra node

% Add controls to the problem using the default bounds (initial values will be specified later)
problem.addOptimVar('controls',repmat(model.controls.xmin,1,nNodes+1), repmat(model.controls.xmax,1,nNodes+1));

% Add duration of the movement 
problem.addOptimVar('dur',0.2, 2);

% Add speed in x direction of the movement. We choose here targetspeed for
% the lower and upper bound to ensure that we have exactly this speed. You
% could also choose other bounds.
problem.addOptimVar('speed',targetSpeed, targetSpeed);

% After adding all the components to X, we can use a previous solution as
% initial guess. In this case, we use the standing solution we produced
% before.
problem.makeinitialguess(initialGuess); 



%% Add objective terms to the problem
Wtracking = 1;
trackingData.resampleData(nNodes);
% 
if model.speed_left==1.8
    var = trackingData.variables;

    grfx = load("data/Walking/grfx.mat");
    grfy = load("data/Walking/grfy.mat");

    rightX = (grfx.speed_18 / 100 / 9.81).';
    rightY = (grfy.speed_18 / 100 / 9.81).';

    first_half = rightX(1:50);   % First 50 samples
    second_half = rightX(51:100); % Last 50 samples
    leftX = [second_half; first_half];
    % 
    first_half = rightY(1:50);   % First 50 samples
    second_half = rightY(51:100); % Last 50 samples
    leftY = [second_half; first_half];

    var(4, 3) = {rightX}; % Ensure same column names);
    var{5, 3} = {rightY}; % Ensure same column names);

    var(11, 3) = {leftX}; % Ensure same column names);
    var{12, 3} = {leftY}; % Ensure same column names);
    % Append the new rows to the existing table
    trackingData = trackingData.setVariables(var);
end




GRFSignals   = {'GRF_x_r', 'GRF_y_r', 'GRF_x_l', 'GRF_y_l'};
AngleSignals = {'hip_flexion_r', 'knee_angle_r', 'ankle_angle_r', 'hip_flexion_l', 'knee_angle_l', 'ankle_angle_l'};
nGRF   = length(GRFSignals);
nAngle = length(AngleSignals);
nAll   = nGRF + nAngle;
problem.addObjective(@trackGRF   , Wtracking*nGRF/nAll  , trackingData.extractData('GRF', GRFSignals));
problem.addObjective(@trackAngles, 0.5*Wtracking*nAngle/nAll, trackingData.extractData('angle', AngleSignals));
%problem.addObjective(@impactTerm, 0.01);

Weffort = 6;
weightsType = 'equal'; 
exponent = 3; 
problem.addObjective(@effortTermMusclesAct, Weffort, weightsType, exponent); 

Wreg = 0.0001;
problem.addObjective(@regTermTreadmill, Wreg);


%% Add constraints to the problem

problem.addConstraint(@dynamicConstraintsBeReal,repmat(model.constraints.fmin,1,nNodes),repmat(model.constraints.fmax,1,nNodes))
problem.addConstraint(@periodicityConstraint,zeros(model.nStates+model.nControls,1),zeros(model.nStates+model.nControls,1),isSymmetric)

%problem.addConstraint(@treadSpeedPeriodicityConstraint,zeros(3,1),zeros(3,1),isSymmetric)
problem.addConstraint(@treadmillSpeedConstraintsSigmoidImproved_sliding_window,repmat([-0.00;-0.00;0;0],1,nNodes),repmat([0.00;0.00;0;0],1,nNodes))
%problem.derivativetest()
%%problem.addConstraint(@treadmillSpeedConstraints_add_var,repmat([-0.00;-0.00],1,nNodes),repmat([0.00;0.00],1,nNodes))

fprintf('passed test')


end
