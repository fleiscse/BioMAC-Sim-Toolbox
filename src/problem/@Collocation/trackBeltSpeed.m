%======================================================================
%> @file trackGRF.m
%> @brief Collocation function to track ground reaction forces
%> @details
%> Details: Collocation::trackGRF()
%>
%> @author Eva Dorschky
%> @date November, 2017
%======================================================================

%======================================================================
%> @brief Matlab function to track ground reaction forces
%>
%> @details 
%> Computes difference between simulated and measured ground
%> reaction forces.
%> Supports data in bodyweight ('BW'), bodyweight percentage ('BW%'), and newton ('N').
%>
%> @param   obj     Collocation class object
%> @param   option  String parsing the demanded output: 'objval' or 'gradient'
%>                  (or 'init' for initialization)
%> @param   X       Double array: State vector containing at least 'states' of the model
%> @param   data    TrackingData: All rows with type 'GRF' will be tracked
%>
%> @retval  output  Objective values for input option 'objval' or vector
%>                  with gradient for input option 'gradient'
%======================================================================
function output = trackBeltSpeed(obj,option,X,data)

fctname = 'trackBeltSpeed';

%% initalization
if strcmp(option,'init')

    % check input parameter
    if ~isfield(obj.idx,'states') % check whether controls are stored in X
        error('Model states are not stored in state vector X.')
    end

    % initialize some variables (faster to get it once)
    variables = data.variables(strcmp(data.variables.type, 'beltspeed'), :); % data with the correct type
    obj.objectiveInit.(fctname).nVars = height(variables); % total number of variables
    obj.objectiveInit.(fctname).measMean = cell2mat(variables.mean'); % mean measured data
  
    obj.objectiveInit.(fctname).measVar = ones(size(obj.objectiveInit.(fctname).measMean)); % use 1 to operate as if we are not dividing by the variance
   
       
    obj.objectiveInit.(fctname).names = variables.name; % names of variables


    % Return a dummy value
    output = NaN;
    return;
end


%% compute demanded output
% get variables from initalization (faster)
nVars = obj.objectiveInit.(fctname).nVars;
measMean = obj.objectiveInit.(fctname).measMean;
measVar = obj.objectiveInit.(fctname).measVar;
names = obj.objectiveInit.(fctname).names;

% initialize some more variables
x = X(obj.idx.states); % extract states
bL = obj.idx.belt_left;
bR = obj.idx.belt_right;


if strcmp(option,'objval') %objective value for tracking GRFs  
    output = 0;
    
    exp_speeds = zeros(obj.nNodes, 2);
   
    exp_speeds(:,1) = X(bL);
    exp_speeds(:,2) = X(bR);

    
    
    for iVar = 1:nVars
        switch names{iVar}
            case 'belt_r'
                simVar = exp_speeds(:, 2);
            case 'belt_l'
                simVar = exp_speeds(:, 1);
          
        end
        output = output + ...
            sum((simVar - measMean(:, iVar)).^2./measVar(:, iVar)) / nVars /obj.nNodes;
      
    end

elseif strcmp(option,'gradient') %gradient for tracking GRFs 
    % track GRF
    output = zeros(size(X));
    exp_speeds = zeros(obj.nNodes, 2);
    exp_speeds(:,1) = X(bL);
    exp_speeds(:,2) = X(bR);

    for iVar = 1:nVars
         switch names{iVar}
            case 'belt_r'
                simVar = exp_speeds(:, 2);
                 output(bR) = output(bR) +...
            2*(simVar - measMean(:, iVar))./measVar(:, iVar)  / nVars / obj.nNodes;
            case 'belt_l'
                simVar = exp_speeds(:, 1);
                output(bL) = output(bL) +...
            2*(simVar - measMean(:, iVar))./measVar(:, iVar)  / nVars / obj.nNodes;
        end
        
       
        
    end
   
else
    error('Unknown option');
end

end