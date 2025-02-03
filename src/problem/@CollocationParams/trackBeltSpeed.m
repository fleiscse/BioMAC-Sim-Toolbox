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
    if ~isfield(obj.idx,'belt_left') % check whether duration is stored in X
        error('beltspeed is not stored in state vector X.')
    end

    

    % Return a dummy value
    output = NaN;
    return;
end



%% compute demanded output
% get variables from initalization (faster)

measuredSpeed = data.';
if strcmp(option,'objval') %objective value for tracking duration  
    compSpeed = X(obj.idx.belt_left); %state variable
    output = mean((measuredSpeed - compSpeed).^2);
   
    
elseif strcmp(option,'gradient') %gradient for tracking duration 
    output = zeros(size(X));
    n = length(measuredSpeed);
    
    
    compSpeed = X(obj.idx.belt_left); %state variable
    grad = 2 / n * (measuredSpeed - compSpeed ) * (-1);
    output(obj.idx.belt_left) = grad;
else
    error('Unknown option');
end

end