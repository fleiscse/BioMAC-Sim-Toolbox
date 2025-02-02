%======================================================================
%> @file periodicityConstraint.m
%> @brief Collocation function to compute periodicity constraint
%> @details
%> Details: Collocation::periodicityConstraint()
%>
%> @author Eva Dorschky
%> @date November, 2017
%======================================================================

%======================================================================
%> @brief Computes constaint violation demanding periodic movement
%>
%> @param obj           Collocation class object
%> @param option        String parsing the demanded output
%> @param X             Double array: State vector containing at least 'states' and 'controls' of
%>                      the model and speed and duration of the periodic movement
%> @param sym           Boolean if movement is symmetric (half period is optimized) or not
%======================================================================
function output = treadSpeedPeriodicityConstraint(obj,option,X,sym)
fprintf('periodic')
%% check input parameter
if ~isfield(obj.idx,'states') || ~isfield(obj.idx,'controls') || ~isfield(obj.idx,'speed') || ~isfield(obj.idx,'dur') % check whether controls are stored in X
    error('State vector X does not contain required states.')
end

if size(X(obj.idx.states),2) ~= (obj.nNodes+1) || size(X(obj.idx.controls),2) ~= (obj.nNodes+1) %> @todo should we add another field to the state vector instead of expanding states to obj.nNodes+1
    error('Model states and controls need to be optimized at collocation nodes + 1')
end
    
%% compute demanded output
nStates = size(obj.idx.states,1);
nControls =  size(obj.idx.controls,1);
nNodes = obj.nNodes
icleft = 1;
icright = 2;

% duration
dur = X(obj.idx.dur);
% speed
speed = X(obj.idx.speed);

%
vBeltLeft = X(obj.idx.belt_left);
vBeltRight = X(obj.idx.belt_right);


if strcmp(option,'confun') %constraints of periodicity constraint
    output = 0;
    
   
    if ~sym
        % state must be periodic, but not necessarily symmetric
        output(icleft) = X(obj.idx.belt_left(:,end)) - X(obj.idx.belt_left(:,1));
        % controls must be periodic
        output(icright) = X(obj.idx.belt_right(:,end)) - X(obj.idx.belt_right(:,1));
        output(3) = 0;
    else
        % state must be periodic, but not necessarily symmetric
        output(icleft) = X(obj.idx.belt_left(:,end)) - X(obj.idx.belt_left(:,1));
        % controls must be periodic
        output(icright) = X(obj.idx.belt_right(:,end)) - X(obj.idx.belt_right(:,1));
        output(3) = X(obj.idx.belt_right(:,end)) - X(obj.idx.belt_left(:,1 + nNodes/2)); %is this correct??
    end
    
elseif strcmp(option,'jacobian') %jacobian of periodicity constraint
    output = spalloc(3,length(X),obj.Jnnz); %> @todo where to get Jnnz
    

    % compute derivatives of the displacement
    displacementxddur   = unitdisplacementx*speed(1);
    displacementxdspeed = unitdisplacementx*dur;
    if numel(speed) > 1
        displacementzddur   = unitdisplacementz*speed(2);
        displacementzdspeed = unitdisplacementz*dur;
    else
        displacementzddur   = 0;
        displacementzdspeed = [];
    end
    
    if ~sym
        output(icx,obj.idx.states(:,end))   = speye(nStates);
        output(icx,obj.idx.states(:,1))     = -speye(nStates);
        output(icx,obj.idx.dur)             = -displacementxddur-displacementzddur;
        output(icx,obj.idx.speed)           = [-displacementxdspeed, -displacementzdspeed];
        output(icu,obj.idx.controls(:,end)) = speye(nControls);
        output(icu,obj.idx.controls(:,1))   = -speye(nControls);
    else
        output(icx,obj.idx.states(:,end))        = speye(nStates);
        output(icx,obj.model.idxSymmetry.xindex) = -obj.model.idxSymmetry.xsign*ones(1,nStates).*speye(nStates);
        output(icx,obj.idx.dur)                  = -displacementxddur-displacementzddur;
        output(icx,obj.idx.speed)                = [-displacementxdspeed, -displacementzdspeed];
        output(icu,obj.idx.controls(:,end))      = speye(nControls);
        output(icu, obj.idx.controls(obj.model.idxSymmetry.uindex,1)) = -obj.model.idxSymmetry.usign*ones(1,nControls).*speye(nControls);
    end
    
else
    error('Unknown option');
end
end

