%======================================================================
%> @file @Collocation/treadmillSpeedConstraints.m
%> @brief Collocation function to compute the treadmill speed
%> @details
%> Details: Collocation::treadmillSpeedConstraints()
%>
%> @author Sophie Fleischmann
%> @date January, 2025
%======================================================================

%======================================================================
%> @brief Computes constaint violation demanding the real treadmill speed
%>
%> @param obj           Collocation class object
%> @param option        String parsing the demanded output
%> @param X             Double array: State vector containing at least 'states' of
%>                      the model and current treadmill speeds
%> @param sym           Boolean if movement is symmetric (half period is optimized) or not
%======================================================================
function output = treadmillSpeedConstraintsParamsSigmoid18_sliding_window(obj,option,X,grfx, grfy, delay)

%% check input parameter
if ~isfield(obj.idx,'belt_left') % check whether controls are stored in X
    error('Model states and left and right belt speeds need to be stored in state vector X.')
end


%% compute demanded output
nNodesDur = obj.nNodesDur + 1; %number of collocation nodes + 1

nconstraintspernode = 1; %one for left and one for right


%delay = X(obj.idx.grf_delay);
Kfy =  X(obj.idx.Kfy);
Kfx = X(obj.idx.Kfx);
Kp = X(obj.idx.Kp);
Kd = X(obj.idx.Kd);

c = 0.012;

if strcmp(option,'confun')
    output = zeros(nconstraintspernode*(nNodesDur-1),1);
    
    % dynamic equations must be zero
        for iNode=1:(nNodesDur-1)

        delayed_index = mod(iNode - delay-2, nNodesDur-1) +1;
        delayed_index2 = mod(iNode - delay+1, nNodesDur-1) +1;
        ic = (1:nconstraintspernode) +  (iNode-1)*nconstraintspernode; %indices of constraints of iNode in c
        
       
        fx2 = grfx(delayed_index2); %TODO: check if this is the same as when I add the 2 states
        fy2 = grfy(delayed_index2);
       
        fx1 = grfx(delayed_index); %TODO: check if this is the same as when I add the 2 states
        fy1 = grfy(delayed_index);

        fy_current= grfy(mod(iNode-1, nNodesDur-1) +1);
             
        v_curr = X(obj.idx.belt_left(iNode));
        v_prev = X(obj.idx.belt_left(mod(iNode - 2, nNodesDur-1)+1)); %only works if n_constraints (and not n+1) --> if I need n+1 points: use an if statement
       
        
        v_left = v_curr + Kfx *((fx2 - fx1)/(3*c)) + Kfy*((fy2 - fy1)/(3*c)) + Kp*(1.8 - v_curr) + Kd * ((-v_curr+ v_prev)/c);

        sigmoid_left = 0.0005 + 1 / (1 + exp(-50 * fy_current+10));
        
        v_left_next = X(obj.idx.belt_left(mod(iNode, nNodesDur - 1) + 1));

        %apply sigmoid: belt speed should be 1.2 if there is no vertical
        %force, else the computed speed
        v_left = sigmoid_left*v_left + (1-sigmoid_left) * 1.8;
        
        diff =  v_left- v_left_next ;
        output(ic) = diff;	% backward Euler discretization
        
        end
    

elseif strcmp(option,'jacobian')
    output = spalloc(nconstraintspernode*(nNodesDur-1),length(X),obj.Jnnz);


    idxVLeft = obj.idx.belt_left; % the indices of the left Belt speed in X
  %  idxDelay = obj.idx.grf_delay;
    idxKfy =  obj.idx.Kfy;
    idxKfx = obj.idx.Kfx;
    idxKp = obj.idx.Kp;
    idxKd = obj.idx.Kd;

   % idxC = obj.idx.c;
    
    
    for iNode = 1:(nNodesDur-1)

        %f_current= obj.model.getGRF(X(obj.idx.states(:,mod(iNode-1, nNodesDur-1) +1 )));
       %sigmoid_left = 0.0005 + 1 / (1 + exp(-50 * lfy_current+10 ));
        %sigmoid_right = 0.0005 + 1 / (1 + exp(-50 * rfy_current+10));
        
       
        ic = (1:nconstraintspernode) +  (iNode-1)*nconstraintspernode;
        delayed_index = mod(iNode - delay-2, nNodesDur-1) +1; %starts at 95 if we have 100 nodes or 96 if we have 101 nodes, gets NEGATIVE derivative
        delayed_index2 = mod(iNode - delay+1, nNodesDur-1) +1; %starts at 96 (for 100 nodes), is gets POSITIVE derivative
        next_index =  mod(iNode, nNodesDur - 1) + 1;

        fx2 = grfx(delayed_index2); %TODO: check if this is the same as when I add the 2 states
        fy2 = grfy(delayed_index2);
       
        fx1 = grfx(delayed_index); %TODO: check if this is the same as when I add the 2 states
        fy1 = grfy(delayed_index);

    
        v_curr = X(obj.idx.belt_left(iNode));
        v_prev = X(obj.idx.belt_left(mod(iNode - 2, nNodesDur-1)+1)); %only works if n_constraints (and not n+1) --> if I need n+1 points: use an if statement
      
        fy_current= grfy(mod(iNode-1, nNodesDur-1) +1);
        sigmoid_left = 0.0005 + 1 / (1 + exp(-50 * fy_current+10));


        %derivative wrt to next speed
        
        output(ic(1), idxVLeft(next_index)) = -1;

        output(ic(1), idxVLeft(iNode)) = sigmoid_left*(1 -Kp - Kd/c); % der from pd Part wrt v(n)
        output(ic(1), idxVLeft(mod(iNode - 2, nNodesDur-1)+1)) =  sigmoid_left*Kd/c; %derivative wrt v(n-1)


        %derivative wrt Kgrf
        output(ic(1), idxKfx) = sigmoid_left*((fx2 - fx1)/(3*c));
        
        %derivative wrt Kfx
        output(ic(1), idxKfy) = sigmoid_left*((fy2 - fy1)/(3*c));

        %derivative wrt Kp, Kd and Kpd
        output(ic(1), idxKp) = sigmoid_left*(1.8- v_curr);
        output(ic(1), idxKd) = sigmoid_left* ((-v_curr+ v_prev)/c);

        %derivative wrt c
       % output(ic(1), idxC) = -1 * Kgrf *((fx2 - fx1)/c^2) - Kgrf*Kfy*((fy2 - fy1)/c^2) - Kpd*Kd * ((-v_curr+ v_prev)/c^2);


     
        
    end
else
    error('Unknown option.');
end
end
