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
function output = treadmillSpeedConstraintsSigmoid(obj,option,X,sym)

%% check input parameter
if ~isfield(obj.idx,'states') || ~isfield(obj.idx,'belt_left') || ~isfield(obj.idx,'belt_right')% check whether controls are stored in X
    error('Model states and left and right belt speeds need to be stored in state vector X.')
end


%% compute demanded output
nNodesDur = obj.nNodesDur; %number of collocation nodes + 1

nconstraintspernode = 2; %one for left and one for right

m = -obj.model.bodymass * obj.model.gravity(2);
delay = obj.model.grf_delay;
Kfy =  obj.model.Kfy;
Kgrf = obj.model.Kgrf;
Kp = obj.model.Kp;
Kd = obj.model.Kd;
Kpd = obj.model.Kpd;
c = obj.model.c;

if strcmp(option,'confun')
    output = zeros(nconstraintspernode*(nNodesDur-1),1);
    
    % dynamic equations must be zero
    for iNode=1:(nNodesDur-1)

        delayed_index = mod(iNode - delay-1, nNodesDur-1) +1;
        delayed_index2 = mod(iNode - delay, nNodesDur-1) +1;
        delayedIx2 = X(obj.idx.states(:,delayed_index2 )); %mroe recent then delayed Idx
        delayedIx = X(obj.idx.states(:,delayed_index));
        ic = (1:nconstraintspernode) +  (iNode-1)*nconstraintspernode; %indices of constraints of iNode in c
        
        grfDelayedIx = obj.model.getGRF(delayedIx);
        grfDelayedIx2 = obj.model.getGRF(delayedIx2);

        rfx2 = grfDelayedIx2(1)*m; %TODO: check if this is the same as when I add the 2 states
        rfy2 = grfDelayedIx2(2)*m;
        lfx2 = grfDelayedIx2(7)*m;
        lfy2 = grfDelayedIx2(8)*m;

        rfx1 = grfDelayedIx(1)*m;
        rfy1 = grfDelayedIx(2)*m;
        lfx1 = grfDelayedIx(7)*m;
        lfy1 = grfDelayedIx(8)*m;

        f_current= obj.model.getGRF(X(obj.idx.states(:,mod(iNode-1, nNodesDur-1) +1 )));
        rfy_current = f_current(2)*m;
        lfy_current= f_current(8)*m;

      
        v_left_curr = X(obj.idx.belt_left(iNode));
        v_left_prev = X(obj.idx.belt_left(mod(iNode - 2, nNodesDur-1)+1)); %only works if n_constraints (and not n+1) --> if I need n+1 points: use an if statement
        v_right_curr = X(obj.idx.belt_right(iNode));
        v_right_prev = X(obj.idx.belt_right(mod(iNode - 2, nNodesDur-1)+1));
        

        v_left = v_left_curr + Kgrf *((lfx2 - lfx1)/c) + Kgrf*Kfy*((lfy2 - lfy1)/c) + Kpd*Kp*(obj.model.speed_left - v_left_curr) + Kpd*Kd * ((-v_left_curr+ v_left_prev)/c);
        v_right = v_right_curr + Kgrf *((rfx2 - rfx1)/c) + Kgrf*Kfy*((rfy2 - rfy1)/c) + Kpd*Kp*(obj.model.speed_right - v_right_curr) + Kpd*Kd * ((-v_right_curr+ v_right_prev)/c);

        sigmoid_left = 0.0000001 + 1 / (1 + exp(-50 * lfy_current+1000));
        sigmoid_right = 0.0000001 + 1 / (1 + exp(-50 * rfy_current+1000));
        
        v_left_next = X(obj.idx.belt_left(mod(iNode, nNodesDur - 1) + 1));
        v_right_next = X(obj.idx.belt_right(mod(iNode, nNodesDur - 1) + 1));

        %apply sigmoid: belt speed should be 1.2 if there is no vertical
        %force, else the computed speed
        v_left = sigmoid_left*v_left + (1-sigmoid_left) * obj.model.speed_left;
        v_right = sigmoid_right*v_right + (1-sigmoid_right) * obj.model.speed_right;

        
        diff =  v_left- v_left_next ;
        diff2 = v_right - v_right_next ;
        output(ic) = [diff;diff2];	% backward Euler discretization
        
        end
    

elseif strcmp(option,'jacobian')
    output = spalloc(nconstraintspernode*(nNodesDur-1),length(X),obj.Jnnz);

    idxFxHeelR = obj.model.extractState('Fx', 'heel_r'); %the index in the state vector
    idxFxHeelRinX = obj.idx.states(idxFxHeelR, :); %all indices of Fx in X (should be len 1xn_nodes)

    idxFxToeR = obj.model.extractState('Fx', 'front_r'); %the index in the state vector
    idxFxToeRinX = obj.idx.states(idxFxToeR, :); %all indices of Fx in X (should be len 1xn_nodes)

    idxFxHeelL = obj.model.extractState('Fx','heel_l'); %the index in the state vector
    idxFxHeelLinX = obj.idx.states(idxFxHeelL, :); %all indices of Fx in X (should be len 1xn_nodes)

    idxFxToeL = obj.model.extractState('Fx', 'front_l'); %the index in the state vector
    idxFxToeLinX = obj.idx.states(idxFxToeL, :); %all indices of Fx in X (should be len 1xn_nodes)


    idxFyHeelR = obj.model.extractState('Fy', 'heel_r'); %the index in the state vector
    idxFyHeelRinX = obj.idx.states(idxFyHeelR, :); %all indices of Fx in X (should be len 1xn_nodes)

    idxFyToeR = obj.model.extractState('Fy', 'front_r'); %the index in the state vector
    idxFyToeRinX = obj.idx.states(idxFyToeR, :); %all indices of Fx in X (should be len 1xn_nodes)

    idxFyHeelL = obj.model.extractState('Fy','heel_l'); %the index in the state vector
    idxFyHeelLinX = obj.idx.states(idxFyHeelL, :); %all indices of Fx in X (should be len 1xn_nodes)

    idxFyToeL = obj.model.extractState('Fy', 'front_l'); %the index in the state vector
    idxFyToeLinX = obj.idx.states(idxFyToeL, :); %all indices of Fx in X (should be len 1xn_nodes)


    idxVLeft = obj.idx.belt_left; % the indices of the left Belt speed in X
    idxVRight = obj.idx.belt_right; % the indices of the right belt speed in X

    
    
    for iNode = 1:(nNodesDur-1)

        f_current= obj.model.getGRF(X(obj.idx.states(:,mod(iNode-1, nNodesDur-1) +1 )));
        rfy_current = f_current(2)*m;
        lfy_current= f_current(8)*m;

        sigmoid_left = 0.0000001 + 1 / (1 + exp(-50 * lfy_current+1000 ));
        sigmoid_right = 0.0000001 + 1 / (1 + exp(-50 * rfy_current+1000));
        
       
        ic = (1:nconstraintspernode) +  (iNode-1)*nconstraintspernode;
        delayed_index = mod(iNode - delay-1, nNodesDur-1) +1; %starts at 95 if we have 100 nodes or 96 if we have 101 nodes, gets NEGATIVE derivative
        delayed_index2 = mod(iNode - delay, nNodesDur-1) +1; %starts at 96 (for 100 nodes), is gets POSITIVE derivative
        next_index =  mod(iNode, nNodesDur - 1) + 1;

        delayedIx2 = X(obj.idx.states(:,delayed_index2 )); %mroe recent then delayed Idx
        delayedIx = X(obj.idx.states(:,delayed_index));
        
        grfDelayedIx = obj.model.getGRF(delayedIx);
        grfDelayedIx2 = obj.model.getGRF(delayedIx2);

        rfx2 = grfDelayedIx2(1)*m; %TODO: check if this is the same as when I add the 2 states
        rfy2 = grfDelayedIx2(2)*m;
        lfx2 = grfDelayedIx2(7)*m;
        lfy2 = grfDelayedIx2(8)*m;

        rfx1 = grfDelayedIx(1)*m;
        rfy1 = grfDelayedIx(2)*m;
        lfx1 = grfDelayedIx(7)*m;
        lfy1 = grfDelayedIx(8)*m;

        v_left_curr = X(obj.idx.belt_left(iNode));
        v_left_prev = X(obj.idx.belt_left(mod(iNode - 2, nNodesDur-1)+1)); %only works if n_constraints (and not n+1) --> if I need n+1 points: use an if statement
        v_right_curr = X(obj.idx.belt_right(iNode));
        v_right_prev = X(obj.idx.belt_right(mod(iNode - 2, nNodesDur-1)+1));
        

        v_left = v_left_curr + Kgrf *((lfx2 - lfx1)/c) + Kgrf*Kfy*((lfy2 - lfy1)/c) + Kpd*Kp*(obj.model.speed_left - v_left_curr) + Kpd*Kd * ((-v_left_curr+ v_left_prev)/c);
        v_right = v_right_curr + Kgrf *((rfx2 - rfx1)/c) + Kgrf*Kfy*((rfy2 - rfy1)/c) + Kpd*Kp*(obj.model.speed_right - v_right_curr) + Kpd*Kd * ((-v_right_curr+ v_right_prev)/c);



        %derivative of left belt wrt left GRFx (at states heel and toe)
        output(ic(1), idxFxToeLinX(delayed_index2)) = sigmoid_left*m*Kgrf/c;
        output(ic(1), idxFxHeelLinX(delayed_index2)) = sigmoid_left*m*Kgrf / c;
        output(ic(1), idxFxToeLinX(delayed_index)) = -sigmoid_left*m*Kgrf/c;
        output(ic(1), idxFxHeelLinX(delayed_index)) = -sigmoid_left*m*Kgrf / c;


        %derivative of left belt wrt left GRFy (at states heel and toe)
        output(ic(1), idxFyToeLinX(delayed_index2)) = sigmoid_left*Kfy * m*Kgrf/c;
        output(ic(1), idxFyHeelLinX(delayed_index2)) = sigmoid_left*Kfy * m*Kgrf / c;
        output(ic(1), idxFyToeLinX(delayed_index)) = -m*sigmoid_left*Kgrf* Kfy/c;
        output(ic(1), idxFyHeelLinX(delayed_index)) = -m*sigmoid_left*Kgrf*Kfy /c;

        %derivative wrt to next speed
        
        output(ic(1), idxVLeft(next_index)) = -1;

        %derivative wrt to current vertical force
        output(ic(1), idxFyHeelLinX(mod(iNode-1, nNodesDur-1) +1)) = sigmoid_left * (1-sigmoid_left) * v_left + (-obj.model.speed_left) * sigmoid_left * (1-sigmoid_left);
        output(ic(1), idxFyToeLinX(mod(iNode-1, nNodesDur-1) +1)) = sigmoid_left * (1-sigmoid_left) * v_left + (-obj.model.speed_left) * sigmoid_left * (1-sigmoid_left);


        %derivative of right belt wrt right GRFx (at states heel and toe)
        output(ic(2), idxFxToeRinX(delayed_index2)) = sigmoid_right*m*Kgrf/c;
        output(ic(2), idxFxHeelRinX(delayed_index2)) = sigmoid_right*m*Kgrf / c;
        output(ic(2), idxFxToeRinX(delayed_index)) = -sigmoid_right*m*Kgrf/c;
        output(ic(2), idxFxHeelRinX(delayed_index)) = -sigmoid_right*m*Kgrf / c;

        %derivative of right belt wrt right GRFy (at states heel and toe)
        output(ic(2), idxFyToeRinX(delayed_index2)) = Kfy * sigmoid_right*m*Kgrf/c;
        output(ic(2), idxFyHeelRinX(delayed_index2)) = Kfy *sigmoid_right* m*Kgrf / c;
        output(ic(2), idxFyToeRinX(delayed_index)) = -m*sigmoid_right*Kgrf* Kfy/c;
        output(ic(2), idxFyHeelRinX(delayed_index)) = -m*sigmoid_right*Kgrf*Kfy /c;



        output(ic(1), idxVLeft(iNode)) = sigmoid_left*(1 -Kpd*Kp - Kpd*Kd/c); % der from pd Part wrt v(n)
        output(ic(1), idxVLeft(mod(iNode - 2, nNodesDur-1)+1)) =  sigmoid_left*Kpd*Kd/c; %derivative wrt v(n-1)

        output(ic(2), idxVRight(iNode)) = sigmoid_right*(1 -Kpd*Kp - Kpd*Kd/c); % der from pd Part wrt v(n)
        output(ic(2), idxVRight(mod(iNode - 2, nNodesDur-1)+1)) =  sigmoid_right*Kpd*Kd/c; %derivative wrt v(n-1)

        %derivative of left belt wrt previous speed
        output(ic(2), idxVRight(next_index)) = -1;

        %derivative wrt to vertical force
        output(ic(2), idxFyHeelRinX(mod(iNode -1, nNodesDur-1) +1)) = sigmoid_right * (1-sigmoid_right) * v_right + (-obj.model.speed_right) * sigmoid_right * (1-sigmoid_right);
        output(ic(2), idxFyToeRinX(mod(iNode -1, nNodesDur-1) +1)) = sigmoid_right * (1-sigmoid_right) * v_right + (-obj.model.speed_right) * sigmoid_right * (1-sigmoid_right);
    end
else
    error('Unknown option.');
end
end
