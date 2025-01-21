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
function output = treadmillSpeedConstraints(obj,option,X,sym)

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
      
        v_left_curr = X(obj.idx.belt_left(iNode));
        v_left_prev = X(obj.idx.belt_left(mod(iNode - 2, nNodesDur-1)+1)); %only works if n_constraints (and not n+1) --> if I need n+1 points: use an if statement
        v_right_curr = X(obj.idx.belt_right(iNode));
        v_right_prev = X(obj.idx.belt_right(mod(iNode - 2, nNodesDur-1)+1));
        

        v_left = v_left_curr + Kgrf *((lfx2 - lfx1)/c) + Kgrf*Kfy*((lfy2 - lfy1)/c) + Kpd*Kp*(obj.model.speed_left - v_left_curr) + Kpd*Kd * ((-v_left_curr+ v_left_prev)/c);
        v_right = v_right_curr + Kgrf *((rfx2 - rfx1)/c) + Kgrf*Kfy*((rfy2 - rfy1)/c) + Kpd*Kp*(obj.model.speed_right - v_right_curr) + Kpd*Kd * ((-v_right_curr+ v_right_prev)/c);
        
        v_left_next = X(obj.idx.belt_left(mod(iNode, nNodesDur - 1) + 1));
        v_right_next = X(obj.idx.belt_right(mod(iNode, nNodesDur - 1) + 1));
        
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
       
        ic = (1:nconstraintspernode) +  (iNode-1)*nconstraintspernode;
        delayed_index = mod(iNode - delay-1, nNodesDur-1) +1; %starts at 95 if we have 100 nodes or 96 if we have 101 nodes, gets NEGATIVE derivative
        delayed_index2 = mod(iNode - delay, nNodesDur-1) +1; %starts at 96 (for 100 nodes), is gets POSITIVE derivative
        next_index =  mod(iNode, nNodesDur - 1) + 1;
        %derivative of left belt wrt left GRFx (at states heel and toe)
        output(ic(1), idxFxToeLinX(delayed_index2)) = m*Kgrf/c;
        output(ic(1), idxFxHeelLinX(delayed_index2)) = m*Kgrf / c;
        output(ic(1), idxFxToeLinX(delayed_index)) = -m*Kgrf/c;
        output(ic(1), idxFxHeelLinX(delayed_index)) = -m*Kgrf / c;


        %derivative of left belt wrt left GRFy (at states heel and toe)
        output(ic(1), idxFyToeLinX(delayed_index2)) = Kfy * m*Kgrf/c;
        output(ic(1), idxFyHeelLinX(delayed_index2)) = Kfy * m*Kgrf / c;
        output(ic(1), idxFyToeLinX(delayed_index)) = -m*Kgrf* Kfy/c;
        output(ic(1), idxFyHeelLinX(delayed_index)) = -m*Kgrf*Kfy /c;

        %derivative wrt to next speed
        
        output(ic(1), idxVLeft(next_index)) = -1;


        %derivative of right belt wrt right GRFx (at states heel and toe)
        output(ic(2), idxFxToeRinX(delayed_index2)) = m*Kgrf/c;
        output(ic(2), idxFxHeelRinX(delayed_index2)) = m*Kgrf / c;
        output(ic(2), idxFxToeRinX(delayed_index)) = -m*Kgrf/c;
        output(ic(2), idxFxHeelRinX(delayed_index)) = -m*Kgrf / c;

        %derivative of right belt wrt right GRFy (at states heel and toe)
        output(ic(2), idxFyToeRinX(delayed_index2)) = Kfy * m*Kgrf/c;
        output(ic(2), idxFyHeelRinX(delayed_index2)) = Kfy * m*Kgrf / c;
        output(ic(2), idxFyToeRinX(delayed_index)) = -m*Kgrf* Kfy/c;
        output(ic(2), idxFyHeelRinX(delayed_index)) = -m*Kgrf*Kfy /c;

%         if delayed_index2==1 %then it is also a function of the extra node
%             output(ic(1), idxFxToeLinX(delayed_index2+nNodesDur-1)) = -Kgrf/c;
%             output(ic(1), idxFxHeelLinX(delayed_index2+nNodesDur-1)) = -Kgrf / c;
%             output(ic(1), idxFyToeLinX(delayed_index2+nNodesDur-1)) = -Kfy * Kgrf/c;
%             output(ic(1), idxFyHeelLinX(delayed_index2+nNodesDur-1)) = -Kfy * Kgrf / c;
% 
%             
%             
%         end
%         if delayed_index==1
%             output(ic(1), idxFxToeLinX(delayed_index+nNodesDur-1)) = -Kgrf/c;
%             output(ic(1), idxFxHeelLinX(delayed_index+nNodesDur-1)) = -Kgrf / c;
%             output(ic(1), idxFyToeLinX(delayed_index+nNodesDur-1)) = -Kgrf* Kfy/c;
%             output(ic(1), idxFyHeelLinX(delayed_index+nNodesDur-1)) = -Kgrf*Kfy /c;
%            
%         end




       

        output(ic(1), idxVLeft(iNode)) = 1 -Kpd*Kp - Kpd*Kd/c; % der from pd Part wrt v(n)
        output(ic(1), idxVLeft(mod(iNode - 2, nNodesDur-1)+1)) =  Kpd*Kd/c; %derivative wrt v(n-1)

        output(ic(2), idxVRight(iNode)) = 1 -Kpd*Kp - Kpd*Kd/c; % der from pd Part wrt v(n)
        output(ic(2), idxVRight(mod(iNode - 2, nNodesDur-1)+1)) =  Kpd*Kd/c; %derivative wrt v(n-1)

        %derivative of left belt wrt previous speed
        output(ic(2), idxVRight(next_index)) = -1;
        
    end
else
    error('Unknown option.');
end
end
