function [diffs, errors, pds, grfs]= investigateBeltControl(resultWalking)

bL = resultWalking.problem.idx.belt_left;
compSpeed = resultWalking.X(bL);
%bR = resultWalking.problem.idx.belt_right;
%compSpeedR = resultWalking.X(bR);
plot(compSpeed)


%model = resultWalking.problem.model;
%m = -model.bodymass * model.gravity(2);
delay = 5%model.grf_delay;
Kfy =  resultWalking.X(101);
Kgrf = resultWalking.X(102);
Kp = resultWalking.X(103);
Kd = resultWalking.X(104);
Kpd = resultWalking.X(105);
c = 0.01% model.c;




grfx = resultWalking.problem.constraintTerms.varargin{1,1};
grfy = resultWalking.problem.constraintTerms.varargin{1,2};
nNodesDur=101;
diffs = zeros(1*(nNodesDur-1),1);
errors = zeros(1*(nNodesDur-1),1);
pds = zeros(1*(nNodesDur-1),1);
grfs = zeros(1*(nNodesDur-1),1);
compspeeds = zeros(1*(nNodesDur-1),1);
fxchanges = zeros(1*(nNodesDur-1),1);
fychanges = zeros(1*(nNodesDur-1),1);
v_right_curr = compSpeed(100);
v_right_prev = compSpeed(99);
    
    % dynamic equations must be zero
    for iNode=1:(nNodesDur-1)

        delayed_index = mod(iNode - delay-1, nNodesDur-1) +1;
        delayed_index2 = mod(iNode - delay, nNodesDur-1) +1;
        
     

        rfx2 = grfx(delayed_index2); %TODO: check if this is the same as when I add the 2 states
        rfy2 = grfy(delayed_index2);
        

        rfx1 = grfx(delayed_index);
        rfy1 = grfy(delayed_index);
      
      %  v_right_curr = compSpeed(iNode);
        %v_right_prev = compSpeed(mod(iNode - 2, nNodesDur-1)+1);
       % v_right_curr = v_right
      %  v_right_prev = v_right

        v_right = v_right_curr + Kgrf *((rfx2 - rfx1)/c) + Kgrf*Kfy*((rfy2 - rfy1)/c) + Kpd*Kp*(1.2 - v_right_curr) + Kpd*Kd * ((-v_right_curr+ v_right_prev)/c);
        fx_change=(rfx2 - rfx1)/c
        fy_change=(rfy2 - rfy1)/c;
        grf_part = Kgrf *((rfx2 - rfx1)/c) + Kgrf*Kfy*((rfy2 - rfy1)/c);
        pd_part = Kpd*Kp*(1.2 - v_right_curr) + Kpd*Kd * ((-v_right_curr+ v_right_prev)/c)
        e = 1.2 - v_right_curr;
        
        v_right_next = compSpeed(mod(iNode, nNodesDur - 1) + 1);
        v_right_prev = v_right_curr;
        v_right_curr = v_right;
        
        compspeeds(iNode)=v_right;
        diff = v_right - v_right_next ;
        diffs(iNode) = diff;	% backward Euler discretization
        pds(iNode) = pd_part;
        grfs(iNode) = grf_part;
        errors(iNode) = e;
        fxchanges(iNode) = fx_change;
        fychanges(iNode) = fy_change;
    end

      % what would beltspeed be if we just "kept going?

    v_right_curr = compSpeedR(2);
    v_right_prev = compSpeedR(1);
    allspeeds = zeros(1*(nNodesDur-1),1);
    allspeeds(1) = v_right_prev
    allspeeds(2) = v_right_curr
    for iNode=3:(nNodesDur-1)

        delayed_index = mod(iNode - delay-1, nNodesDur-1) +1;
        delayed_index2 = mod(iNode - delay, nNodesDur-1) +1;
        
     

        rfx2 = grfx(delayed_index2)*m; %TODO: check if this is the same as when I add the 2 states
        rfy2 = grfy(delayed_index2)*m;
        

        rfx1 = grfx(delayed_index)*m;
        rfy1 = grfy(delayed_index2)*m;
      
        
        


 
        v_right = v_right_curr + Kgrf *((rfx2 - rfx1)/c) + Kgrf*Kfy*((rfy2 - rfy1)/c) + Kpd*Kp*(model.speed_right - v_right_curr) + Kpd*Kd * ((-v_right_curr+ v_right_prev)/c);
        allspeeds(iNode)=v_right;
        v_right_prev = v_right_curr;
        v_right_curr = v_right;
        
        
       
    end


        
    

    




end
