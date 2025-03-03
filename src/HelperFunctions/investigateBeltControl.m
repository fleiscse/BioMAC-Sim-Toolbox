function [diffs, errors, pds, grfs]= investigateBeltControl(resultWalking)

bL = resultWalking.problem.idx.belt_left;
compSpeed = resultWalking.X(bL);
bR = resultWalking.problem.idx.belt_right;
compSpeedR = resultWalking.X(bR);
plot(compSpeed)

grfx = resultWalking.X


model = resultWalking.problem.model;
m = -model.bodymass * model.gravity(2);
delay = model.grf_delay;
Kfy =  model.Kfy;
Kgrf = model.Kgrf;
Kp = model.Kp;
Kd = model.Kd;
Kpd = model.Kpd;
c = model.c;



data = resultWalking.problem.extractData(resultWalking.X);
row = strcmp(data.name, 'GRF_x_r');
grfx = data.sim{row};
row = strcmp(data.name, 'GRF_y_r');
grfy = data.sim{row};
nNodesDur=101;
diffs = zeros(1*(nNodesDur-1),1);
errors = zeros(1*(nNodesDur-1),1);
pds = zeros(1*(nNodesDur-1),1);
grfs = zeros(1*(nNodesDur-1),1);
compspeeds = zeros(1*(nNodesDur-1),1);
fxchanges = zeros(1*(nNodesDur-1),1);
fychanges = zeros(1*(nNodesDur-1),1);

    
    % dynamic equations must be zero
    for iNode=1:(nNodesDur-1)
        X = resultWalking.X;

        delayed_index = mod(iNode - delay-1, nNodesDur-1) +1;
        delayed_index2 = mod(iNode - delay, nNodesDur-1) +1;

        delayed_index = mod(iNode - delay-1, nNodesDur-1) +1;
        delayed_index2 = mod(iNode - delay, nNodesDur-1) +1;
        delayedIx2 = X(resultWalking.problem.idx.states(:,delayed_index2 )); %mroe recent then delayed Idx
        delayedIx = X(resultWalking.problem.idx.states(:,delayed_index));
        
     

        rfx2 = grfx(delayed_index2)*m; %TODO: check if this is the same as when I add the 2 states
        rfy2 = grfy(delayed_index2)*m;
        

        rfx1 = grfx(delayed_index)*m;
        rfy1 = grfy(delayed_index2)*m;

        grfDelayedIx = resultWalking.problem.model.getGRF(delayedIx);
        grfDelayedIx2 = resultWalking.problem.model.getGRF(delayedIx2);

        rfx2 = grfDelayedIx2(1)*m; %TODO: check if this is the same as when I add the 2 states
        rfy2 = grfDelayedIx2(2)*m;
        lfx2 = grfDelayedIx2(7)*m;
        lfy2 = grfDelayedIx2(8)*m;

        rfx1 = grfDelayedIx(1)*m;
        rfy1 = grfDelayedIx(2)*m;
        lfx1 = grfDelayedIx(7)*m;
        lfy1 = grfDelayedIx(8)*m;
      
        v_right_curr = compSpeedR(iNode);
        v_right_prev = compSpeedR(mod(iNode - 2, nNodesDur-1)+1);
        


 
        v_right = v_right_curr + Kgrf *((rfx2 - rfx1)/c) + Kgrf*Kfy*((rfy2 - rfy1)/c) + Kpd*Kp*(model.speed_right - v_right_curr) + Kpd*Kd * ((-v_right_curr+ v_right_prev)/c);
        fx_change=(rfx2 - rfx1)/c;
        fy_change=(rfy2 - rfy1)/c;
        grf_part = Kgrf *((rfx2 - rfx1)/c) + Kgrf*Kfy*((rfy2 - rfy1)/c);
        pd_part = Kpd*Kp*(model.speed_right - v_right_curr) + Kpd*Kd * ((-v_right_curr+ v_right_prev)/c);
        e = model.speed_right - v_right_curr;
        
        v_right_next = compSpeedR(mod(iNode, nNodesDur - 1) + 1);
        
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
