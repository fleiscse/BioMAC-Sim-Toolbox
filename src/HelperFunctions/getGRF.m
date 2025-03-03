function [fx, fy] = getGRF(resultWalking)

m = - resultWalking.problem.model.bodymass * resultWalking.problem.model.gravity(2);


settings.GRF   = {'GRF_x_r', 'GRF_y_r'};
res = resultWalking.problem.extractData(resultWalking.X, settings, [], 1);
idxFx = find(strcmp(res.name, 'GRF_x_r') ==1);
idxFy = find(strcmp(res.name, 'GRF_y_r') ==1);
fx = res.sim{idxFx} * m;
fy = res.sim{idxFy} * m ;


end
