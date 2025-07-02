load("Winter_normal.mat");
exp = load("treadmillWalking18.mat"); %OR OTHER FILE

tread = dataStruct;
tread.variables.mean{7} = 1.8; %speed;

if tread.variables.mean{7} == 1.2
    tread.variables.mean{6} = 1.11;
    tread.variables.var{6} = 0.014;
else
    tread.variables.name{6} = 0.95;
    tread.variables.var{6} = 0.009;
end


tread.studyName = 'BeReal';
tread.participantHeight = 1.9;
tread.participantMass = 100;
tread.movementType = 'Walking';
tread.participantName = 'S01';
tread.movementDescription = 'Treadmill Walking';


tread.variables.mean{1} = deg2rad(exp.mean.rhip');
tread.variables.mean{2} = deg2rad(-exp.mean.rknee');
tread.variables.mean{3} = deg2rad(exp.mean.rank');

tread.variables.mean{4} = exp.mean.fx2'/100/9.81;
tread.variables.mean{5} = exp.mean.fy2'/100/9.81;

tread.variables.mean{8} = deg2rad(exp.mean.lhip');
tread.variables.mean{9} = deg2rad(-exp.mean.lknee');
tread.variables.mean{10} = deg2rad(exp.mean.lankle');

tread.variables.mean{11} = exp.mean.fx'/100/9.81;
tread.variables.mean{12} = exp.mean.fy'/100/9.81;


tread.variables.var{1} = deg2rad(exp.std.rhip');
tread.variables.var{2} = deg2rad(-exp.std.rknee');
tread.variables.var{3} = deg2rad(exp.std.rank');

tread.variables.var{4} = exp.std.fx2'/100/9.81;
tread.variables.var{5} = exp.std.fy2'/100/9.81;

tread.variables.var{8} = deg2rad(exp.std.lhip');
tread.variables.var{9} = deg2rad(-exp.std.lknee');
tread.variables.var{10} =deg2rad(exp.std.lankle');

tread.variables.var{11} = exp.std.fx'/100/9.81;
tread.variables.var{12} = exp.std.fy'/100/9.81;

tread.movementEvents.index(3) = 101;
tread.movementEvents.index(2) = 51;

dataStruct = tread;

save("treadmill_fast.mat", "dataStruct");


