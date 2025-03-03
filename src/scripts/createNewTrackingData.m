data = load("data\Walking\Winter_normal.mat");

grfx = load("data\Walking\grfx.mat");
grfy = load("data\Walking\grfy.mat");

rightX = (grfx.speed_18 / 100 / 9.81).';
rightY = (grfy.speed_18 / 100 / 9.81).';

first_half = rightX(1:50);   % First 50 samples
second_half = rightX(51:100); % Last 50 samples
leftX = [second_half, first_half];
% 
first_half = rightY(1:50);   % First 50 samples
second_half = rightY(51:100); % Last 50 samples
leftY = [second_half, first_half];

data.dataStruct.variables.mean{4} = rightX;
data.dataStruct.variables.mean{5} = rightY;


data.dataStruct.variables.mean{11} = leftX;
data.dataStruct.variables.mean{12} = leftY;

save('data/Walking/WalkingFast.mat', "data")