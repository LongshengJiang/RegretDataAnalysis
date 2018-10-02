function Plot_Q(filename, OptimalThetaVec)
% Plot_Q(filename, OptimalThetaVec) shows the curve of the Q-function in
% RTx.
% Inputs:
%       filename refers to the ID code of the dataset.
%       OptimalThetaVec is the final optimal parameters learned in function
%       Learn_Random(filename)
% For example,
%       Plot_Q('Sample', [0.0146 0.3914 0.5096 0.6095 0.9967])
%       shows the shape of the Q-function of the sample dataset. 


% Add path to the supplementary programs
addpath('SupplementaryProgram');

figure;
hold on;

% Create member name of the data family and load the file
TrialNum = ['_' num2str(1) ];
load(['../Database/' filename TrialNum '.mat']);
load(['../Database/' filename TrialNum 'OutcomeScaleMeasure.mat']);

% Assign values to parameters
%Alpha    = 1;
Q_scale  = -1;
slope_W  =  1;

% Get the x,y coordinates for each point
[ Q_x1, Q_y1 ] = GetQ( ResultArray,slope_W,Q_scale, 11, OptimalThetaVec);

% Display Q_y1
disp('Q_y1= ');
disp(Q_y1);

% Plot the figure for each data file
plot(Q_x1(1:end),Q_y1(1:end),'-o');
    
grid on;
ylabel('Q(\delta)')
xlabel('\delta')
title('Q function ');

% Remove path to the supplementary programs
rmpath('SupplementaryProgram');