clear;
close all;

% Add path to the supplementary programs
addpath('SupplementaryProgram');

% Refer to the specific data set
filename = 'YG'; 

% Put the optimal parameters of the probability function here
ThetaVec = [  0.0146    0.3914    0.5096    0.6095    0.9967];

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
ColumnNum = 11;

% Retrieve the question information. The information is related to
    % eliciting Q-function, thus is named DataQ. At this point
    DataQ = [ResultArray{3, ColumnNum}.Question(2:4);
         ResultArray{4,ColumnNum}.Question(2:4);
         ResultArray{6,ColumnNum}.Question(2:4);
         ResultArray{7,ColumnNum}.Question(2:4);
         ResultArray{8,ColumnNum}.Question(2:4);
         ResultArray{9,ColumnNum}.Question(2:4);
         ResultArray{10,ColumnNum}.Question(2:4);
         ResultArray{11,ColumnNum}.Question(2:4)];
% Get the x,y coordinates for each point
[ Q_x1, Q_y1 ] = GetQ_LearnP( DataQ, ThetaVec);

% Plot the figure for each data file
plot(Q_x1(1:end),Q_y1(1:end),'-o');
    
grid on;
ylabel('Q(\delta)')
xlabel('\delta')
title('Q function ');

% Curve fitting
% The equation in the following box is the elicited Q function
%     ______________________________________________________
%    |   Q(delta_c) = delta_c + s_r*sinh(s_l*delta_c)       |
%    |  where s_r = 0.9135, s_l = 4.732                     |
%    |______________________________________________________|
%
% use x as the horizontal variable while y as the vertical variable.
%     x = linspace(-0.9*slope_W, 0);
%     s_r = 0.9135;
%     s_l = 4.732;
%     y = x + s_r * sinh(s_l * x);
%     plot(x, y, 'b--');
%     
%     legend('Q-function data points', 'Q(\delta_c) = \delta_c + 0.9135 sinh(4.732 \delta_c)');
%     legend('location', 'southeast');
% hold off;

% Remove path to the supplementary programs
rmpath('SupplementaryProgram');