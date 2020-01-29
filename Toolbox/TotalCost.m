clc;
clear;
close all;

% Refer to the specific data set
filename = 'YG'; 

% Create member name of the data family and load the file
TrialNum = ['_' num2str(1) ];
load(['../Database/' filename TrialNum '.mat']);
load(['../Database/' filename TrialNum 'OutcomeScaleMeasure.mat']);

Sum = 0;
for i = 1:12
    for j = 1:10
        Sum = Sum + ResultArray{i,j}.MoneyCost;
    end
end

Sum_norm = Sum/MoneyScale;

disp(Sum);
disp(Sum_norm);