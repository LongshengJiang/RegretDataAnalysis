function [ Q_x, Q_y ] = GetQ_LearnP( DataQ, ThetaVec)

% Sort the Qdata in order according to the second conlumn
SortedDataQ = sortrows(DataQ,  2);

% Separate the Sorted DataQ into to catagories.
SortedDataQ_A = SortedDataQ(1:4,:); % the severe outcome of Robot option is -100
SortedDataQ_B = SortedDataQ(5:8,:); % the severe outcome of Robot option is -90

% Preset the Q-value for C = 0.5 as -1. 
Q_05 = -1;

% Calcuate the Q-value for C = 0.4.
Prob = SortedDataQ_B(1, 3);
Q_04 = NextQ(Prob, Q_05, ThetaVec);

% Calcuate the Q-value for C = 0.6.
Prob = SortedDataQ_A(1, 3);
Q_06 = NextQ(Prob, Q_04, ThetaVec);

% Calcuate the Q-value for C = 0.3.
Prob = SortedDataQ_B(2, 3);
Q_03 = NextQ(Prob, Q_06, ThetaVec);

% Calcuate the Q-value for C = 0.7.
Prob = SortedDataQ_A(2, 3);
Q_07 = NextQ(Prob, Q_03, ThetaVec);

% Calcuate the Q-value for C = 0.2.
Prob = SortedDataQ_B(3, 3);
Q_02 = NextQ(Prob, Q_07, ThetaVec);

% Calcuate the Q-value for C = 0.8.
Prob = SortedDataQ_A(3, 3);
Q_08 = NextQ(Prob, Q_02, ThetaVec);

% Calcuate the Q-value for C = 0.1.
Prob = SortedDataQ_B(4, 3);
Q_01 = NextQ(Prob, Q_08, ThetaVec);

% Calcuate the Q-value for C = 0.9.
Prob = SortedDataQ_A(4, 3);
Q_09 = NextQ(Prob, Q_01, ThetaVec);

% Construct the coordinates of points on Q-function.
Q_y = [Q_09, Q_08, Q_07, Q_06, Q_05, Q_04, Q_03, Q_02, Q_01, 0];
Q_x = -1* [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0];   
end

function [Q_Next] = NextQ(Prob, Q_Previous, ThetaVec)
% NextQ is used to calculate the next Q-coordinate giving the previous one.

% Introplate the weighted probability.
Pi_Prob = interp1([0, 0.25, 0.5, 0.75, 1], ThetaVec, Prob);

% Calculate the Q-coordinate.
Q_Next = Pi_Prob/(1-Pi_Prob) * Q_Previous;
end
