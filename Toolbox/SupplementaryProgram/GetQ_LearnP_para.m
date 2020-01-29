function [ Q_x, Q_y ] = GetQ_LearnP_para( DataQ, beta_vec)

% Sort the Qdata in order according to the second conlumn
SortedDataQ = sortrows(DataQ,  2);

% Separate the Sorted DataQ into to catagories.
SortedDataQ_A = SortedDataQ(1:4,:); % the severe outcome of Robot option is -100
SortedDataQ_B = SortedDataQ(5:8,:); % the severe outcome of Robot option is -90

% Preset the Q-value for C = 0.5 as -1. 
Q_05 = -1;

% Calcuate the Q-value for C = 0.4.
Prob = SortedDataQ_B(1, 3);
Q_04 = NextQ(Prob, Q_05, beta_vec);

% Calcuate the Q-value for C = 0.6.
Prob = SortedDataQ_A(1, 3);
Q_06 = NextQ(Prob, Q_04, beta_vec);

% Calcuate the Q-value for C = 0.3.
Prob = SortedDataQ_B(2, 3);
Q_03 = NextQ(Prob, Q_06, beta_vec);

% Calcuate the Q-value for C = 0.7.
Prob = SortedDataQ_A(2, 3);
Q_07 = NextQ(Prob, Q_03, beta_vec);

% Calcuate the Q-value for C = 0.2.
Prob = SortedDataQ_B(3, 3);
Q_02 = NextQ(Prob, Q_07, beta_vec);

% Calcuate the Q-value for C = 0.8.
Prob = SortedDataQ_A(3, 3);
Q_08 = NextQ(Prob, Q_02, beta_vec);

% Calcuate the Q-value for C = 0.1.
Prob = SortedDataQ_B(4, 3);
Q_01 = NextQ(Prob, Q_08, beta_vec);

% Calcuate the Q-value for C = 0.9.
Prob = SortedDataQ_A(4, 3);
Q_09 = NextQ(Prob, Q_01, beta_vec);

% Construct the coordinates of points on Q-function.
Q_y = [Q_09, Q_08, Q_07, Q_06, Q_05, Q_04, Q_03, Q_02, Q_01, 0];
Q_x = -1* [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0];   
end

function [Q_Next] = NextQ(Prob, Q_Previous, beta_vec)
% NextQ is used to calculate the next Q-coordinate giving the previous one.

% We use the parametric w-function:
%           w = beta3 + beta4 * exp( -beta1 * ( -log( p ) )^beta2 )
%
beta1 = beta_vec(1);
beta2 = beta_vec(2);
beta3 = beta_vec(3);
beta4 = beta_vec(4);
Pi_Prob =  min( beta3 + beta4* exp( -beta1 * ( -log(Prob) )^beta2 ), 0.999999999);

% Calculate the Q-coordinate.
Q_Next = Pi_Prob/(1-Pi_Prob) * Q_Previous;
end
