function TrainingPerformance(filename, OptimalThetaVec)
% TrainingPerformance(filename, OptimalThetaVec) shows the
% classification results on the training data. 
% Inputs:
%       filename refers to the ID code of the dataset.
%       OptimalThetaVec is the final optimal parameters learned in function
%       Learn_Random(filename)
% For example,
%       TrainingPerformance('Sample', [0.0146 0.3914 0.5096 0.6095 0.9967])
%       shows the classification performance on the training data of the 
%       sample dataset. 

% Add path to the supplementary programs
addpath('SupplementaryProgram/');

% Load the data file
TrialNum = ['_' num2str(1)];
load(['../Database/' filename TrialNum '.mat']);
load(['../Database/' filename TrialNum 'OutcomeScaleMeasure.mat']);

% Extract the information for later use
[QuestionMat, FinalChoice, DegreeTruth] = GetAll(ResultArray);

% Calculate the net advantage for each validation question using Regret
% theory 
    % The equation in the following box is the elicited C function
    %     ______________________________________________________
    %    |      C(x_star) = slope_C * x_star,                   |
    %    | where slope_C = 1                                    |
    %    |______________________________________________________|
    % where x_star is the scaled outcome which is defined by x_star = x/x_s, in
    % which x is the original money amout and x_s is the amount the user feels
    % significantly important.

    % Define parameters
    slope_C =1;
    % Calculate the term C(X_H_star) as c_H_Vec
    c_H_Vec = slope_C * QuestionMat(2,:)/MoneyScale;
    
    % Calculate the term C(X_R_star) as c_R_Vec
    c_R_Vec = slope_C * QuestionMat(3,:)/MoneyScale;
    
    % Get the degree of truth
    % Get the degree of truth for Prefering Robot
    PrefRobot_Vec = DegreeTruth(1,:);

    % Get the degree of truth for Prefering Robot
    PrefHuman_Vec = DegreeTruth(2,:);
    
 % Get the real choices from the participant's response
    ChooseRobot_Real = PrefRobot_Vec > 0;
    ChooseHuman_Real = PrefHuman_Vec > 0;
    EitherOK_Real = PrefRobot_Vec == 0 & PrefHuman_Vec == 0;
 

% ===========================================================
% Expected Value Criterion
% ===========================================================

% Calculate the net advantage based on EV
x_R_EV_Vec = QuestionMat(3,:)./ MoneyScale;
x_H_EV_Vec = QuestionMat(2,:) ./ MoneyScale;


% Calculate the term p as Prob_Vec 
Prob_Vec = QuestionMat(4,:);

e_EV_Vec = ( 1 - Prob_Vec ) .* x_R_EV_Vec - x_H_EV_Vec;

% Evaluate the prediction performance    
    % Get the predictions by the EV model
    ChooseRobot_EV = e_EV_Vec > 0;
    ChooseHuman_EV = e_EV_Vec < 0;
    EitherOK_EV    = e_EV_Vec == 0;

    % Calculate the amount of correct predictions for EV
        % Calculate the amount of correct predictions in choosing Robot
        AmountCorrectRobot_EV = ...
               sum( ChooseRobot_EV == ChooseRobot_Real & ChooseRobot_EV == 1 );

        % Calculate the amount of correct predictions in choosing Human
        AmountCorrectHuman_EV = ...
               sum( ChooseHuman_EV == ChooseHuman_Real & ChooseHuman_EV == 1 );  

        % Calculate the amount of correct predictions in choosing Either is okay.
        AmountCorrectEitherOK_EV = ...
               sum( EitherOK_EV == EitherOK_Real & EitherOK_EV == 1 );  

        % Calculate the total amount of correct predictions
        AmountCorrect_EV = AmountCorrectRobot_EV + AmountCorrectHuman_EV...
                           + AmountCorrectEitherOK_EV;    
    
% ===========================================================
% Regret Theory
% ===========================================================
    
% Linear interpolation form
% Linear interpolation
    % Load sample points
    % Get the x,y coordinates for each point 
    % Extract data from the database
    ColumnNum = 11;
    
    DataQ = [ResultArray{3, ColumnNum}.Question(2:4);
         ResultArray{4,ColumnNum}.Question(2:4);
         ResultArray{6,ColumnNum}.Question(2:4);
         ResultArray{7,ColumnNum}.Question(2:4);
         ResultArray{8,ColumnNum}.Question(2:4);
         ResultArray{9,ColumnNum}.Question(2:4);
         ResultArray{10,ColumnNum}.Question(2:4);
         ResultArray{11,ColumnNum}.Question(2:4)];
         
    % Sort the data in order
    SortedDataQ = sortrows(DataQ,  2);
    

% Define RT performance matrix, the 1st column is Theta1, the 2nd column is
% Theta2, the 3rd column is Theta3, and the 4th column is the amount of
% correct predictions
RT_Performance = -1*ones(1, 6);

% Evaluate the classification performance using the optimal parameters

[ Q_x1, Q_y1 ] = GetQ_LearnP( DataQ, OptimalThetaVec);

% Calculate the term Q ( C(0) - C(x_H_star)) as Q_0H
Q_0H_Vec = -1 * interp1(Q_x1, Q_y1, c_H_Vec, 'linear', 'extrap');

% Linear interpolation
Q_RH_Vec = interp1(Q_x1, Q_y1, c_R_Vec - c_H_Vec, 'linear', 'extrap');

% Get the decision weight from probability weight function
Pi_Prob_Vec = interp1([0, 0.25, 0.5, 0.75, 1], OptimalThetaVec, Prob_Vec);

% The net advantage of Robot to Human is caculated by the following
% equation.
%     ______________________________________________________
%    |  e_Q =  p * Q ( C(0) - C(x_H_star))                  |
%    |         +(1-p)* Q(C(x_R_star) - C(x_H_star))         |
%    |______________________________________________________|
%

% Calculate the net advantage e_Q for all questions as e_Q_Vec
e_Q_Vec = Pi_Prob_Vec .* Q_0H_Vec ...
          + ( 1 - Pi_Prob_Vec ) .* Q_RH_Vec;

% Evaluate the prediction performance of RT
% Get the predictions by the RT model
ChooseRobot_RT = e_Q_Vec > 0;
ChooseHuman_RT = e_Q_Vec < 0;
EitherOK_RT    = e_Q_Vec == 0;

% Calculate the amount of correct predictions for RT
% Calculate the amount of correct predictions in choosing Robot
AmountCorrectRobot_RT = ...
       sum( ChooseRobot_RT == ChooseRobot_Real & ChooseRobot_RT == 1 );

% Calculate the amount of correct predictions in choosing Human
AmountCorrectHuman_RT = ...
       sum( ChooseHuman_RT == ChooseHuman_Real & ChooseHuman_RT == 1 );  

% Calculate the amount of correct predictions in choosing Either is okay.
AmountCorrectEitherOK_RT = ...
       sum( EitherOK_RT == EitherOK_Real & EitherOK_RT == 1 );  

% Calculate the total amount of correct predictions
AmountCorrect_RT = AmountCorrectRobot_RT + AmountCorrectHuman_RT...
                   + AmountCorrectEitherOK_RT;

% Store the information in to the RT performance matrix
RT_Performance = [OptimalThetaVec, AmountCorrect_RT]; 

% Choose only the sets of parameters giving the separation performance
% above the threshold. 
ThetaVec_star = RT_Performance(1:end-1);
Num_Theta_star = size(ThetaVec_star, 1);

% Display
disp('Amount of correc prediction by EV:');
disp(AmountCorrect_EV);
disp('Amount of correct prediction by RTx:');
disp(AmountCorrect_RT);
disp('Optimal parameters');    
disp(ThetaVec_star);    

