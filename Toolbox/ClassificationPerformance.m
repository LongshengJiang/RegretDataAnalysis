clear;
close all;

% Add path to the supplementary programs
addpath('SupplementaryProgram/');

% Refer to the specific data set.
filename = 'YG'; 
% Put the optimal parameters of the probability function here
ThetaVec = [0.1790    0.3822    0.5054    0.5361    0.8975];

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
% Regret Theory with extension non-parametric
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
    
% Evaluate the classification performance using the optimal parameters

[ Q_x1, Q_y1 ] = GetQ_LearnP( SortedDataQ, ThetaVec);

% Calculate the term Q ( C(0) - C(x_H_star)) as Q_0H
Q_0H_Vec = -1 * interp1(Q_x1, Q_y1, c_H_Vec, 'linear', 'extrap');

% Linear interpolation
Q_RH_Vec = interp1(Q_x1, Q_y1, c_R_Vec - c_H_Vec, 'linear', 'extrap');

% Get the decision weight from probability weight function
Pi_Prob_Vec = interp1([0, 0.25, 0.5, 0.75, 1], ThetaVec, Prob_Vec);

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
ChooseRobot_RTx = e_Q_Vec > 0;
ChooseHuman_RTx = e_Q_Vec < 0;
EitherOK_RTx    = e_Q_Vec == 0;

% Calculate the amount of correct predictions for RT
% Calculate the amount of correct predictions in choosing Robot
AmountCorrectRobot_RTx = ...
       sum( ChooseRobot_RTx == ChooseRobot_Real & ChooseRobot_RTx == 1 );

% Calculate the amount of correct predictions in choosing Human
AmountCorrectHuman_RTx = ...
       sum( ChooseHuman_RTx == ChooseHuman_Real & ChooseHuman_RTx == 1 );  

% Calculate the amount of correct predictions in choosing Either is okay.
AmountCorrectEitherOK_RTx = ...
       sum( EitherOK_RTx == EitherOK_Real & EitherOK_RTx == 1 );  

% Calculate the total amount of correct predictions
AmountCorrect_RTx = AmountCorrectRobot_RTx + AmountCorrectHuman_RTx...
                   + AmountCorrectEitherOK_RTx;


% ===========================================================
% Regret Theory with extension parametric
% ===========================================================

% Curve fit to obtain the parametric w-function
    % Specify the x-coordinates
    P_vec = 0:0.25:1;
    % Define the parametric model:
    %   w = beta3 + beta4 * exp( -beta1 * (-log(P_vec)) ^ beta2 )
    % Coefficent range:
    %   beta1 is in [0.05, 20]
    %   beta2 is in [0.05, 20]
    %   beta3 is in [0, 1]
    %   beta4 is in [0, 1]
    w_fittype = fittype('beta3 + beta4 * exp( -beta1 * (-log(P_vec)) ^ beta2 )',...
                            'independent','P_vec', 'dependent', 'ThetaVec');
    % Fit the w-function
    [w_FunctionFit, GoF_w] = fit(P_vec', ThetaVec', w_fittype, ...
                                 'StartPoint', [1, 1, 0, 1], ...
                                 'Lower', [0.05, 0.05, 0, 0], ...
                                 'Upper', [20, 20, 1, 1]);
    disp(w_FunctionFit);
    disp(GoF_w);
    % Plot the w-function
    plot(w_FunctionFit, P_vec, ThetaVec);
% Fit the parametric q function
% Define the fit type
Q_fittype = fittype('s_t * Q_x1 + s_r * sinh( s_l*Q_x1 )',...
                    'independent','Q_x1', 'dependent', 'Q_y1');
% Curve fit the q function
[Q_FunctionFit, GoF_Q] = fit(Q_x1', Q_y1', Q_fittype,...
                             'StartPoint', [1, 1, 1],...
                             'Lower', [0, 0, 0],...
                             'Upper', [10, 20, 500]);
disp(Q_FunctionFit);
disp(GoF_Q);
 % Plot the q-function
   figure;
   plot(Q_FunctionFit, Q_x1, Q_y1);
                                                  
% Calculate the term Q ( C(0) - C(x_H_star)) as Q_0H
Q_0H_Vec_para =  Q_FunctionFit( -c_H_Vec );
% Calculate the term Q ( C(x_R) - C(x_H_star)) as Q_0H
Q_RH_Vec_para = Q_FunctionFit( c_R_Vec - c_H_Vec );

% Compute the probability weights for RTx
Pi_Prob_Vec_para = w_FunctionFit( Prob_Vec );

% The net advantage of Robot to Human is caculated by the following
% equation.
%     ______________________________________________________
%    |  e_Q =  p * Q ( C(0) - C(x_H_star))                  |
%    |         +(1-p)* Q(C(x_R_star) - C(x_H_star))         |
%    |______________________________________________________|
%

% Calculate the net advantage e_Q for all questions as e_Q_Vec
e_Q_Vec_para = Pi_Prob_Vec_para .* Q_0H_Vec_para ...
          + ( 1 - Pi_Prob_Vec_para ) .* Q_RH_Vec_para;


% Evaluate the prediction performance of RT
% Get the predictions by the RT model
ChooseRobot_RTx_para = e_Q_Vec_para > 0;
ChooseHuman_RTx_para = e_Q_Vec_para < 0;
EitherOK_RTx_para    = e_Q_Vec_para == 0;

% Calculate the amount of correct predictions for RT
% Calculate the amount of correct predictions in choosing Robot
AmountCorrectRobot_RTx_para = ...
       sum( ChooseRobot_RTx_para' == ChooseRobot_Real & ChooseRobot_RTx_para' == 1 );

% Calculate the amount of correct predictions in choosing Human
AmountCorrectHuman_RTx_para = ...
       sum( ChooseHuman_RTx_para' == ChooseHuman_Real & ChooseHuman_RTx_para' == 1 );  

% Calculate the amount of correct predictions in choosing Either is okay.
AmountCorrectEitherOK_RTx_para = ...
       sum( EitherOK_RTx_para' == EitherOK_Real & EitherOK_RTx_para' == 1 );  

% Calculate the total amount of correct predictions
AmountCorrect_RTx_para = AmountCorrectRobot_RTx_para + AmountCorrectHuman_RTx_para...
                   + AmountCorrectEitherOK_RTx_para;

% Display
disp('Optimal parameters');    
disp(ThetaVec);    
disp('Amount of correc prediction by EV:');
disp(AmountCorrect_EV);
disp('Amount of correct prediction by RTx_nonpara:');
disp(AmountCorrect_RTx);
disp('Amount of correct prediction by RTx_para:');
disp(AmountCorrect_RTx_para);


