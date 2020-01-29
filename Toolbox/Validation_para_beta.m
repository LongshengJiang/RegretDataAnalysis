
clear;
close all;

% Add path to the supplementary programs
addpath('SupplementaryProgram/');
                    
% Refer to the specific data set.
filename = 'YG'; 
% Put the optimal parameters of the probability function here
beta_vec = [ 1.0611    1.1662    0.1263    0.7585];

        
% Load the data file
TrialNum = ['_' num2str(1)];
load(['../Database/' filename TrialNum '.mat']);
load(['../Database/' filename TrialNum 'OutcomeScaleMeasure.mat']);

% Extract the information for later use
[QuestionSetupArray, FinalChoiceArray, DegreeTruthArray] = GetV(ResultArray);

QuestionMat = QuestionSetupArray{1};

% Calculate the net advantage for each validation question using Regret
% theory
   
    % Define parameter
    slope_C = 1;
    % Calculate the term C(X_H_star) as c_H_Vec
    c_H_Vec = slope_C * QuestionMat(2,:)/MoneyScale;
    
    % Calculate the term C(X_R_star) as c_R_Vec
    c_R_Vec = slope_C * QuestionMat(3,:)/MoneyScale;

% Linear interpolation form
    % The requried question config is stored in column 11 in the result
    % array structure.
    ColumnNum = 11;
    
    % Load sample points
        DataQ = [ResultArray{3, ColumnNum}.Question(2:4);
         ResultArray{4,ColumnNum}.Question(2:4);
         ResultArray{6,ColumnNum}.Question(2:4);
         ResultArray{7,ColumnNum}.Question(2:4);
         ResultArray{8,ColumnNum}.Question(2:4);
         ResultArray{9,ColumnNum}.Question(2:4);
         ResultArray{10,ColumnNum}.Question(2:4);
         ResultArray{11,ColumnNum}.Question(2:4)];

% =========================================================================
%                       Prediction
% =========================================================================
     
% -------------------------------------------------------------------------
%                           RT
% -------------------------------------------------------------------------  
    
    % Get the x,y coordinates of the points for RT
    [ Q_x1_RT, Q_y1_RT ] = GetQ( DataQ );
    
    % Calculate the term Q ( C(0) - C(x_H_star)) as Q_0H
    Q_0H_Vec_RT = -1 * interp1(Q_x1_RT, Q_y1_RT, c_H_Vec, 'linear', 'extrap');
    % Calculate the term Q ( C(x_R) - C(x_H_star)) as Q_0H
    Q_RH_Vec_RT = interp1(Q_x1_RT, Q_y1_RT, c_R_Vec - c_H_Vec,  'linear', 'extrap');
    
    % Retrieve the term p as Prob_Vec 
    Prob_Vec = QuestionMat(4,:);
    
    % Calculate the net advantage e_Q for all questions as e_Q_Vec_RTx
    %
    % The net advantage of Robot to Human is caculated by the following
    % equation.
    %     ______________________________________________________
    %    |  e_Q =  p * Q ( C(0) - C(x_H_star))                  |
    %    |         +(1-p)* Q(C(x_R_star) - C(x_H_star))         |
    %    |______________________________________________________|
    %
    
    e_Q_Vec_RT = Prob_Vec .* Q_0H_Vec_RT ...
              + ( 1 - Prob_Vec ) .* Q_RH_Vec_RT;
    
        % Form the prediction by crispy regret theory
    ChoicePredict_RT = cell(1,10);
    for i = 1: 10
        
        % Predict the choice by following the logic of regret theory
        if e_Q_Vec_RT(i) > 0
            ChoicePredict_RT{i} = 'R';
        elseif e_Q_Vec_RT(i) < 0
            ChoicePredict_RT{i} = 'H';
        else
            ChoicePredict_RT{i} = 'I';
        end
           
    end
          
          
% -------------------------------------------------------------------------
%                           RTx_para
% -------------------------------------------------------------------------
    % Obtain the parametric w-function
    % Define the parametric model:
    %   w = beta3 + beta4 * exp( -beta1 * (-log(P_vec)) ^ beta2 )
    w_FunctionFit = @(p, beta) ...
                        beta(3) + beta(4) .* exp( -beta(1) .* ( -log(p) ).^ beta(2) );
   
    % Plot the w-function
    plot(linspace(0,1,30), w_FunctionFit(linspace(0,1,30), beta_vec));
    hold on;
    plot([0,1], [0,1], '--k');
    
    % Get the x,y coordinates of the points for RTx
    % We use the method used in the Learn_Random script to compute the
    % points on the q function.
    [ Q_x1_RTx, Q_y1_RTx ] = GetQ_LearnP_para( DataQ, beta_vec);
    
        % Fit the parametric q function
        % Define the fit type
        Q_fittype = fittype('s_t * Q_x1 + s_r * sinh( s_l*Q_x1 )',...
                            'independent','Q_x1', 'dependent', 'Q_y1');
        % Curve fit the q function
        [Q_FunctionFit, GoF_Q] = fit(Q_x1_RTx', Q_y1_RTx', Q_fittype,...
                                     'StartPoint', [1, 1, 1],...
                                     'Lower', [0, 0, 0],...
                                     'Upper', [10, 20, 500]);
   % Plot the q-function
   figure;
   plot(Q_FunctionFit, Q_x1_RTx, Q_y1_RTx);
       
    % Calculate the term Q ( C(0) - C(x_H_star)) as Q_0H
    Q_0H_Vec_RTx_column =  Q_FunctionFit( -c_H_Vec );
    Q_0H_Vec_RTx = Q_0H_Vec_RTx_column'; 
    % Calculate the term Q ( C(x_R) - C(x_H_star)) as Q_0H
    Q_RH_Vec_RTx_column = Q_FunctionFit( c_R_Vec - c_H_Vec );
    Q_RH_Vec_RTx = Q_RH_Vec_RTx_column';
    % Compute the probability weights for RTx
    Pi_Prob_Vec_RTx = w_FunctionFit( Prob_Vec, beta_vec);
    
    % Calculate the net advantage e_Q for all questions as e_Q_Vec_RTx
    %
    % The net advantage of Robot to Human is caculated by the following
    % equation.
    %     ______________________________________________________
    %    |  e_Q =  p * Q ( C(0) - C(x_H_star))                  |
    %    |         +(1-p)* Q(C(x_R_star) - C(x_H_star))         |
    %    |______________________________________________________|
    %
    
    e_Q_Vec_RTx = Pi_Prob_Vec_RTx .* Q_0H_Vec_RTx ...
              + ( 1 - Pi_Prob_Vec_RTx ) .* Q_RH_Vec_RTx;
          
          
    % Form the prediction by crispy regret theory with extension
    ChoicePredict_RTx = cell(1,10);
    for i = 1: 10
        
        % Predict the choices
        if e_Q_Vec_RTx(i) > 0
            ChoicePredict_RTx{i} = 'R';
        elseif e_Q_Vec_RTx(i) < 0
            ChoicePredict_RTx{i} = 'H';
        else
            ChoicePredict_RTx{i} = 'I';
        end
           
    end
          
    
% -------------------------------------------------------------------------
%                           EV
% -------------------------------------------------------------------------
 
% Calculate the net advantage by using the expected value criterion
    % Extract the moneny value x_R as x_R_Vec 
    x_R_Vec = QuestionMat(3,:)/MoneyScale;
    
    % Extract the money value x_H as x_H_Vec
    x_H_Vec = QuestionMat(2,:)/MoneyScale;
    
    % Calculate the net advantage e_EV for all questions as e_EV_Vec
    e_EV_Vec = (1-Prob_Vec) .* x_R_Vec - x_H_Vec;
    
    % Form the prediction by crispy expected value criterion
    ChoicePredict_EV = cell(1,10);
    for i = 1: 10
        
        % Predict the choice by following the logic of regret theory
        if e_EV_Vec(i) > 0
            ChoicePredict_EV{i} = 'R';
        elseif e_EV_Vec(i) < 0
            ChoicePredict_EV{i} = 'H';
        else
            ChoicePredict_EV{i} = 'I';
        end
           
    end

    
% -------------------------------------------------------------------------
%                           Human subjects
% -------------------------------------------------------------------------  

% Extract the actual choices from the degree of truth
    % Define an array to store the extracted final choices
   FinalChoiceArray = cell(1,2);
    
    % For each of the two sets of validation data, we determine the choices
    for i = 1:2
        % Define the temperary choice storage.
        FinalChoice = cell(1,10);
        DegreeTruth = DegreeTruthArray{i};
            for j = 1:10
                % If the degree of truth to prefer robot is positive
                if DegreeTruth(1, j) > 0
                    % the choice is preferring robot.
                    FinalChoice{j} = 'R';

                % Or if the degree of truth to prefer human is positive
                elseif DegreeTruth(2, j) > 0
                     % the choice is preferring human.
                    FinalChoice{j} = 'H';
                % Otherwise, the choice is equally like
                else
                    FinalChoice{j} = 'I';
                end
            end
        FinalChoiceArray{i} = FinalChoice;
    end    
      
% Extract the actual choices
FinalChoice1st = FinalChoiceArray{1};
FinalChoice2nd = FinalChoiceArray{2};
    

% =========================================================================
%                           Validation
% ========================================================================= 

% -------------------------------------------------------------------------
%                           RT
% -------------------------------------------------------------------------
% Calculate the predicting accuracy of RTx
    % For the first round choice of the user
    NumAccuratePred_RT_1st = 0;
    for i = 1:10
         if strcmp(ChoicePredict_RT{i}, FinalChoice1st{i})
             NumAccuratePred_RT_1st = NumAccuratePred_RT_1st + 1; 
         end
    end
    
    % For the second round choice of the user
    NumAccuratePred_RT_2nd = 0;
    for i = 1:10
         if strcmp(ChoicePredict_RT{i}, FinalChoice2nd{i})
             NumAccuratePred_RT_2nd = NumAccuratePred_RT_2nd + 1; 
         end
    end
    
    % For both
    NumAccuratePred_RT_Both = 0;
    for i = 1:10
         if strcmp(ChoicePredict_RT{i}, FinalChoice2nd{i}) && ...
                 strcmp(ChoicePredict_RT{i}, FinalChoice1st{i})
             NumAccuratePred_RT_Both = NumAccuratePred_RT_Both + 1; 
         end
    end


% -------------------------------------------------------------------------
%                           RTx_para
% -------------------------------------------------------------------------

% Calculate the predicting accuracy of RTx
    % For the first round choice of the user
    NumAccuratePred_RTx_1st = 0;
    for i = 1:10
         if strcmp(ChoicePredict_RTx{i}, FinalChoice1st{i})
             NumAccuratePred_RTx_1st = NumAccuratePred_RTx_1st + 1; 
         end
    end
    
    % For the second round choice of the user
    NumAccuratePred_RTx_2nd = 0;
    for i = 1:10
         if strcmp(ChoicePredict_RTx{i}, FinalChoice2nd{i})
             NumAccuratePred_RTx_2nd = NumAccuratePred_RTx_2nd + 1; 
         end
    end
    
    % For both
    NumAccuratePred_RTx_Both = 0;
    for i = 1:10
         if strcmp(ChoicePredict_RTx{i}, FinalChoice2nd{i}) && ...
                 strcmp(ChoicePredict_RTx{i}, FinalChoice1st{i})
             NumAccuratePred_RTx_Both = NumAccuratePred_RTx_Both + 1; 
         end
    end
    

% -------------------------------------------------------------------------
%                           EV
% -------------------------------------------------------------------------    

% Calculate the predicting accuracy of EV
    % For the first round choice of the user
    NumAccuratePred_EV_1st = 0;
    for i = 1:10
         if strcmp(ChoicePredict_EV{i}, FinalChoice1st{i})
             NumAccuratePred_EV_1st = NumAccuratePred_EV_1st + 1; 
         end
    end
    
    % For the second round choice of the user
    NumAccuratePred_EV_2nd = 0;
    for i = 1:10
         if strcmp(ChoicePredict_EV{i}, FinalChoice2nd{i})
             NumAccuratePred_EV_2nd = NumAccuratePred_EV_2nd + 1; 
         end
    end    
    
    % For both
    NumAccuratePred_EV_Both = 0;
    for i = 1:10
         if strcmp(ChoicePredict_EV{i}, FinalChoice2nd{i}) && ...
                 strcmp(ChoicePredict_EV{i}, FinalChoice1st{i})
             NumAccuratePred_EV_Both = NumAccuratePred_EV_Both + 1; 
         end
    end
    

% -------------------------------------------------------------------------
%                           the subject herself
% ------------------------------------------------------------------------- 

% Calculate the consistent answers provided by the participant
    % For the second round choice of the user
    NumConsistent = 0;
    for i = 1:10
         if strcmp(FinalChoice1st{i}, FinalChoice2nd{i})
             NumConsistent = NumConsistent + 1; 
         end
    end 
    
% =========================================================================
%                           Display
% =========================================================================
% Display
Newline = char(10);

disp(['The prediction from RT', Newline]);
disp(ChoicePredict_RT);

disp(['The prediction from RTx', Newline]);
disp(ChoicePredict_RTx);

disp(['The prediction from EV', Newline]);
disp(ChoicePredict_EV);

disp(['The first round choice of the user.', Newline]);
disp(FinalChoice1st);

disp(['The second round choice of the user.', Newline]);
disp(FinalChoice2nd);

disp(['The number of accurate prediction for the 1st data by RT.', Newline]);
disp(NumAccuratePred_RT_1st);

disp(['The number of accurate prediction for the 1st data by RTx.', Newline]);
disp(NumAccuratePred_RTx_1st);

disp(['The number of accurate prediction for the 1st data by EV.', Newline]);
disp(NumAccuratePred_EV_1st);

disp(['The number of accurate prediction for the 2nd data by RT.', Newline]);
disp(NumAccuratePred_RT_2nd);

disp(['The number of accurate prediction for the 2nd data by RTx.', Newline]);
disp(NumAccuratePred_RTx_2nd);

disp(['The number of accurate prediction for the 2nd data by EV.', Newline]);
disp(NumAccuratePred_EV_2nd);

disp(['The number of consistent answers from the participant', Newline]);
disp(NumConsistent);

disp(['The number of accurate prediction for the consistent data by RT.', Newline]);
disp(NumAccuratePred_RT_Both);

disp(['The number of accurate prediction for the consistent data by RTx.', Newline]);
disp(NumAccuratePred_RTx_Both);

disp(['The number of accurate prediction for the consistent data by EV.', Newline]);
disp(NumAccuratePred_EV_Both);

% Remove path to the supplementary programs
rmpath('SupplementaryProgram/');