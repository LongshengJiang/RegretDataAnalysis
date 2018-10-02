function ValidationPerformance(filename, OptimalThetaVec)
% ValidationPerformance(filename, OptimalThetaVec) shows the
% classification results on the validation data. 
% Inputs:
%       filename refers to the ID code of the dataset.
%       OptimalThetaVec is the final optimal parameters learned in function
%       Learn_Random(filename)
% For example,
%       ValidationPerformance('Sample', [0.0146 0.3914 0.5096 0.6095 0.9967])
%       shows the classification performance on the validation data of 
%       the sample dataset. 


% Add path to the supplementary programs
addpath('SupplementaryProgram/');

% Load the data file
TrialNum = ['_' num2str(1)];
load(['../Database/' filename TrialNum '.mat']);
load(['../Database/' filename TrialNum 'OutcomeScaleMeasure.mat']);

% Extract the information for later use
[QuestionSetupArray, FinalChoiceArray, DegreeTruthArray] = GetV(ResultArray);

QuestionMat = QuestionSetupArray{1};

% Calculate the net advantage for each validation question using Regret
% theory
    
    % 
    % The equation in the following box is the elicited C function
    %     ______________________________________________________
    %    |      C(x_star) = slope_C * x_star                    |
    %    |   where slope_C = 0.9628                             |
    %    |______________________________________________________|
    % where x_star is the scaled outcome which is defined by x_star = x/x_s, in
    % which x is the original money amout and x_s is the amount the user feels
    % significantly important.

    % Define parameter
    slope_C = 1;
    % Calculate the term C(X_H_star) as c_H_Vec
    c_H_Vec = slope_C * QuestionMat(2,:)/MoneyScale;
    
    % Calculate the term C(X_R_star) as c_R_Vec
    c_R_Vec = slope_C * QuestionMat(3,:)/MoneyScale;

    % Instruction:
    %
    % The equation in the following box is the fitted Q-function
    %     ______________________________________________________
    %    |   Q(delta_c) = delta_c + s_r*sinh(s_l*delta_c)       | 
    %    |   where s_r = 2.727, s_l = 2.909                     |
    %    |______________________________________________________|
    

% Linear interpolation form
    % Load sample points
    % Get the x,y coordinates for each point
    [ Q_x1, Q_y1 ] = GetQ( ResultArray,slope_C,-1,11, OptimalThetaVec);
    
    % Calculate the term Q ( C(0) - C(x_H_star)) as Q_0H
    Q_0H_Vec = -1 * interp1(Q_x1, Q_y1, c_H_Vec, 'linear', 'extrap');

    % Linear interpolation
    Q_RH_Vec = interp1(Q_x1, Q_y1, c_R_Vec - c_H_Vec,  'linear', 'extrap');

    % Calculate the term p as Prob_Vec 
    Prob_Vec = QuestionMat(4,:);
    
    Pi_Prob_Vec = interp1([0, 0.25, 0.5, 0.75 1], OptimalThetaVec, Prob_Vec);
    
    % Calculate the net advantage e_Q for all questions as e_Q_Vec
    %
    % The net advantage of Robot to Human is caculated by the following
    % equation.
    %     ______________________________________________________
    %    |  e_Q =  p * Q ( C(0) - C(x_H_star))                  |
    %    |         +(1-p)* Q(C(x_R_star) - C(x_H_star))         |
    %    |______________________________________________________|
    %
    
    e_Q_Vec = Pi_Prob_Vec .* Q_0H_Vec ...
              + ( 1 - Pi_Prob_Vec ) .* Q_RH_Vec;
          
          
    % Form the prediction by crispy regret theory
    ChoicePredict_RT = cell(1,10);
    for i = 1: 10
        
        % Predict the choice by following the logic of regret theory
        if e_Q_Vec(i) > 0
            ChoicePredict_RT{i} = 'R';
        elseif e_Q_Vec(i) < 0
            ChoicePredict_RT{i} = 'H';
        else
            ChoicePredict_RT{i} = 'I';
        end
           
    end
          
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
    
% Calculate the predicting accuracy of RT
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
    
% Calculate the consistent answers provided by the participant
    % For the second round choice of the user
    NumConsistent = 0;
    for i = 1:10
         if strcmp(FinalChoice1st{i}, FinalChoice2nd{i})
             NumConsistent = NumConsistent + 1; 
         end
    end 
    
    
% Display
Newline = char(10);
disp(['The prediction from Regret Theory with extension.', Newline]);
disp(ChoicePredict_RT);

disp(['The prediction from Expected Value Criterion.', Newline]);
disp(ChoicePredict_EV);

disp(['The first round choice of the user.', Newline]);
disp(FinalChoice1st);

disp(['The second round choice of the user.', Newline]);
disp(FinalChoice2nd);

disp(['The number of accurate prediction for the 1st data by RTx.', Newline]);
disp(NumAccuratePred_RT_1st);

disp(['The number of accurate prediction for the 1st data by EV.', Newline]);
disp(NumAccuratePred_EV_1st);

disp(['The number of accurate prediction for the 2nd data by RTx.', Newline]);
disp(NumAccuratePred_RT_2nd);

disp(['The number of accurate prediction for the 2nd data by EV.', Newline]);
disp(NumAccuratePred_EV_2nd);

disp(['The number of consistent answers from the participant', Newline]);
disp(NumConsistent);

disp(['The number of accurate prediction for the consistent data by RTx.', Newline]);
disp(NumAccuratePred_RT_Both);

disp(['The number of accurate prediction for the consistent data by EV.', Newline]);
disp(NumAccuratePred_EV_Both);

% Remove path to the supplementary programs
rmpath('SupplementaryProgram/');