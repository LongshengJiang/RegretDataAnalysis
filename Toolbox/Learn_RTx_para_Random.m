clear;
close all;

tic;

% =========================================================================
% Load data 
% =========================================================================

% Add path to the supplementary programs
addpath('SupplementaryProgram/');

% Refer to the specific data set.
filename = 'YG'; 
    
% Load the data file
TrialNum = ['_' num2str(1)];
load(['../Database/' filename TrialNum '.mat']);
load(['../Database/' filename TrialNum 'OutcomeScaleMeasure.mat']);

% Extract the information for later use
[QuestionMat, FinalChoice, DegreeTruth] = GetAll(ResultArray);

% Calculate the choiceless utilities

    % Calculate the term C(X_H_star) as c_H_Vec
    c_H_Vec =  QuestionMat(2,:)/MoneyScale;
    % Calculate the term C(X_R_star) as c_R_Vec
    c_R_Vec =  QuestionMat(3,:)/MoneyScale;
    % Retrieve the term p as Prob_Vec 
    Prob_Vec = QuestionMat(4,:);
    
% Retrieve the real choices by the subjects.
    % Get the degree of truth for Prefering Robot
    PrefRobot_Vec = DegreeTruth(1,:);
    PrefHuman_Vec = DegreeTruth(2,:);
    
    % Get the real choices
    ChooseRobot_Real = PrefRobot_Vec > 0;
    ChooseHuman_Real = PrefHuman_Vec > 0;
    EitherOK_Real = PrefRobot_Vec == 0 & PrefHuman_Vec == 0;
 

% Retrieve the question configuration given the indifferent condition 
    % The requried question config is stored in column 11 in the result
    % array structure.
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

     
% =========================================================================
% Random search for the best probability weighting function
% =========================================================================


% Define the total points to be searched.
SearchCap = 200000;
Loosen    = 0;

% Define RT performance matrix, the first 2 column is beta_vec, while the
% is the amount of correct predictions
RTx_Performance = -1*ones(SearchCap, 5);

% Search  for the optimal parameters
for i = 1:SearchCap
    % Randomly select a value for the parameter beta1
        % Uniformly randomly select a value for log_beta1 from 4*[-0.5, 0.5].
        log_beta1 = 4 * ( rand - 1/2 );
        beta1 = exp( log_beta1 );

        % Uniformly randomly select a value for log_beta2 from 4* [-0.5, 0.5].
        log_beta2 = 4 * ( rand - 1/2 );
        beta2 = exp( log_beta2 );
      
        % Uniformly randomly select a value for beta3 from [0, 0.5].
        beta3 = rand*0.5;
       
        % Uniformly randomly select a value for beta3 from [0.5, 1].
        beta4 = min( rand*0.5 +0.5, 1-beta3 );

        % Pack the beta parameters into a vector
        beta_vec = [beta1, beta2, beta3, beta4];
    
    % Calculate the coordinates of points on the Q-function    
    [ Q_x1, Q_y1 ] = GetQ_LearnP_para( DataQ, beta_vec);
    
    % Interpolate the term Q ( C(0) - C(x_H)) as Q_0H. The
    % property Q(-delta) = -Q(delta) is used here. 
    Q_0H_Vec = -1 * interp1(Q_x1, Q_y1, c_H_Vec, 'linear', 'extrap');

    % Interpolate the term Q ( C(R) - C(x_H)) as Q_RH
    Q_RH_Vec = interp1(Q_x1, Q_y1, c_R_Vec - c_H_Vec, 'linear', 'extrap');

%     % do regression on the indifference points (Q_x1, Q_y1)
%     % We use the parametric q function:
%     %           q = alpha1 * sinh(alpha2 * Q_x1) + alpha3 * Q_x1 
%     %
%     Q_fittype = fittype('alpha1 * sinh( alpha2 * Q_x1 ) + alpha3 * Q_x1 ',...
%                             'independent','Q_x1', 'dependent', 'Q_y1');
%     
%     Q_FunctionFit = fit(Q_x1', Q_y1', Q_fittype,...
%                                  'StartPoint', [1, 1, 1],...
%                                  'Lower', [0, 0, 0],...
%                                  'Upper', [10, 20, 500]);                    
%     
%     % Compute the term Q ( C(0) - C(x_H)) as Q_0H. 
%     Q_0H_Vec_col = Q_FunctionFit( -c_H_Vec );
%     Q_0H_Vec = Q_0H_Vec_col';
%     
%     % Interpolate the term Q ( C(R) - C(x_H)) as Q_RH
%     Q_RH_Vec_col = Q_FunctionFit( c_R_Vec - c_H_Vec);
%     Q_RH_Vec = Q_RH_Vec_col';
    
    
    % Get the decision weight from probability weight function
    % We use the parametric w-function:
    %           w = exp( -beta1 * ( -log( p ) )^beta2 )
    w_FunctionFit = @(p, beta) ...
                        beta(3) + beta(4) .* exp( -beta(1) .* ( -log(p) ).^ beta(2) );
    
    Pi_Prob_Vec = w_FunctionFit( Prob_Vec, beta_vec);

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
    RTx_Performance(i,:) =...
        [beta_vec, AmountCorrect_RT];

   
end

% Find the optimal parameters
MaxAmount_RTx = max(RTx_Performance(:, end))-Loosen;
MaxAmountIndx = find(RTx_Performance(:, end) >= MaxAmount_RTx);

beta_vec_star = RTx_Performance(MaxAmountIndx, 1:end-1);
Num_Theta_star = size(beta_vec_star, 1);

% Display
disp('Maximum Amount of correct prediction by RT:');
disp(MaxAmount_RTx);
disp('Optimal parameters');    
disp(beta_vec_star);    

% ==================================================================
% Among the optimal parameters find the set giving best Q-shape
% ==================================================================

Q_fittype = fittype('alpha1 * sinh( alpha2 * Q_x1 ) + alpha3 * Q_x1 ',...
                            'independent','Q_x1', 'dependent', 'Q_y1');

Q_fitCell = cell(Num_Theta_star, 1);
Q_fit_rsquare = -1*ones(Num_Theta_star, 1);

warning('off', 'all');

disp('Num_Theta_star:');
disp(Num_Theta_star);

% Wait the user to continue
waitforbuttonpress;


for i = 1: Num_Theta_star
    
    disp([num2str(i) '/' num2str(Num_Theta_star)]);
    
    [ Q_x1, Q_y1 ] = GetQ_LearnP_para( DataQ, beta_vec_star(i,:));
    
    [Q_FunctionFit, GoF] = fit(Q_x1', Q_y1', Q_fittype,...
                                 'StartPoint', [1, 1, 1],...
                                 'Lower', [0, 0, 0],...
                                 'Upper', [20, 10, 500]);
    
    Q_fitCell{i} = Q_FunctionFit;
    Q_fit_rsquare(i) = GoF.rsquare;
    
end

% Find the optiimal of the optimal parameters
BestFit_rsquare = max(Q_fit_rsquare);
BestFitIndx_rsquare= find(Q_fit_rsquare == BestFit_rsquare);

beta_vec_2stars = beta_vec_star(BestFitIndx_rsquare, :);
Num_Theta_2stars = size(beta_vec_2stars, 1);

% Plot the optimal w function
w_Plot_p = linspace(0, 1, 50);
plot(w_Plot_p, w_FunctionFit(w_Plot_p, beta_vec_2stars));
hold on;
plot([0,1], [0,1], '--k');
hold off;
% Plot the optimal q function
figure;
[ Q_x1, Q_y1 ] = GetQ_LearnP_para( DataQ, beta_vec_2stars);
plot(Q_fitCell{BestFitIndx_rsquare}, Q_x1, Q_y1);

 % Display
disp('Maximum Amount of correct prediction by RT:');
disp(MaxAmount_RTx);
disp(Q_fitCell{BestFitIndx_rsquare});
disp('R-square of the best fit');
disp(BestFit_rsquare);
disp('Final Optimal parameters');    
disp(beta_vec_2stars);

% =========================================================================
% The peroformance of EV
% =========================================================================
% Calculate the net advantage based on EV
x_R_EV_Vec = QuestionMat(3,:)./ MoneyScale;
x_H_EV_Vec = QuestionMat(2,:) ./ MoneyScale;

% Calculate the net advantage
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

% Display
disp('Amount of correc prediction by EV:');
disp(AmountCorrect_EV);
   

warning('on', 'all'); 

toc;