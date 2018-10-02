function Learn_Random(filename)
% Learn_Random(filename) uses random search technique to find the optimal
% probability weighting function: OptimalThetaVec. 
% OptimalThetaVec is important for other tools in this toolbox, so make
% sure to run this function first and record the OptimalThetaVec. 
% Input: 
%       filename refers to the ID of the dataset.
% For example: 
%       Learn_Random('Sample')

tic;

% =========================================================================
% Random search for the best probability weighting function
% =========================================================================
% Define the total points to be searched.
SearchCap = 200000;

% =========================================================================
% Load data 
% =========================================================================

% Add path to the supplementary programs
addpath('SupplementaryProgram/');

% Refer to the specific data set.
%filename = 'Sample'; 
    
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

     

% Define and initialize a 5-D vector as weighted probabilities for original
% probability vector [0, 0.25, 0.5, 0.75, 1].
OrigProbVec = [0, 0.25, 0.5, 0.75, 1];
ThetaVec = -1*ones(1,5);
    
% Define RT performance matrix, the first 5 column is ThetaVec, while the
% is the amount of correct predictions
RT_Performance = -1*ones(SearchCap, 6);

% Search  for the optimal parameters
for i = 1:SearchCap
    
    % Randomly select a ThetaVec
        % Uniformly randomly select the first element of ThetaVec.
        Theta000 = rand;
        ThetaVec(1) = Theta000;

        % Use the monotonicity of probaility weighting function
        % and uniformaly randomly select the second element of ThetaVec.
        Theta025 = rand * (1 - Theta000) + Theta000;
        ThetaVec(2) = Theta025;
    
        % Use the monotonicity of probaility weighting function
        % and uniformaly randomly select the third element of ThetaVec.
        Theta050 = rand * (1 - Theta025) + Theta025;
        ThetaVec(3) = Theta050;
        
        % Use the monotonicity of probaility weighting function
        % and uniformaly randomly select the fourth element of ThetaVec.
        Theta075 = rand * (1 - Theta050) + Theta050;
        ThetaVec(4) = Theta075;
        
        % Use the monotonicity of probaility weighting function
        % and uniformaly randomly select the fourth element of ThetaVec.
        Theta100 = rand * (1 - Theta075) + Theta075;
        ThetaVec(5) = Theta100;
    
    % Error proof
    if min(ThetaVec) < 0
        error('ThetaVec is not defined.')
    end
    
    % Calculate the coordinates of points on the Q-function    
    [ Q_x1, Q_y1 ] = GetQ_LearnP( DataQ, ThetaVec);

    % Interpolate the term Q ( C(0) - C(x_H)) as Q_0H. The
    % property Q(-delta) = -Q(delta) is used here. 
    Q_0H_Vec = -1 * interp1(Q_x1, Q_y1, c_H_Vec, 'linear', 'extrap');

    % Interpolate the term Q ( C(R) - C(x_H)) as Q_RH
    Q_RH_Vec = interp1(Q_x1, Q_y1, c_R_Vec - c_H_Vec, 'linear', 'extrap');

    % Get the decision weight from probability weight function
    Pi_Prob_Vec = interp1(OrigProbVec, ThetaVec, Prob_Vec);

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
    RT_Performance(i,:) =...
        [ThetaVec, AmountCorrect_RT];

   
end

% Find the optimal parameters
MaxAmount_RT = max(RT_Performance(:, end));
MaxAmountIndx = find(RT_Performance(:, end) == MaxAmount_RT);

ThetaVec_star = RT_Performance(MaxAmountIndx, 1:end-1);
Num_Theta_star = size(ThetaVec_star, 1);

% Display
disp('Maximum Amount of correct prediction by RTx:');
disp(MaxAmount_RT);
disp('Intermediate optimal parameters');    
disp(ThetaVec_star);    

% ==================================================================
% Among the optimal parameters find the set giving best Q-shape
% ==================================================================

Q_fittype = fittype('Q_x1+s_r*exp(s_l*Q_x1)',...
                            'independent','Q_x1', 'dependent', 'Q_y1');

Q_fitCell = cell(Num_Theta_star, 1);
Q_fit_rsquare = -1*ones(Num_Theta_star, 1);

warning('off', 'all');

for i = 1: Num_Theta_star
    
    [ Q_x1, Q_y1 ] = GetQ_LearnP( DataQ, ThetaVec_star(i,:));
    
    [Q_FunctionFit, GoF] = fit(Q_x1', Q_y1', Q_fittype);
    
    Q_fitCell{i} = Q_FunctionFit;
    Q_fit_rsquare(i) = GoF.rsquare;
    
end

% Find the optiimal of the optimal parameters
BestFit_rsquare = max(Q_fit_rsquare);
BestFitIndx_rsquare= find(Q_fit_rsquare == BestFit_rsquare);

ThetaVec_2stars = ThetaVec_star(BestFitIndx_rsquare, :);
Num_Theta_2stars = size(ThetaVec_2stars, 1);

% Display
disp('R-square of the best fit');
disp(BestFit_rsquare);
disp('Final Optimal parameters (OptimalThetaVec)');    
disp(ThetaVec_2stars);

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
end