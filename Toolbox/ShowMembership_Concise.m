clear;
clc;
close all;

% Add path to the supplementary programs
addpath('SupplementaryProgram/');

% Refer to the specific data set.
filename = 'YG'; 
% Put the optimal parameters of the probability function here
ThetaVec = [  0.0146    0.3914    0.5096    0.6095    0.9967];

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
    %    | where slope_C = 0.9628                               |
    %    |______________________________________________________|
    % where x_star is the scaled outcome which is defined by x_star = x/x_s, in
    % which x is the original money amout and x_s is the amount the user feels
    % significantly important.

    % Define parameters
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
% Linear interpolation
    % Load sample points
    % Get the x,y coordinates for each point
    [ Q_x1, Q_y1 ] = GetQ( ResultArray,slope_C,-1,11, ThetaVec);
    
    % Calculate the term Q ( C(0) - C(x_H_star)) as Q_0H
    Q_0H_Vec = -1 * interp1(Q_x1, Q_y1, c_H_Vec, 'linear', 'extrap');

% Linear interpolation
    Q_RH_Vec = interp1(Q_x1, Q_y1, c_R_Vec - c_H_Vec, 'linear', 'extrap');

    
    % Calculate the term p as Prob_Vec 
    Prob_Vec = QuestionMat(4,:);
    % Pi_Prob_Vec = Prob_Vec;
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

% Get the degree of truth
    % Get the degree of truth for Prefering Robot
    PrefRobot_Vec = DegreeTruth(1,:);

    % Get the degree of truth for Prefering Robot
    PrefHuman_Vec = DegreeTruth(2,:);

    % Get the degree of truth for Indifference
    Indif_Vec = DegreeTruth(3,:);
    
    % Contruct a matrix to store the membership info. This matrix is of
    % dimension 100 by 4. 
    MembershipMat = [e_Q_Vec', PrefRobot_Vec',PrefHuman_Vec', Indif_Vec'];
    % Sorted the matrix to make e_Q_Vec increasing
    MembershipMatSort = sortrows(MembershipMat);
    
    % Get individual truth values
    x = MembershipMatSort(:,1);
    Membership_R = MembershipMatSort(:,2);
    Membership_H = MembershipMatSort(:,3);
    Membership_I = MembershipMatSort(:,4);
 
    
% Find the points corresponding to R, H, I, respectively
    % Find the index of the data points where 'R' is chosen.
    % Construct the logic vector for R
    LogicVec_R = zeros(1,100);
    for i = 1:100
        if strcmpi(FinalChoice{i}, 'R')
            LogicVec_R(i) = 1;
        end
    end
    % Get the index of data point choosing 'R'
    Indx_R = find(LogicVec_R);

    % Find the index of the data points where 'H' is chosen.
    % Construct the logic vector for H
    LogicVec_H = zeros(1,100);
    for i = 1:100
        if strcmpi(FinalChoice{i}, 'H')
            LogicVec_H(i) = 1;
        end
    end
    % Get the index of data point choosing 'H'
    Indx_H = find(LogicVec_H);
    
     % Find the index of the data points where 'I' is chosen.
        % Construct the logic vector for I
        LogicVec_I = zeros(1,100);
        for i = 1:100
            if strcmpi(FinalChoice{i}, 'I')
                LogicVec_I(i) = 1;
            end
        end
        % Get the index of data point choosing 'H'
        Indx_I = find(LogicVec_I);

        
% =====================================================================     
% Plot the membership function of RT
% =====================================================================

% Plots of Prefering Robot
xrange = [-1.5, 1.5];

subplot(2,1,1);
 
% Plots of Indifferent

    % Plot all the points where the choice is 'R' 
    % Add +0.025 to offset the data points, avoiding obstructions.
    plot(e_Q_Vec(Indx_R), Indif_Vec(Indx_R)+0.025 , 'rs', 'Markerfacecolor', 'g');
    hold on;
    
    % Plot all the points where the choice is 'H'
    plot(e_Q_Vec(Indx_H), Indif_Vec(Indx_H)-0.025, 'bo', 'Markerfacecolor', 'g');
    
    % Plot all the points where the choice is 'I'
    plot(e_Q_Vec(Indx_I), Indif_Vec(Indx_I), 'g^', 'Markerfacecolor', 'g');
    
    % Plot the decision making line
    plot([0, 0], [-0.35, 1.35], '--k', 'LineWidth', 2);
    
    % Plot editing
    xlim(xrange);
    ylim([-0.2, 1.35]);  
    set(gca, 'ytick', [0, 0.25, 0.5, 0.75, 1]);
    grid on;
    
    title('Membership function of Indifference (RT)');
    xlabel('e_Q');
    ylabel('\mu_I');
    
    
    hold off; 
    

% =====================================================================     
% Plot the membership function of EV
% ===================================================================== 
  
% Calculate the net advantage based on EV
x_R_EV_Vec = QuestionMat(3,:)./ MoneyScale;
x_H_EV_Vec = QuestionMat(2,:) ./ MoneyScale;

e_EV_Vec = ( 1 - Prob_Vec ) .* x_R_EV_Vec - x_H_EV_Vec;

ScaleFactor = 6;
e_EV_Vec = ScaleFactor * e_EV_Vec;

% Define the degree of truth in EV, i.e. e_EV_Vec>0 --> PrefRobot_EV=1
%                                                       PrefHuman_EV=0
%                                                       Indif_EV    =0
%
%                                        e_EV_Vec<0 --> PrefRobot_EV=0
%                                                       PrefHuman_EV=1
%                                                       Indif_EV    =0
%
%                                        e_EV_Vec=0 --> PrefRobot_EV=0
%                                                       PrefHuman_EV=0
%                                                       Indif_EV    =1

    % Calculate PrefRobot_EV_Vec
    Indx_e_EV_Positive = find(e_EV_Vec > 0);
    PrefRobot_EV_Vec = zeros(100, 1);
    PrefRobot_EV_Vec(Indx_e_EV_Positive) = 1;

    % Calculate PrefHuman_EV_Vec
    Indx_e_EV_Negative = find(e_EV_Vec < 0);
    PrefHuman_EV_Vec = zeros(100,1);
    PrefHuman_EV_Vec(Indx_e_EV_Negative) = 1;

    % Calculate Indif_EV_Vec
    Indx_e_EV_zero = find(e_EV_Vec == 0);
    Indif_EV_Vec = zeros(100,1);
    Indif_EV_Vec(Indx_e_EV_zero) = 1;

% Plots of Prefering Robot
xrange = [-1.5, 1.5];

% Plots of Indifferent
subplot(2,1,2);

    % Plot the membership grades of indifference
    % Plot all the points where the choice is 'R'    
    plot(e_EV_Vec(Indx_R), Indif_EV_Vec(Indx_R)+0.05, 'rs', 'Markerfacecolor', 'g');
    hold on;
    
    % Plot all the points where the choice is 'H'
    plot(e_EV_Vec(Indx_H), Indif_EV_Vec(Indx_H)-0.05, 'bo', 'Markerfacecolor', 'g');

    % Plot all the points where the choice is 'I'
    plot(e_EV_Vec(Indx_I), Indif_EV_Vec(Indx_I), 'g^', 'Markerfacecolor', 'g');

    % Plot the decision making line
    plot([0, 0], [-0.35, 1.35], '--k', 'LineWidth', 2);
    
    % Plot editing
    xlim(xrange);
    ylim([-0.2, 1.35]);  
    set(gca, 'ytick', [0, 1]);
    grid on;

    title('Membership function of Indifference (EV)');
    xlabel('e_{EV}');
    ylabel('\mu_I');


    hold off; 

% Print the graph to Gallery
print(['../Gallery/ShowMemebership' filename TrialNum], '-dpng');
print(['../Gallery/ShowMemebership' filename TrialNum], '-depsc');

% Remove path to the supplementary programs
rmpath('SupplementaryProgram/');  