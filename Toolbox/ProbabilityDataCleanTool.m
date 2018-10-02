function ProbabilityDataCleanTool(filename, RowIndx)
% OutcomeDataCleanTool(filename, RowIndx) visualizes the change of
% probabilities within module 3 ~ 4 and 6 ~ 11. 
% Inputs: 
%       filename refers to the ID code of the dataset,
%       RowIndx refers to the row of ResultArray which stores the regret
%       data. 
% For example: 
%       ProbabilityDataCleanTool('Sample', 6) shows the plot of the
%       probabilities within the 6th module. 


% Load the data file
TrialNum = ['_' num2str(1)];
load(['../Database/' filename TrialNum '.mat']);


Probability = zeros(1,11);
Chosen = cell(1,11);
DegreeOfChooseR = -1*ones(1,11);
DegreeOfChooseH = -1*ones(1,11);
DegreeOfChooseI = -1*ones(1,11);


for i = 1:11
    Probability(i) = ResultArray{RowIndx,i}.Question(4);
    Chosen{i} = ResultArray{RowIndx,i}.PreviousChoice;
    DegreeOfChooseR(i) = ResultArray{RowIndx,i}.DegreeOfTruth(1);
    DegreeOfChooseH(i) = ResultArray{RowIndx,i}.DegreeOfTruth(2);
    DegreeOfChooseI(i) = ResultArray{RowIndx,i}.DegreeOfTruth(3);
end
plot(Probability, '--o')
shg;