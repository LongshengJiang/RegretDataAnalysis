function OutcomeDataCleanTool(filename, RowIndx)
% OutcomeDataCleanTool(filename, RowIndx) visualizes the change of outcomes
% within module 1 and 2. 
% Inputs: 
%       filename refers to the ID code of the dataset,
%       RowIndx refers to the row of ResultArray which stores the regret
%       data. 
% For example: 
%       OutcomeDataCleanTool('Sample', 2) shows the plot of the
%       outcomes within the 2nd module. 


% Load the data file
TrialNum = ['_' num2str(1)];
load(['../Database/' filename TrialNum '.mat']);


Outcome = zeros(1,11);
Chosen = cell(1,11);
DegreeOfChooseR = -1*ones(1,11);
DegreeOfChooseH = -1*ones(1,11);
DegreeOfChooseI = -1*ones(1,11);


for i = 1:11
    Outcome(i) = ResultArray{RowIndx,i}.Question(2);
    Chosen{i} = ResultArray{RowIndx,i}.PreviousChoice;
    DegreeOfChooseR(i) = ResultArray{RowIndx,i}.DegreeOfTruth(1);
    DegreeOfChooseH(i) = ResultArray{RowIndx,i}.DegreeOfTruth(2);
    DegreeOfChooseI(i) = ResultArray{RowIndx,i}.DegreeOfTruth(3);
end
plot(Outcome, '--o')
shg
disp('Choice:');
disp(Chosen(1:10))
disp('Degree of Choosing R:');
disp(DegreeOfChooseR(1:10));

disp('Degree of Choosing H:');
disp(DegreeOfChooseH(1:10));

disp('Degree of Choosing I:');
disp(DegreeOfChooseI(1:10));

disp('Question Setup:');
Question = ResultArray{RowIndx,1}.Question(1:4);
disp(Question)