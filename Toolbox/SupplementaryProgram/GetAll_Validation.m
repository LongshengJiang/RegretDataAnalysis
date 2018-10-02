function [QuestionSetupMatAll, FinalChoiceAll, DegreeTruthMatAll] = GetAll_Validation(ResultArray, CaseID)


% Define the empty matrix for storing the degree of truth of all questions,
% except questions of type V
DegreeTruthMatAll = -1*ones(5,10);

% Define the cell array to store the final choices by the user of all questions,
% except questions of type V
FinalChoiceAll = cell(1, 10);

% Define the empty matrix for storing the question setup of  of all questions,
% except questions of type V
QuestionSetupMatAll = -1*ones(4,10);

% Create the array of validation questions
% of type V
switch CaseID
    % When the CaseID is 1, retrieve the first validation set.
    case 1
        ResultArray_Validation = ResultArray(5, 1:10);
    % When the CaseID is 2, retrieve the second validation set.    
    case 2
        ResultArray_Validation = ResultArray(12, 1:10);
end

% For each row of ResultArray_QW
for i = 1:1
    
    % For each column of ResultArray_QW
    for j = 1:10
    
    % Calculate the linear index constructed by i and j    
    Indx = (i-1)*10 + j;
    % extracting info
    DegreeTruthMatAll(:,Indx) = ResultArray_Validation{i,j}.DegreeOfTruth;
    FinalChoiceAll{Indx} = ResultArray_Validation{i, j}.PreviousChoice;
    QuestionSetupMatAll(:,Indx) = ResultArray_Validation{i, j}.Question';

    end
end

end