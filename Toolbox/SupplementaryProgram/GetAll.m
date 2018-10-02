function [QuestionSetupMatAll, FinalChoiceAll, DegreeTruthMatAll] = GetAll(ResultArray)


% Define the empty matrix for storing the degree of truth of all questions,
% except questions of type V
DegreeTruthMatAll = -1*ones(5,100);

% Define the cell array to store the final choices by the user of all questions,
% except questions of type V
FinalChoiceAll = cell(1, 100);

% Define the empty matrix for storing the question setup of  of all questions,
% except questions of type V
QuestionSetupMatAll = -1*ones(4,100);

% Create the array of questions which contain results, excluding questions
% of type V
ResultArray_QW = [ResultArray(1:4, 1:10);
                  ResultArray(6:11, 1:10)];

% For each row of ResultArray_QW
for i = 1:10
    
    % For each column of ResultArray_QW
    for j = 1:10
    
    % Calculate the linear index constructed by i and j    
    Indx = (i-1)*10 + j;
    % extracting info
    DegreeTruthMatAll(:,Indx) = ResultArray_QW{i,j}.DegreeOfTruth;
    QuestionSetupMatAll(:,Indx) = ResultArray_QW{i, j}.Question';
                                                                            %FinalChoiceAll{Indx} = ResultArray_QW{i, j}.PreviousChoice;
     % If the degree of truth to prefer robot is positive
    if ResultArray_QW{i,j}.DegreeOfTruth(1) > 0
        % the choice is preferring robot.
        FinalChoiceAll{Indx} = 'R';

    % Or if the degree of truth to prefer human is positive
    elseif ResultArray_QW{i,j}.DegreeOfTruth(2) > 0
         % the choice is preferring human.
        FinalChoiceAll{Indx} = 'H';
    % Otherwise, the choice is equally like
    else
        FinalChoiceAll{Indx} = 'I';
    end
    
    end
end

end