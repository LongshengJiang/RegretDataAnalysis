function [QuestionSetupArray, FinalChoiceArray, DegreeTruthArray] = GetV(ResultArray)


% Define the empty matrix for storing the degree of truth of the 1st and the 2nd 
% validation question set.
DegreeTruthMat1st = -1*ones(5,10);
DegreeTruthMat2nd = -1*ones(5,10);

% Define the cell array to store the final choices by the user for the 1st
% and the 2nd validation question set.
FinalChoice1st = cell(1, 10);
FinalChoice2nd = cell(1, 10);

% Define the empty matrix for storing the question setup of the 1st and the
% 2nd validation question set.
QuestionSetupMat1st = -1*ones(4,10);
QuestionSetupMat2nd = -1*ones(4,10);

for i = 1:10
    
    % extracting info form the 1st validation set
    DegreeTruthMat1st(:,i) = ResultArray{5,i}.DegreeOfTruth;
    FinalChoice1st{i} = ResultArray{5, i}.PreviousChoice;
    QuestionSetupMat1st(:,i) = ResultArray{5, i}.Question';

     % extracting info form the 2nd validation set
    DegreeTruthMat2nd(:,i) = ResultArray{12,i}.DegreeOfTruth;
    FinalChoice2nd{i} = ResultArray{12, i}.PreviousChoice;
    QuestionSetupMat2nd(:,i) = ResultArray{12, i}.Question';
end



DegreeTruthArray = {DegreeTruthMat1st, DegreeTruthMat2nd};
FinalChoiceArray = {FinalChoice1st, FinalChoice2nd};
QuestionSetupArray = {QuestionSetupMat1st, QuestionSetupMat2nd};


end