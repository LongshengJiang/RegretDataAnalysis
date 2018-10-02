function [ Q_x, Q_y ] = GetQ( ResultArray , slope_W, Q_scale, ColumnNum, ThetaVec)


DataQ = [ResultArray{3, ColumnNum}.Question(2:4);
         ResultArray{4,ColumnNum}.Question(2:4);
         ResultArray{6,ColumnNum}.Question(2:4);
         ResultArray{7,ColumnNum}.Question(2:4);
         ResultArray{8,ColumnNum}.Question(2:4);
         ResultArray{9,ColumnNum}.Question(2:4);
         ResultArray{10,ColumnNum}.Question(2:4);
         ResultArray{11,ColumnNum}.Question(2:4)];

SortedDataQ = sortrows(DataQ,  2);

SortedDataQ_A = SortedDataQ(1:4,:);
SortedDataQ_B = SortedDataQ(5:8,:);

%Q_y        = Q_scale * SortedDataQ(:,3) ./(1-SortedDataQ(:,3));
Q_05 = Q_scale;

%
Prob = SortedDataQ_B(1, 3);
Q_04 = NextQ(Prob, Q_05, ThetaVec);

%
Prob = SortedDataQ_A(1, 3);
Q_06 = NextQ(Prob, Q_04, ThetaVec);

%
Prob = SortedDataQ_B(2, 3);
Q_03 = NextQ(Prob, Q_06, ThetaVec);

%
Prob = SortedDataQ_A(2, 3);
Q_07 = NextQ(Prob, Q_03, ThetaVec);

%
Prob = SortedDataQ_B(3, 3);
Q_02 = NextQ(Prob, Q_07, ThetaVec);

%
Prob = SortedDataQ_A(3, 3);
Q_08 = NextQ(Prob, Q_02, ThetaVec);

%
Prob = SortedDataQ_B(4, 3);
Q_01 = NextQ(Prob, Q_08, ThetaVec);

%
Prob = SortedDataQ_A(4, 3);
Q_09 = NextQ(Prob, Q_01, ThetaVec);

%__%
Q_y = [Q_09, Q_08, Q_07, Q_06, Q_05, Q_04, Q_03, Q_02, Q_01, 0];
Q_x = -slope_W* [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0];   

end

function [Q_Next] = NextQ(Prob, Q_Previous, ThetaVec)

Pi_Prob = interp1([0, 0.25, 0.5, 0.75, 1], ThetaVec, Prob);

Q_Next = Pi_Prob/(1-Pi_Prob) * Q_Previous;
end
