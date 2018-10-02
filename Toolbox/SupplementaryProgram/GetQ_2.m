function [ Q_x, Q_y ] = GetQ_2( ResultArray , slope_W, Q_scale, ColumnNum)


DataQ = [ResultArray{3, ColumnNum}.Question(2:4);
         ResultArray{4,ColumnNum}.Question(2:4);
         ResultArray{6,ColumnNum}.Question(2:4);
         ResultArray{7,ColumnNum}.Question(2:4);
         ResultArray{8,ColumnNum}.Question(2:4);
         ResultArray{9,ColumnNum}.Question(2:4);
         ResultArray{10,ColumnNum}.Question(2:4);
         ResultArray{11,ColumnNum}.Question(2:4)];
         
%          SetQ{1,10}.Question(2:4);
%          SetQ{1,11}.Question(2:4);
%          SetQ{1,12}.Question(2:4);
%          SetQ{1,13}.Question(2:4);
%          SetQ{1,14}.Question(2:4)];

SortedDataQ = sortrows(DataQ,  2);

SortedDataQ_A = SortedDataQ(1:4,:);
SortedDataQ_B = SortedDataQ(5:8,:);

%Q_y        = Q_scale * SortedDataQ(:,3) ./(1-SortedDataQ(:,3));
Q_05_positive = Q_scale;

Prob = SortedDataQ_B(1, 3);
% Q_04_negative = NextQ(Prob, Q_05_positive);
Q_04_negative = -Prob/(1-Prob) * Q_05_positive;

Q_04_positive = interp1([0, 0.5], [0, 1], 0.4);

Prob = SortedDataQ_A(1, 3);
% % Q_06 = Pi_Prob/(1-Pi_Prob) * Q_04;
Q_06_negative = -Prob/(1-Prob) * Q_04_positive; 
Q_06_positive = interp1([0, 0.5], [0,1], 0.6, 'linear', 'extrap');

%
Prob = SortedDataQ_B(2, 3);
Q_03_negative = -Prob/(1-Prob) * Q_06_positive;
Q_03_positive = interp1([0, 0.5], [0,1], 0.3, 'linear', 'extrap');

%
Prob = SortedDataQ_A(2, 3);
Q_07_negative = -Prob/(1-Prob) * Q_03_positive;
Q_07_positive = interp1([0, 0.5], [0,1], 0.7, 'linear', 'extrap');

%
Prob = SortedDataQ_B(3, 3);
Q_02_negative = -Prob/(1-Prob) * Q_07_positive;
Q_02_positive = interp1([0, 0.5], [0,1], 0.2, 'linear', 'extrap');

%
Prob = SortedDataQ_A(3, 3);
Q_08_negative = -Prob/(1-Prob) * Q_02_positive;
Q_08_positive = interp1([0, 0.5], [0,1], 0.8, 'linear', 'extrap');

%
Prob = SortedDataQ_B(4, 3);
Q_01_negative = -Prob/(1-Prob) * Q_08_positive;
Q_01_positive = interp1([0, 0.5], [0,1], 0.1, 'linear', 'extrap');

%
Prob = SortedDataQ_A(4, 3);
Q_09_negative = -Prob/(1-Prob) * Q_01_positive;
% % Q_09 = Prob/(1-Prob) * Q_01;
% Q_09 = NextQ(Prob, Q_01);

%___%
Q_y = [Q_09_negative, Q_08_negative, Q_07_negative, Q_06_negative, Q_04_negative, Q_03_negative, Q_02_negative, Q_01_negative, 0, Q_01_positive, Q_02_positive, Q_03_positive, Q_04_positive, Q_05_positive, Q_06_positive, Q_07_positive, Q_08_positive];
%Q_y = [Q_09_negative, Q_08_negative, Q_07_negative, Q_06_negative, Q_04_negative, Q_03_negative, Q_02_negative, Q_01_negative, 0];
%Q_x = [-0.9, -0.8, -0.7, -0.6, -0.4, -0.3, -0.2, -0.1, 0];  
Q_x = [-0.9, -0.8, -0.7, -0.6, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5,0.6, 0.7, 0.8];  

%___%


% 
% Prob = SortedDataQ_B(4, 3);
% % Q_01 = Prob/(1-Prob) * Q_08;
% Q_01 = NextQ(Prob, Q_08);
% 
% Prob = SortedDataQ_A(4, 3);
% % Q_09 = Prob/(1-Prob) * Q_01;
% Q_09 = NextQ(Prob, Q_01);
% 
% 
% Q_x = -slope_W* [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0];   

end

function [Q_Next] = NextQ(Prob, Q_Previous)

Pi_Prob = Prob;
Q_Next = - Pi_Prob/(1-Pi_Prob) * Q_Previous;
end
