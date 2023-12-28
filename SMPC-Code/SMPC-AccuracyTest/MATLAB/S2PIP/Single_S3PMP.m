% clc;
clear;
format longE
ProgammeStart = clock;
% =======================================================================================================================================================================
%**********************************************************【初级版本-单次试验】*******************************************************************
% 定义初始参数，包含矩阵维度N,阶数范围minEP和maxEP，＼原始数据首位范围FirstNumMin，FirstNumMax及实验遍历次数Num
N = 10;
minEp = -16;
maxEp = -minEp;
FirstNumMin = 1;
FirstNumMax = 1;
Num = 300;
% =======================================================================================================================================================================
% 获取并备份Num次实验的结果,RLT_Error相对误差矩阵备份，ABS_Error绝对误差矩阵备份，Theory_V理论V=Va+Vb矩阵备份，Real_V实际V=Va+Vb矩阵备份，A_copy原始Alice参与方数据矩阵备
% 份，B_copy原始Bob参与方数据矩阵备份，MRE_Error最大相对误差备份，SAPE_Error绝对值百分比误差和备份，Fnorm_Error是F范数备份
RLT_Error = [];
ABS_Error = [];
Theory_V = [];
Real_V = [];
A_copy = [];
B_copy = [];
C_copy = [];
MRE_Error = zeros(1,Num,'double');
SAPE_Error = zeros(1,Num,'double');
Fnorm_Error = zeros(1,Num,'double');
% 循环Num次实验
for i = 1:Num
    [A_result, B_result, C_result, Re_Err_result, Ab_Err_result, MRE_result, MAPE_result, SAPE_result, F_norm, Theory_result, Real_result] = ...
        Secure_3PMP(N, minEp, maxEp, FirstNumMin, FirstNumMax);
    
    Theory_V = cat(3,Theory_V,Theory_result);
    Real_V = cat(3,Real_V,Real_result);
    
    A_copy = cat(3,A_copy,A_result);
    B_copy = cat(3,B_copy,B_result);
    C_copy = cat(3,C_copy,C_result);
    
    RLT_Error = cat(3,RLT_Error,Re_Err_result);
    ABS_Error = cat(3,ABS_Error,Ab_Err_result);
    
    MRE_Error(i) = MRE_result;
    SAPE_Error(i) = SAPE_result;
    Fnorm_Error(i) = F_norm;
end

% =======================================================================================================================================================================
% 求解Num次实验的最大相对误差MRE,并显示对应误差出现的实验索引次数index_TestNum，对应存在MRE的相对误差矩阵RLT，返回在该矩阵中
% 坐标(x,y),显示对应MRE元素的理论值Theory_Value、异常值Outlier_Value，对应MRE的绝对误差ABS_max，理论的乘法结果Theory_Matrix，
% 实际的乘法结果Outlier_Matrix，并显示对应实验索引的相对误差矩阵RelativeError
[MRE, index_TestNum] = max(MRE_Error);
RLT = RLT_Error(:,:,index_TestNum);
[x, y] = find(RLT == MRE);
for i=1:length(x)
    Theory_Value = Theory_V(x(i), y(i), index_TestNum);
    Outlier_Value = Real_V(x(i), y(i), index_TestNum);
    ABS_max = ABS_Error(x(i),y(i),index_TestNum);
end
A = A_copy(:, :, index_TestNum);
B = B_copy(:, :, index_TestNum);
C = C_copy(:, :, index_TestNum);

Theory_Matrix = Theory_V(:, :, index_TestNum);
Outlier_Matrix = Real_V(:, :, index_TestNum);
% 查看第frame次实验的相对误差分布
frame = index_TestNum;
RelativeError = RLT_Error(:,:,frame);

% 将每次的百分比绝对误差和SAPE_Error累加除以所有实验包含的元素个数
MAPE = sum(SAPE_Error(:))/(Num*N*N);

% 分位数求解Q1,Q2,Q3,Q4,Q5
p = [0,25,50,75,100];
Q = prctile(Fnorm_Error, p);

fprintf("MRE  = \n     %16.15E\nMAPE = \n     %16.15E\n", MRE,MAPE);
fprintf("Q    = \n     %16.15E     %16.15E     %16.15E     %16.15E     %16.15E\n",Q(1),Q(2),Q(3),Q(4),Q(5));
% =======================================================================================================================================================================
% 将结果写入文本文件
if (minEp + maxEp)<=0
    delta = minEp;
else
    delta = maxEp;
end

filename = strcat('【Ep=(', num2str(minEp),',',num2str(maxEp),'),Dim=(',...
    num2str(N),',',num2str(N),')】S3PMP_ExperimentResults.txt');
fileID = fopen(filename,'w');
fprintf(fileID,'*****************************************************************************************************************************\n');
fprintf(fileID,'The Test based on Dim=(%2d,%2d) and the delta is(E%3d ->E%3d),the repeat time=%d\r\n', N, N, minEp, maxEp, Num);
fprintf(fileID,'=============================================================================================================================\n');
fprintf(fileID,'%-43s%-57s%-44s\r\n','MRE','Index of MRE Number','MAPE');
fprintf(fileID,'-----------------------------------------------------------------------------------------------------------------------------\n');
fprintf(fileID,'%-50.15E%-50d%-50.15E\r\n',MRE, index_TestNum, MAPE);
fprintf(fileID,'=============================================================================================================================\n');
fprintf(fileID,'%-25s%-25s%-25s%-25s%-25s\r\n','Q1','Q2','Q3','Q4','Q5');
fprintf(fileID,'-----------------------------------------------------------------------------------------------------------------------------\n');
fprintf(fileID,'%-25.15E%-25.15E%-25.15E%-25.15E%-25.15E\n',Q);
fclose(fileID);

ProgammeEnd=clock;
TotalTime = etime(ProgammeEnd, ProgammeStart);
fprintf('\nThe Programme Execution Time is %.5fs\n',TotalTime);
