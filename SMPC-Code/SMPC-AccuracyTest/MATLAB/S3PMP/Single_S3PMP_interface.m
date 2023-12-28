clc;
clear;
format longE
ProgammeStart = clock;
% =======================================================================================================================================================================
%**********************************************************�������汾-�������顿*******************************************************************
% �����ʼ��������������ά��N,������ΧminEP��maxEP����ԭʼ������λ��ΧFirstNumMin��FirstNumMax��ʵ���������Num
N = 10;
minEp = -16;
maxEp = -minEp;
FirstNumMin = 1;
FirstNumMax = 1;
Num = 500;
% =======================================================================================================================================================================
% ��ȡ������Num��ʵ��Ľ��,RLT_Error��������󱸷ݣ�ABS_Error���������󱸷ݣ�Theory_V����V=Va+Vb���󱸷ݣ�Real_Vʵ��V=Va+Vb���󱸷ݣ�A_copyԭʼAlice���뷽���ݾ���
% �ݣ�B_copyԭʼBob���뷽���ݾ��󱸷ݣ�MRE_Error���������ݣ�SAPE_Error����ֵ�ٷֱ����ͱ��ݣ�Fnorm_Error��F��������
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

% ѭ��Num��ʵ��
for i = 1:Num
    A = RandMatrixGen3(N,minEp, maxEp, FirstNumMin, FirstNumMax);
    B = RandMatrixGen3(N,minEp, maxEp, FirstNumMin, FirstNumMax);
    C = RandMatrixGen3(N,minEp, maxEp, FirstNumMin, FirstNumMax);
    
    [Va, Vb, Vc, A_result, B_result, C_result, Diff_A, Diff_B, Diff_C, Re_Err_result, Ab_Err_result, MRE_result, ...
        MAPE_result, SAPE_result, F_norm, Theory_result, Real_result, Error_Rate] = S3PMP(A, B, C);
    
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
% ���Num��ʵ������������MRE,����ʾ��Ӧ�����ֵ�ʵ����������index_TestNum����Ӧ����MRE�����������RLT�������ڸþ�����
% ����(x,y),��ʾ��ӦMREԪ�ص�����ֵTheory_Value���쳣ֵOutlier_Value����ӦMRE�ľ������ABS_max�����۵ĳ˷����Theory_Matrix��
% ʵ�ʵĳ˷����Outlier_Matrix������ʾ��Ӧʵ�����������������RelativeError
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
% �鿴��frame��ʵ���������ֲ�
frame = index_TestNum;
RelativeError = RLT_Error(:,:,frame);

% ��ÿ�εİٷֱȾ�������SAPE_Error�ۼӳ�������ʵ�������Ԫ�ظ���
MAPE = sum(SAPE_Error(:))/(Num*N*N);

% ��λ�����Q1,Q2,Q3,Q4,Q5
p = [0,25,50,75,100];
Q = prctile(Fnorm_Error, p);

fprintf("MRE  = \n     %16.15E\nMAPE = \n     %16.15E\n", MRE,MAPE);
fprintf("Q    = \n     %16.15E     %16.15E     %16.15E     %16.15E     %16.15E\n",Q(1),Q(2),Q(3),Q(4),Q(5));
% =======================================================================================================================================================================
% �����д���ı��ļ�
if (minEp + maxEp)<=0
    delta = minEp;
else
    delta = maxEp;
end

filename = strcat('��Ep=(', num2str(minEp),',',num2str(maxEp),'),Dim=(',...
    num2str(N),',',num2str(N),')��S3PMP_ExperimentResults.txt');
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
