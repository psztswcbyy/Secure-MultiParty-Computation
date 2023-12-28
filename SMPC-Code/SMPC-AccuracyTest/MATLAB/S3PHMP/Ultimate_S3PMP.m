clc;
clear;
format longE
% =======================================================================================================================================================================
%*******************************************************�����հ汾���������顿*******************************************************************************************
% �����ʼ��������������ά������[DimMin,DimMax],�����Ͻ緶Χ[EpMin,EpMax]��ԭʼ���ݸ�λԪ����ֵ��Χ[FirstNumMin��FirstNumMax]
% �����½緶Χ[-EpMax,-EpMin],����ʵ���������Num_SingTest��ά��ʵ��������Num_Dim������ʵ��������Num_Ep��
% �����MRE���������ȫ��ʵ����MRE_full����MAPE��ֵ���԰ٷֱ����ȫ��ʵ����MAPE_full����Q��λ��ȫ��ʵ����Q_full��
ProgammeStart = clock;
% ��������ʼ��ģ��
fig = uifigure('Name','Results');
d = uiprogressdlg(fig,'Title','Please Wait',...
    'Message','Computation begins execution','Cancelable','on');
drawnow
pause(.8)

DimMin = 10;
DimMax = 50;
DimUnit = 10;

% �ٶԳ�������԰���EpMin=0,EpMax = 16; EpUnit = 1; maxEp = EpMin:EpUnit:EpMax; minEp = -maxEp; j = (maxEp-EpMin)/EpUnit + 1;
% ����������԰���EpMin=0,EpMax = 16; EpUnit = 1; maxEp = EpMin:EpUnit:EpMax; minEp = 0; j = (maxEp-EpMin)/EpUnit + 1;
% �۸�������԰���EpMin=0,EpMax = -16; EpUnit = -1; minEp = EpMin:EpUnit:EpMax; maxEp = 0; j = (minEp-EpMin)/EpUnit + 1;


EpMin = 0;
EpMax = 16;
EpUnit = 1;

FirstNumMin = 1;
FirstNumMax = 1;

Num_SingleTest = 100;
Num_Dim = ((DimMax - DimMin)/DimUnit + 1);
Num_Ep = ((EpMax-EpMin)/EpUnit + 1);

MRE_full = zeros(Num_Dim, Num_Ep, 'double');
MAPE_full = zeros(Num_Dim, Num_Ep, 'double');
Q_full = zeros(Num_Dim, 5, Num_Ep, 'double');
steps = Num_Dim * Num_Ep;

for N = DimMin:DimUnit:DimMax
    for maxEp = EpMin:EpUnit:EpMax
%         minEp = -maxEp;
        minEp = 0;
%         maxEp = 0;
        
        % ����step����
        step = (N - DimMin)/DimUnit * Num_Ep + (maxEp - EpMin)/EpUnit+1;
        % step = (N - DimMin)/DimUnit * Num_Ep + (minEp - EpMin)/EpUnit+1;

        RLT_Error = [];
        ABS_Error = [];
        Theory_V = [];
        Real_V = [];
        A_copy = [];
        B_copy = [];
        C_copy = [];
        
        MRE_Error = zeros(1,Num_SingleTest,'double');
        SAPE_Error = zeros(1,Num_SingleTest,'double');
        Fnorm_Error = zeros(1,Num_SingleTest,'double');
        
        
        for i = 1:Num_SingleTest
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
        % ���Num��ʵ������������MRE,����ʾ��Ӧ�����ֵ�ʵ����������index_TestNum����Ӧ����MRE�����������RLT�������ڸþ�����
        % ����(x,y),��ʾ��ӦMREԪ�ص�����ֵTheory_Value���쳣ֵOutlier_Value����ӦMRE�ľ������ABS_max�����۵ĳ˷����Theory_Matrix��
        % ʵ�ʵĳ˷����Outlier_Matrix������ʾ��Ӧʵ�����������������RelativeError
        [MRE, index_TestNum] = max(MRE_Error);
        RLT = RLT_Error(:,:,index_TestNum);
        [x, y] = find(RLT == MRE);
        for i = 1:length(x)
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
        MAPE = sum(SAPE_Error(:))/(Num_SingleTest*N*N);
        
        % ��λ�����Q1,Q2,Q3,Q4,Q5
        p = [0,25,50,75,100];
        Q = prctile(Fnorm_Error, p);
        
        %         fprintf("MRE  = \n     %16.15E\nMAPE = \n     %16.15E\n", MRE,MAPE);
        %         fprintf("Q    = \n     %16.15E     %16.15E     %16.15E     %16.15E     %16.15E\n",Q(1),Q(2),Q(3),Q(4),Q(5));
        
        i = (N - DimMin)/DimUnit + 1;
        j = (maxEp-EpMin)/EpUnit + 1;
        %         j = (minEp-EpMin)/EpUnit + 1;
        MRE_full(i,j) = MRE;
        MAPE_full(i,j) = MAPE;
        for k=1:5
            Q_full(i,k,j) = Q(k);
        end
        
        % ȡ����ť���
        if d.CancelRequested
            break
        end
        
        % ���ȸ���ģ��
        d.Value = step/steps;
        d.Message = sprintf('Computation Process in %.1f%%',100*d.Value);
        pause(.001)
    end
end

% =======================================================================================================================================================================
% �����д��EXCEL�ļ���
filename = strcat('��Ep=(', num2str(minEp),',',num2str(maxEp),')Dim=(',...
    num2str(DimMin),',',num2str(DimMax),')��S3PMP_ExperimentResults.xlsx');
xlswrite(filename, MRE_full,1,'C4');
xlswrite(filename, MAPE_full,1,'C12');
for i=1:Num_Ep
    X_range = strcat('C',num2str(20+(i-1)*7));
    xlswrite(filename, Q_full(:,:,i),1,X_range);
end

% ʱ���¼ģ��
ProgammeEnd=clock;
TotalTime = etime(ProgammeEnd, ProgammeStart);
fprintf('\nThe Programme Execution Time is %.5fs\n',TotalTime);
% �������ر�
close(d)
delete(fig)

