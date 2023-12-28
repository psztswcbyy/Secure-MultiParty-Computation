clc;
clear;
format longE
% =======================================================================================================================================================================
%*******************************************************【最终版本—多轮试验】*******************************************************************************************
% 定义初始参数，包含矩阵维度区间[DimMin,DimMax],阶数上界范围[EpMin,EpMax]，原始数据个位元素数值范围[FirstNumMin，FirstNumMax]
% 阶数下界范围[-EpMax,-EpMin],单次实验测试轮数Num_SingTest，维度实验条件数Num_Dim，阶数实验条件数Num_Ep。
% 输出①MRE最大相对误差全局实验结果MRE_full；②MAPE均值绝对百分比误差全局实验结果MAPE_full；③Q分位数全局实验结果Q_full。
ProgammeStart = clock;
% 进度条初始化模块
fig = uifigure('Name','Results');
d = uiprogressdlg(fig,'Title','Please Wait',...
    'Message','Computation begins execution','Cancelable','on');
drawnow
pause(.8)

DimMin = 10;
DimMax = 50;
DimUnit = 10;

% ①对称区间测试案例EpMin=0,EpMax = 16; EpUnit = 1; maxEp = EpMin:EpUnit:EpMax; minEp = -maxEp; j = (maxEp-EpMin)/EpUnit + 1;
% ②正区间测试案例EpMin=0,EpMax = 16; EpUnit = 1; maxEp = EpMin:EpUnit:EpMax; minEp = 0; j = (maxEp-EpMin)/EpUnit + 1;
% ③负区间测试案例EpMin=0,EpMax = -16; EpUnit = -1; minEp = EpMin:EpUnit:EpMax; maxEp = 0; j = (minEp-EpMin)/EpUnit + 1;


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
        
        % 进度step声明
        step = (N - DimMin)/DimUnit * Num_Ep + (maxEp - EpMin)/EpUnit+1;
        % step = (N - DimMin)/DimUnit * Num_Ep + (minEp - EpMin)/EpUnit+1;
        
        RLT_Error = [];
        ABS_Error = [];
        Theory_V = [];
        Real_V = [];
        A_copy = [];
        B_copy = [];
        
        MRE_Error = zeros(1,Num_SingleTest,'double');
        SAPE_Error = zeros(1,Num_SingleTest,'double');
        Fnorm_Error = zeros(1,Num_SingleTest,'double');
        
        
        for i = 1:Num_SingleTest
            A = RandMatrixGen2Condition(N,minEp, maxEp, FirstNumMin, FirstNumMax);
            B = RandMatrixGen2Condition(N,minEp, maxEp, FirstNumMin, FirstNumMax);
            
            [Va, Vb, A_result, B_result, Cond_T, Re_Err_result, Ab_Err_result, MRE_result, MAPE_result, SAPE_result, F_norm, ...
                Theory_result, Real_result, Error_Rate] = S2PIP_Distinguish(A, B);
            
            Theory_V = cat(3,Theory_V,Theory_result);
            Real_V = cat(3,Real_V,Real_result);
            
            A_copy = cat(3,A_copy,A_result);
            B_copy = cat(3,B_copy,B_result);
            
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
        for i = 1:length(x)
            Theory_Value = Theory_V(x(i), y(i), index_TestNum);
            Outlier_Value = Real_V(x(i), y(i), index_TestNum);
            ABS_max = ABS_Error(x(i),y(i),index_TestNum);
        end
        A = A_copy(:, :, index_TestNum);
        B = B_copy(:, :, index_TestNum);
        
        Theory_Matrix = Theory_V(:, :, index_TestNum);
        Outlier_Matrix = Real_V(:, :, index_TestNum);
        % 查看第frame次实验的相对误差分布
        frame = index_TestNum;
        RelativeError = RLT_Error(:,:,frame);
        
        % 将每次的百分比绝对误差和SAPE_Error累加除以所有实验包含的元素个数
        MAPE = sum(SAPE_Error(:))/(Num_SingleTest*N*N);
        
        % 分位数求解Q1,Q2,Q3,Q4,Q5
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
        % 取消按钮检查
        if d.CancelRequested
            break
        end
        
        % 进度更新模块
        d.Value = step/steps;
        d.Message = sprintf('Computation Process in %.1f%%',100*d.Value);
        pause(.001)
    end
end

% =======================================================================================================================================================================
% 将结果写入EXCEL文件中
filename = strcat('【Ep=(', num2str(minEp),',',num2str(maxEp),')Dim=(',...
    num2str(DimMin),',',num2str(DimMax),')】S3PMP_ExperimentResults.xlsx');
xlswrite(filename, MRE_full,1,'C4');
xlswrite(filename, MAPE_full,1,'C12');
for i=1:Num_Ep
    X_range = strcat('C',num2str(20+(i-1)*7));
    xlswrite(filename, Q_full(:,:,i),1,X_range);
end

% 时间记录模块
ProgammeEnd=clock;
TotalTime = etime(ProgammeEnd, ProgammeStart);
fprintf('\nThe Programme Execution Time is %.5fs\n',TotalTime);
% 进度条关闭
close(d)
delete(fig)

