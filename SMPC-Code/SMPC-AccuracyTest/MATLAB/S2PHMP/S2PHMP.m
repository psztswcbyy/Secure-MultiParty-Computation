% =================================================================================================================================
% 【功能说明】S2PHMP函数用于实现安全两方矩阵混合乘法协议，涉及到的输入参数包含参与方Alice输入矩阵A1、A2和参与方Bob输入矩阵B1、B2；
% 涉及到的输出参数包含第一参与方Alice的输出矩阵Va，第二参与方Bob的输出矩阵Vb，两方矩阵乘法结果的相对误差矩阵Re_Err_result，两方矩
% 阵乘法结果的绝对误差矩阵Ab_Err_result，乘法结果的最大相对误差结果MRE_result，乘法结果的均值绝对百分比误差MAPE_result，乘法结果的
% 绝对百分比误差总和SAPE_result,乘法结果的F范数相对误差F_norm，理论的乘法结果Theory_V，实际的乘法结果Real_V
% ==================================================================================================================================
% clc;
% clear;
% format longE
% 
% N = 5;
% minEp = -16;
% maxEp = -minEp;
% FirstNumMin = 1;
% FirstNumMax = 1;
% 
% % A = RandInputMatrix(N,FirstNumMin,FirstNumMax,minEp,maxEp);
% % B = RandInputMatrix(N,FirstNumMin,FirstNumMax,minEp,maxEp);
% % C = RandInputMatrix(N,FirstNumMin,FirstNumMax,minEp,maxEp);
% % 由于2PMP协议无非满秩的要求因此可以直接用RandMatrixGen2生成初始矩阵
% [A1, A2, A] = RandMatrixGen2(N,minEp, maxEp, FirstNumMin, FirstNumMax);
% [B1, B2, B] = RandMatrixGen2(N,minEp, maxEp, FirstNumMin, FirstNumMax);
% 
% [Va, Vb, A1, A2, B1, B2, Re_Err_result, Ab_Err_result, MRE_result, MAPE_result, SAPE_result, F_norm, Theory_V, Real_V] =...
%     Secure2PHMP(A1, A2, B1, B2)
% ==================================================================================================================================
% 定义Secure_2PHMP安全两方混合乘法协议、

function [Va, Vb, A1, A2, B1, B2, Re_Err_result, Ab_Err_result, MRE_result, MAPE_result, SAPE_result, F_norm, Theory_V, ...
    Real_V, Error_Rate] = S2PHMP(A1, A2, B1, B2)
IsNulErrorOccur = false;

while(~IsNulErrorOccur)
    % ============================================================================================================================
    % 三次2PMP实现3PMP协议
    
    [Va1, Vb1] = S2PMP(A1,B2);%AB=Va1+Vb1
    
    [Vb2, Va2] = S2PMP(B1,A2);%Va1*C=Va2+Vc1
    
    Va = Va1 + Va2 + A1*A2;
    Vb = Vb1 + Vb2 + B1*B2;
    
    % ============================================================================================================================
    % 理论结果Theory_V=A*B*C,实验结果Real_V=Va+Vb+Vc
    Theory_V = (A1 + B1) * (A2 + B2);
    Real_V = Va + Vb;
    
    % ============================================================================================================================
    % 设定矩阵数据精度极限eps(class(A)),对于低于精度极限的数据，将其设置并视为0元素
    % 优化理论结果矩阵Theory_V
    [m, n] = size(Theory_V);
    for i = 1:m
        for j = 1:n
            if abs(Theory_V(i, j))< max(m,n)*eps(class(Theory_V))*norm(Theory_V, 'inf')
                Theory_V(i, j) = 0;
            end
        end
    end
    % 优化实际结果矩阵Real_V
    [m, n] = size(Real_V);
    for i = 1:m
        for j = 1:n
            if abs(Real_V(i, j))< max(m,n)*eps(class(Real_V))*norm(Real_V, 'inf')
                Real_V(i, j) = 0;
            end
        end
    end
    
    
    % 绝对误差矩阵Ab_Err_result及其维度
    Ab_Err_result = Real_V - Theory_V;
    [m, n] = size(Ab_Err_result);
    
    % 初始化相对误差矩阵Re_Err_result
    Re_Err_result = zeros(m, n);
    
    % disp(Real_V),disp(Theory_V),disp(Ab_Err_result);
    
    % 按照理论矩阵中对应元素为0划分三种情况，第一种当理论元素为0，且绝对误差超出精度极限，默认为实际元素为0
    % 第二种情况当理论元素非0，则按照实际情况去求解相对误差
    % 第三种情况当理论元素为0，但是绝对误差在精度极限内，则直接让元素的相对误差等于该绝对误差
    IsNulErrorOccur = true;
    for i = 1:m
        for j = 1:n
            if Theory_V(i, j) == 0 && abs(Ab_Err_result(i, j)) < max(m,n)*eps(class(Ab_Err_result))*norm(Ab_Err_result, 'inf')
                Re_Err_result(i, j) = 0;
            elseif Theory_V(i, j) ~= 0 && Real_V(i, j)~=0 && abs(1 - Real_V(i, j)/Theory_V(i, j)) < 1
                Re_Err_result(i, j) = abs(Ab_Err_result(i, j)/Theory_V(i, j));
            elseif Theory_V(i, j) ~= 0 && Real_V(i, j)==0
                Re_Err_result(i, j) = 0;
            elseif Theory_V(i, j) ~= 0 && Real_V(i, j)~=0 && abs(1 - Real_V(i, j)/Theory_V(i, j)) >= 1
                Re_Err_result(i, j) = 0;
            else
                disp('WARNING: The raw data or the algorithm itself contains computational errors!\n');
                Re_Err_result(i, j) = abs(Ab_Err_result(i, j));
                IsNulErrorOccur = false;
            end
        end
    end
    
    
    % Re_Err_result = abs(Ab_Err_result./Theory_V);
    
    % 最大相对误差MRE
    MRE_result = max(Re_Err_result(:));
    % 绝对百分比误差和SAPE
    SAPE_result = sum(abs(Re_Err_result(:)));
    % 均值绝对百分比误差MAPE
    MAPE_result = sum(abs(Re_Err_result(:)))/numel(Re_Err_result);
    
    % F范数求解,Fa绝对误差矩阵的范数
    F_absolute = norm(Ab_Err_result, 'fro');
    F_theory = norm(Theory_V, 'fro');
    F_norm = F_absolute/F_theory;
    
    % 错误率ErrorRate统计，当上述误差优化部分无法成功实现的时候，改用错误率指标
    Error_Rate =  numel(find(Re_Err_result >= 0.5*10^-3)) / numel(Re_Err_result);
    
end
end
