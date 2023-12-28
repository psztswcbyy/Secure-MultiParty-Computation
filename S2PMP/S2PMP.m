% ============================================================================================================================
% Secure_2PMP函数用于实现安全两方矩阵乘法协议，涉及到的输入参数包含参与计算的方阵维度N，参与方阵的元素阶数最小值minEp，参与方
% 阵的元素阶数最大值maxEp，参与方阵元素个位数最小值FirstNumMin，参与方阵元素个位数最大值FirstNumMax；涉及到的输出参数包含第一
% 参与方Alice的初始矩阵A，第二参与方Bob的初始矩阵B，两方矩阵乘法结果的相对误差矩阵Re_Err_result，两方矩阵乘法结果的绝对误差矩
% 阵Ab_Err_result，乘法结果的最大相对误差结果MRE_result，乘法结果的均值绝对百分比误差MAPE_result，乘法结果的绝对百分比误差总和
% SAPE_result,乘法结果的F范数相对误差F_norm，理论的乘法结果Theory_V，实际的乘法结果Real_V
% ============================================================================================================================
% clc;
% clear;
% format longE
% % ============================================================================================================================
% minEp = -10;
% maxEp = 10;
% FirstNumMin = 1;
% FirstNumMax = 1;
% % ============================================================================================================================
% % 以下是对于任意两个维度的矩阵测试
% DimA = randi([2,10],1,2);
% A_row = DimA(1);
% A_col = DimA(2);
% B_row = A_col;
% B_col = randi([2,10],1);
% % ============================================================================================================================
% % 以下是对于方阵的测试
% Com_Dim = 5;
% A_row = Com_Dim;
% A_col = Com_Dim;
% B_row = Com_Dim;
% B_col = Com_Dim;
% % ============================================================================================================================
%
% Ori_A = rand(A_row,A_col,'double') + randi([FirstNumMin,FirstNumMax],A_row,A_col);
% Exp_A = 10.^(randi([minEp, maxEp],A_row,A_col));
% A = Ori_A.*Exp_A;
%
% Ori_B = rand(B_row,B_col,'double') + randi([FirstNumMin,FirstNumMax],B_row,B_col);
% Exp_B = 10.^(randi([minEp, maxEp],B_row,B_col));
% B = Ori_B.*Exp_B;
%
% [Va, Vb, A, B, Re_Err_result, Ab_Err_result, MRE_result, MAPE_result, SAPE_result, F_norm, Theory_V, Real_V] =...
%     S_2PMP(A, B)
% ============================================================================================================================
% 定义Secure_2PMP安全两方乘法协议
function [Va, Vb, A, B, Re_Err_result, Ab_Err_result, MRE_result, MAPE_result, SAPE_result, F_norm, Theory_V, Real_V] =...
    S2PMP(A, B)
IsNulErrorOccur = false;

while(~IsNulErrorOccur)
    % ============================================================================================================================
    % 根据输入矩阵A和B生成随机矩阵RA和RB，混淆矩阵A_hat和B_hat
    % ============================================================================================================================
    [Ra, A_hat] = RMG_2PMP(A);
    [Rb, B_hat] = RMG_2PMP(B);
    
    % ============================================================================================================================
    % 生成随机矩阵ra和rb以及Bob的私有随机矩阵Vb
    [~, ra] = RMG_2PMP(Ra*Rb);
    [~, Vb]  =RMG_2PMP(ra);
    rb = Ra * Rb - ra;
    
    % Bob生成T矩阵
    T = A_hat * B + (rb - Vb);
    
    % 随机矩阵Va
    Va = T + ra - (Ra * B_hat);
    
    % 理论结果Theory_V=A*B,实验结果Real_V=Va+Vb
    Theory_V = A * B;
    Real_V = Va + Vb;
    
    % ============================================================================================================================
    % 矩阵精度优化模块
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
    % ============================================================================================================================
    % 误差参数测试模块
    % ============================================================================================================================
    
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
    
end
end
