% ============================================================================================================================
% Secure_3PMP函数用于实现安全三方矩阵乘法协议，涉及到的输入参数包含参与计算的方阵维度N，参与方阵的元素阶数最小值minEp，参与方
% 阵的元素阶数最大值maxEp，参与方阵元素个位数最小值FirstNumMin，参与方阵元素个位数最大值FirstNumMax；涉及到的输出参数包含第一
% 参与方Alice的初始矩阵A，第二参与方Bob的初始矩阵B，第三参与方Carol的初始矩阵C，三方矩阵乘法结果的相对误差矩阵Re_Err_result，
% 三方矩阵乘法结果的绝对误差矩阵Ab_Err_result，乘法结果的最大相对误差结果MRE_result，乘法结果的均值绝对百分比误差MAPE_result，
% 乘法结果的绝对百分比误差总和SAPE_result，乘法结果的F范数相对误差F_norm，理论的乘法结果Theory_V，实际的乘法结果Real_V
% ============================================================================================================================
% clc;
% clear;
% format longE
%
% N = 5;
% minEp = -8;
% maxEp = 8;
% FirstNumMin = 1;
% FirstNumMax = 1;
%
% % A = RandInputMatrix(N,FirstNumMin,FirstNumMax,minEp,maxEp);
% % B = RandInputMatrix(N,FirstNumMin,FirstNumMax,minEp,maxEp);
% % C = RandInputMatrix(N,FirstNumMin,FirstNumMax,minEp,maxEp);
% 由于2PMP协议无非满秩的要求因此可以直接用RandMatrixGen2生成初始矩阵
% A = RandMatrixGen2(N,minEp, maxEp, FirstNumMin, FirstNumMax);
% B = RandMatrixGen2(N,minEp, maxEp, FirstNumMin, FirstNumMax);
% C = RandMatrixGen2(N,minEp, maxEp, FirstNumMin, FirstNumMax);
%
% [Va, Vb, Vc, A, B, C, Re_Err_result, Ab_Err_result, MRE_result, MAPE_result, SAPE_result, F_norm, Theory_V, Real_V] =...
%     Secure3x2PMP(A, B, C)
% ============================================================================================================================
% 定义Secure_3PMP安全三方乘法协议
function [Va, Vb, Vc, A, B, C, Re_Err_result, Ab_Err_result, MRE_result, MAPE_result, SAPE_result, F_norm, Theory_V, Real_V] =...
    Secure_3x2PMP(A, B, C)
IsNulErrorOccur = false;

while(~IsNulErrorOccur)
    % ============================================================================================================================
    % 三次2PMP实现3PMP协议
    
    [Va1, Vb1] = S2PMP(A,B);%AB=Va1+Vb1
    
    [Va2, Vc1] = S2PMP(Va1,C);%Va1*C=Va2+Vc1
    
    [Vb2, Vc2] = S2PMP(Vb1,C);%Vb1*C=Vb2+Vc2
    
    Va = Va2;
    Vb = Vb2;
    Vc = Vc1 + Vc2;
    
    % ============================================================================================================================
    % 理论结果Theory_V=A*B*C,实验结果Real_V=Va+Vb+Vc
    Theory_V = A * B * C;
    Real_V = Va + Vb + Vc;
    
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
    
end
end

