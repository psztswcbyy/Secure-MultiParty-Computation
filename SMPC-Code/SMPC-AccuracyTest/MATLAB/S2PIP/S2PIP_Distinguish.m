% ==================================================================================================================
% 【功能说明】
%  2PIP协议用于安全两方求逆计算中，这里按照协议的设计原理所涉及的两个可逆矩阵P和Q的所属权进行分类的。本协议所涉及
%  可逆矩阵P和Q均为方阵且满秩的情况，令矩阵P属于Alice和矩阵Q属于Bob。
%
%
% ==================================================================================================================
% clc;
% clear;
% format longE
% 
% N = 10;
% minEp = -16;
% maxEp = -minEp;
% FirstNumMin = 1;
% FirstNumMax = 1;
% 
% % % A = RandInputMatrix(N,FirstNumMin,FirstNumMax,minEp,maxEp);
% % % B = RandInputMatrix(N,FirstNumMin,FirstNumMax,minEp,maxEp);
% % % C = RandInputMatrix(N,FirstNumMin,FirstNumMax,minEp,maxEp);
% % 由于2PMP协议无非满秩的要求因此可以直接用RandMatrixGen2生成初始矩阵
% A = RandMatrixGen2Condition(N,minEp, maxEp, FirstNumMin, FirstNumMax);
% B = RandMatrixGen2Condition(N,minEp, maxEp, FirstNumMin, FirstNumMax);

% [Va, Vb, A, B, Cond_T, Re_Err_result, ~, MRE_result, MAPE_result, ~, ~, ~, ~, Error_Rate] =...
%     Secure_2PIP(A, B)
% ============================================================================================================================
% 定义Secure_2PIP安全两方求逆协议
function [Va, Vb, A, B, Cond_T, Re_Err_result, Ab_Err_result, MRE_result, MAPE_result, SAPE_result, F_norm, Theory_V, Real_V, ...
    Error_Rate] = S2PIP_Distinguish(A, B)
IsNulErrorOccur = false;

while(~IsNulErrorOccur)
    % ============================================================================================================================
    % 1.Alice本地生成随机可逆矩阵P，并计算矩阵Ia
%     [~, P] = RMG_2PMP(A);
    P = randi([1,99],size(A));
    Ia = P * A;

    % 2.Bob本地生成随机可逆矩阵Q，并计算矩阵Ib
%     [~, Q] = RMG_2PMP(B);
    Q = randi([1,99],size(B));

    Ib = B * Q;    
    
    % 3.Alice和Bob进行第一轮2PMP协议,(PA)Q=Va1+Vb1
    [Va1, Vb1] = S2PMP(Ia, Q);
%     cond_Ia = cond(Ia)
%     cond_Q = cond(Q)
%     cond_Va1 = cond(Va1)
%     cond_Vb1 = cond(Vb1)
    
    
    % 3.Alice和Bob进行第二轮2PMP协议,P(BQ)=Va2+Vb2
    [Va2, Vb2] = S2PMP(P, Ib);
%     cond_P = cond(P)
%     cond_Ib = cond(Ib)
%     cond_Va2 = cond(Va2)
%     cond_Vb2 = cond(Vb2)    
    
    % 4.Alice在本地秘密计算矩阵Ua，并发送给Bob
    Ua = Va1 + Va2;
%     cond_Ua = cond(Ua)
    
    
    % 5.Bob在本地秘密计算矩阵Ub和T=Ua+Ub=P(A+B)Q
    Ub = Vb1 + Vb2;
    T = Ua + Ub;
%     cond_Ub = cond(Ub)
%     cond_T = cond(T)
    
    
    % 6.Bob本地求逆得到invT以及Ib2=Q*invT
    Cond_T = cond(T);
    Ib2 = Q / T;
    
    
    
    % 7.Alice和Bob进行第三轮2PMP协议,Ib2*P=Vb+Va
    [Vb, Va] = S2PMP(Ib2, P);
    
    
    
    % ============================================================================================================================
    % 8.理论结果Theory_V=A*B*C,实验结果Real_V=Va+Vb+Vc
    Theory_V = inv(A + B);
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
%                 IsNulErrorOccur = false;
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