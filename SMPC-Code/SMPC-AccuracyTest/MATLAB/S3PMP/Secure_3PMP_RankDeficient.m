% ============================================================================================================================
% Secure_3PMP_RankDeficient函数用于实现只满足非满秩约束但阶数不对齐版本的安全三方矩阵乘法协议，涉及到的输入参数包含参与
% 计算的方阵维度N，参与方阵的元素阶数最小值minEp，参与方阵的元素阶数最大值maxEp，参与方阵元素个位数最小值FirstNumMin，参
% 与方阵元素个位数最大值FirstNumMax；涉及到的输出参数包含第一参与方Alice的初始矩阵A，第二参与方Bob的初始矩阵B，第三参与
% 方Carol的初始矩阵C，三方矩阵乘法结果的相对误差矩阵Re_Err_result，三方矩阵乘法结果的绝对误差矩阵Ab_Err_result，乘法结果
% 的最大相对误差结果MRE_result，乘法结果的均值绝对百分比误差MAPE_result，乘法结果的绝对百分比误差总和SAPE_result，乘法结
% 果的F范数相对误差F_norm，理论的乘法结果Theory_V，实际的乘法结果Real_V
% ============================================================================================================================
% clc;
% clear;
% format longE
% N = 5;
% minEp = -8;
% maxEp = 8;
% FirstNumMin = 1;
% FirstNumMax = 1;
% 
% [A, B, C, Re_Err_result, Ab_Err_result, MRE_result, MAPE_result, ~, ~, ~, ~,Error_Rate] =... 
%     Secure_3PMPRankDeficient(N, minEp, maxEp, FirstNumMin, FirstNumMax)
% ============================================================================================================================
% 定义Secure_3PMP安全三方乘法协议
function [A, B, C, Re_Err_result, Ab_Err_result, MRE_result, MAPE_result, SAPE_result, F_norm, Theory_V, Real_V, Error_Rate] =...
    Secure_3PMP_RankDeficient(N, minEp, maxEp, FirstNumMin, FirstNumMax)

IsNulErrorOccur = false;

while(~IsNulErrorOccur)
    % ============================================================================================================================
    % 1.生成原始数据矩阵A、B、C，随机矩阵RA、RB、Rc，非满秩混淆矩阵A_hat、B_hat、C_hat
    [A, Ra, A_hat] = RandMatrixGen3s(N, minEp, maxEp, FirstNumMin, FirstNumMax);
    [B, Rb, B_hat] = RandMatrixGen3s(N, minEp, maxEp, FirstNumMin, FirstNumMax);
    [C, Rc, C_hat] = RandMatrixGen3s(N, minEp, maxEp, FirstNumMin, FirstNumMax);
    % ============================================================================================================================
    % 2.Bob节点计算Mb、w1、w2、r1、r2，并将w1、r1发送给Carol节点，w2、r2发送给Alice节点，
    Mb = A_hat * Rb * C_hat;
    w1 = A_hat * B_hat;
    r1 = A_hat * Rb;
    
    w2 = B_hat* C_hat;
    r2 = Rb * C_hat;
    % ============================================================================================================================
    % 3.Alice节点本地计算Sa，Ma
    Sa = Ra * r2;
    Ma = A * w2;
    % ============================================================================================================================
    % 4.Carol节点本地计算Sc，Mc
    Sc = r1 * Rc;
    Mc = w1 * Rc;
    % ============================================================================================================================
    % 5.Bob节点本地实现满秩分解，B_hat=B1*B2，B1发给Alice节点，B2发给Carol节点
    [B1, B2] = fullrank_decomposition(B_hat,rank(B_hat));
    % ============================================================================================================================
    % 6.Alice生成随机矩阵ra，Va，计算Ta，t1发给BOb
    [ra, ~, Va] = RandMatrixGen3s(N, minEp, maxEp, FirstNumMin, FirstNumMax);
%     [~, ra] = RMG_Unconstrained(Ra * Rb * Rc);
%     [~, Va] = RMG_Unconstrained(ra);
    Ta = Ma + Sa - Va - ra;
    t1 = Ra * B1;
    % ============================================================================================================================
    % 7.Carol计算t2发给Bob
    t2 = B2 * Rc;
    % ============================================================================================================================
    % 8.Bob生成随机矩阵rb，Vb，计算Mb、Sb发给Carol
    [rb, ~, Vb] = RandMatrixGen3s(N, minEp, maxEp, FirstNumMin, FirstNumMax);
    Sb = t1 * t2;
%     [~, rb] = RMG_Unconstrained(Ra * Rb * Rc);
    rc = Ra * Rb * Rc - ra - rb;
%     [~, Vb] = RMG_Unconstrained(rb);
    Tb = Ta + Sb - Mb - Vb - rb;
    % ============================================================================================================================
    % 9.Carol计算Vc
    Vc = Tb - Mc + Sc - rc;
    
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
    
    % 错误率ErrorRate统计，当上述误差优化部分无法成功实现的时候，改用错误率指标
    Error_Rate =  numel(find(Re_Err_result >= 0.5*10^-3)) / numel(Re_Err_result);
end
end

