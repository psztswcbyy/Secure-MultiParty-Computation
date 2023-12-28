% ============================================================================================================================
% S3PMP是安全三方矩阵乘法协议的接口形式，涉及到的输入参数包含第一参与方Alice的初始矩阵A，第二参与方Bob的初始矩阵B，第三参与方
% Carol的初始矩阵C；涉及到的输出参数包含第一参与方Alice的初始矩阵A，第二参与方Bob的初始矩阵B，第三参与方Carol的初始矩阵C，三
% 方矩阵乘法结果的相对误差矩阵Re_Err_result，三方矩阵乘法结果的绝对误差矩阵Ab_Err_result，乘法结果的最大相对误差结果MRE_result，
% 乘法结果的均值绝对百分比误差MAPE_result，乘法结果的绝对百分比误差总和SAPE_result，乘法结果的F范数相对误差F_norm，理论的乘法
% 结果Theory_V，实际的乘法结果Real_V，以及最新更新版本的参数Diff表示原始数据矩阵A和混淆后的矩阵A_hat之间的阶数差值最大值。
% ============================================================================================================================
% clc;
% clear;
% format longE
% % ============================================================================================================================
% % 以下是对于方阵的测试
% N = 10;
% minEp = -8;
% maxEp = -minEp;
% FirstNumMin = 1;
% FirstNumMax = 1;
% %
% % A = RandInputMatrix(N,FirstNumMin,FirstNumMax,minEp,maxEp);
% % B = RandInputMatrix(N,FirstNumMin,FirstNumMax,minEp,maxEp);
% % C = RandInputMatrix(N,FirstNumMin,FirstNumMax,minEp,maxEp);
% 
% A = RandMatrixGen3(N,minEp, maxEp, FirstNumMin, FirstNumMax);
% B = RandMatrixGen3(N,minEp, maxEp, FirstNumMin, FirstNumMax);
% C = RandMatrixGen3(N,minEp, maxEp, FirstNumMin, FirstNumMax);
% 
% [~, ~, ~, ~, ~, ~, Diff_A, Diff_B, Diff_C, ~, ~, MRE_result, MAPE_result,...
%     ~, ~, ~, ~,Error_Rate] = S_3PMP(A, B, C)


% ============================================================================================================================
% 定义Secure_3PMP安全三方乘法协议
function [Va, Vb, Vc, A, B, C, Diff_A, Diff_B, Diff_C, Re_Err_result, Ab_Err_result, MRE_result, MAPE_result,...
    SAPE_result, F_norm, Theory_V, Real_V, Error_Rate] = S3PMP(A, B, C)

IsNulErrorOccur = false;

while(~IsNulErrorOccur)
    % ============================================================================================================================
    % 1.生成原始数据矩阵A、B、C，随机矩阵RA、RB、Rc，非满秩混淆矩阵A_hat、B_hat、C_hat
    [Ra, A_hat, Diff_A] = RMG_Enhanced(A);
    [Rb, B_hat, Diff_B] = RMG_Enhanced(B);
    [Rc, C_hat, Diff_C] = RMG_Enhanced(C);
    % % 对照组实验1，RMG_Unconstrained是指数完全对齐,但是不考虑矩阵非满秩的情况；
    % [Ra, A_hat] = RMG_Unconstrained(A);
    % [Rb, B_hat] = RMG_Unconstrained(B);
    % [Rc, C_hat] = RMG_Unconstrained(C);
    % Diff_A = 0;
    % Diff_B = 0;
    % Diff_C = 0;
    % % 对照组实验2，RMG_3PMP是指数指数只有一个没对齐但是对应矩阵时非满秩的情况；
    % [Ra, A_hat] = RMG_3PMP(A);
    % [Rb, B_hat] = RMG_3PMP(B);
    % [Rc, C_hat] = RMG_3PMP(C);
    % Diff_A = 0;
    % Diff_B = 0;
    % Diff_C = floor(log10(abs(C_hat(1,end)-C(1,end))));
    
    % ============================================================================================================================
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
%     Cond_B1 = cond(B1)
%     Cond_B2 = cond(B2)
%     Cond_B = cond(B)
    % ============================================================================================================================
    % 6.Alice生成随机矩阵ra，Va，计算Ta，t1发给BOb
    [~, ra] = RMG_Unconstrained(Ra * Rb * Rc);
    [~, Va] = RMG_Unconstrained(ra);
    Ta = Ma + Sa - Va - ra;
    t1 = Ra * B1;
    % ============================================================================================================================
    % 7.Carol计算t2发给Bob
    t2 = B2 * Rc;
    % ============================================================================================================================
    % 8.Bob生成随机矩阵rb，Vb，计算Mb、Sb发给Carol
    Sb = t1 * t2;
    [~, rb] = RMG_Unconstrained(Ra * Rb * Rc);
    rc = Ra * Rb * Rc - ra - rb;
    [~, Vb] = RMG_Unconstrained(rb);
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


