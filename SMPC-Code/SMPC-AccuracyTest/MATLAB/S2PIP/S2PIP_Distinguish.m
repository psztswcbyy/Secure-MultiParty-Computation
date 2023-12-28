% ==================================================================================================================
% ������˵����
%  2PIPЭ�����ڰ�ȫ������������У����ﰴ��Э������ԭ�����漰�������������P��Q������Ȩ���з���ġ���Э�����漰
%  �������P��Q��Ϊ���������ȵ�����������P����Alice�;���Q����Bob��
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
% % ����2PMPЭ���޷����ȵ�Ҫ����˿���ֱ����RandMatrixGen2���ɳ�ʼ����
% A = RandMatrixGen2Condition(N,minEp, maxEp, FirstNumMin, FirstNumMax);
% B = RandMatrixGen2Condition(N,minEp, maxEp, FirstNumMin, FirstNumMax);

% [Va, Vb, A, B, Cond_T, Re_Err_result, ~, MRE_result, MAPE_result, ~, ~, ~, ~, Error_Rate] =...
%     Secure_2PIP(A, B)
% ============================================================================================================================
% ����Secure_2PIP��ȫ��������Э��
function [Va, Vb, A, B, Cond_T, Re_Err_result, Ab_Err_result, MRE_result, MAPE_result, SAPE_result, F_norm, Theory_V, Real_V, ...
    Error_Rate] = S2PIP_Distinguish(A, B)
IsNulErrorOccur = false;

while(~IsNulErrorOccur)
    % ============================================================================================================================
    % 1.Alice������������������P�����������Ia
%     [~, P] = RMG_2PMP(A);
    P = randi([1,99],size(A));
    Ia = P * A;

    % 2.Bob������������������Q�����������Ib
%     [~, Q] = RMG_2PMP(B);
    Q = randi([1,99],size(B));

    Ib = B * Q;    
    
    % 3.Alice��Bob���е�һ��2PMPЭ��,(PA)Q=Va1+Vb1
    [Va1, Vb1] = S2PMP(Ia, Q);
%     cond_Ia = cond(Ia)
%     cond_Q = cond(Q)
%     cond_Va1 = cond(Va1)
%     cond_Vb1 = cond(Vb1)
    
    
    % 3.Alice��Bob���еڶ���2PMPЭ��,P(BQ)=Va2+Vb2
    [Va2, Vb2] = S2PMP(P, Ib);
%     cond_P = cond(P)
%     cond_Ib = cond(Ib)
%     cond_Va2 = cond(Va2)
%     cond_Vb2 = cond(Vb2)    
    
    % 4.Alice�ڱ������ܼ������Ua�������͸�Bob
    Ua = Va1 + Va2;
%     cond_Ua = cond(Ua)
    
    
    % 5.Bob�ڱ������ܼ������Ub��T=Ua+Ub=P(A+B)Q
    Ub = Vb1 + Vb2;
    T = Ua + Ub;
%     cond_Ub = cond(Ub)
%     cond_T = cond(T)
    
    
    % 6.Bob��������õ�invT�Լ�Ib2=Q*invT
    Cond_T = cond(T);
    Ib2 = Q / T;
    
    
    
    % 7.Alice��Bob���е�����2PMPЭ��,Ib2*P=Vb+Va
    [Vb, Va] = S2PMP(Ib2, P);
    
    
    
    % ============================================================================================================================
    % 8.���۽��Theory_V=A*B*C,ʵ����Real_V=Va+Vb+Vc
    Theory_V = inv(A + B);
    Real_V = Va + Vb;
    
    % ============================================================================================================================
    % �趨�������ݾ��ȼ���eps(class(A)),���ڵ��ھ��ȼ��޵����ݣ��������ò���Ϊ0Ԫ��
    % �Ż����۽������Theory_V
    [m, n] = size(Theory_V);
    for i = 1:m
        for j = 1:n
            if abs(Theory_V(i, j))< max(m,n)*eps(class(Theory_V))*norm(Theory_V, 'inf')
                Theory_V(i, j) = 0;
            end
        end
    end
    % �Ż�ʵ�ʽ������Real_V
    [m, n] = size(Real_V);
    for i = 1:m
        for j = 1:n
            if abs(Real_V(i, j))< max(m,n)*eps(class(Real_V))*norm(Real_V, 'inf')
                Real_V(i, j) = 0;
            end
        end
    end
    
    
    % ����������Ab_Err_result����ά��
    Ab_Err_result = Real_V - Theory_V;
    [m, n] = size(Ab_Err_result);
    
    % ��ʼ�����������Re_Err_result
    Re_Err_result = zeros(m, n);
    
    % disp(Real_V),disp(Theory_V),disp(Ab_Err_result);
    
    % �������۾����ж�ӦԪ��Ϊ0���������������һ�ֵ�����Ԫ��Ϊ0���Ҿ����������ȼ��ޣ�Ĭ��Ϊʵ��Ԫ��Ϊ0
    % �ڶ������������Ԫ�ط�0������ʵ�����ȥ���������
    % ���������������Ԫ��Ϊ0�����Ǿ�������ھ��ȼ����ڣ���ֱ����Ԫ�ص���������ڸþ������
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
    
    % ���������MRE
    MRE_result = max(Re_Err_result(:));
    % ���԰ٷֱ�����SAPE
    SAPE_result = sum(abs(Re_Err_result(:)));
    % ��ֵ���԰ٷֱ����MAPE
    MAPE_result = sum(abs(Re_Err_result(:)))/numel(Re_Err_result);
    
    % F�������,Fa����������ķ���
    F_absolute = norm(Ab_Err_result, 'fro');
    F_theory = norm(Theory_V, 'fro');
    F_norm = F_absolute/F_theory;
    
    % ������ErrorRateͳ�ƣ�����������Ż������޷��ɹ�ʵ�ֵ�ʱ�򣬸��ô�����ָ��    
    Error_Rate =  numel(find(Re_Err_result >= 0.5*10^-3)) / numel(Re_Err_result);
    
end
end