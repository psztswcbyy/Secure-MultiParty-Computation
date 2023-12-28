% ============================================================================================================================
% Secure_2PMP��������ʵ�ְ�ȫ��������˷�Э�飬�漰����������������������ķ���ά��N�����뷽���Ԫ�ؽ�����СֵminEp�����뷽
% ���Ԫ�ؽ������ֵmaxEp�����뷽��Ԫ�ظ�λ����СֵFirstNumMin�����뷽��Ԫ�ظ�λ�����ֵFirstNumMax���漰�����������������һ
% ���뷽Alice�ĳ�ʼ����A���ڶ����뷽Bob�ĳ�ʼ����B����������˷���������������Re_Err_result����������˷�����ľ�������
% ��Ab_Err_result���˷�����������������MRE_result���˷�����ľ�ֵ���԰ٷֱ����MAPE_result���˷�����ľ��԰ٷֱ�����ܺ�
% SAPE_result,�˷������F����������F_norm�����۵ĳ˷����Theory_V��ʵ�ʵĳ˷����Real_V
% ============================================================================================================================
% clc;
% clear;
% format longE
% N = 5;
% minEp = -10;
% maxEp = 10;
% FirstNumMin = 1;
% FirstNumMax = 1;
% [A, B, Re_Err_result, Ab_Err_result, MRE_result, MAPE_result, SAPE_result, F_norm, Theory_V, Real_V] =...
%     Secure2PMP(N, minEp, maxEp, FirstNumMin, FirstNumMax)
% ============================================================================================================================
% ����Secure_2PMP��ȫ�����˷�Э��
function [A, B, Re_Err_result, Ab_Err_result, MRE_result, MAPE_result, SAPE_result, F_norm, Theory_V, Real_V] =...
    Secure_2PMP(N, minEp, maxEp, FirstNumMin, FirstNumMax)
IsNulErrorOccur = false;

while(~IsNulErrorOccur)
    % ============================================================================================================================
    % ����ԭʼ���ݾ���A��B���������RA��RB����������A_hat��B_hat
    % ============================================================================================================================
    [A, Ra, A_hat] = RandMatrixGen2(N, minEp, maxEp, FirstNumMin, FirstNumMax);
    [B, Rb, B_hat] = RandMatrixGen2(N, minEp, maxEp, FirstNumMin, FirstNumMax);
    
    % ============================================================================================================================
    % ������չ�ɴ��ģ������������˼·������ֱ��ͨ���÷�ʽ���ɷ����Ⱦ��󣬲��߱�˫���ʼ�������
    % ============================================================================================================================
    % [A, Ra, A_hat] = RandBigMatrixGen2(N, minEp, maxEp, FirstNumMin, FirstNumMax);
    % [B, Rb, B_hat] = RandBigMatrixGen2(N, minEp, maxEp, FirstNumMin, FirstNumMax);
    % [~, ra, Vb] = RandBigMatrixGen2(N, minEp, maxEp, FirstNumMin, FirstNumMax);
    
    % ============================================================================================================================
    % �����������ra��rb�Լ�Bob��˽���������Vb
    % [ra, Vb] = RandMatrixGen2(N, minEp, maxEp, FirstNumMin, FirstNumMax);%��仰����������һ���ؼ����������������봦����Լ��󽵵����
    [~, ra] = RMG_2PMP(Ra*Rb);%�������RMG_2PMPЭ�������������Ч���ǳ���
    [~, Vb]  =RMG_2PMP(ra);
    rb = Ra * Rb - ra;
    
    % Bob����T����
    T = A_hat * B + (rb - Vb);
    
    
    % �������Va
    Va = T + ra - (Ra * B_hat);
    
    
    % ���۽��Theory_V=A*B,ʵ����Real_V=Va+Vb
    Theory_V = A * B;
    Real_V = Va + Vb;
    
    % ============================================================================================================================
    % % �������
    % OAM_T = floor(log2(abs(T))) + randi([-1,1],2*N,2*N);
    % T_bad = (rand(2*N)+1).*(2.^OAM_T)+T;
    % Va_bad = T_bad + ra - (Ra * B_hat);
    % Real_V_bad = Va_bad + Vb;
    % Theory_V
    % Real_V_bad
    
    % ============================================================================================================================
    % ���󾫶��Ż�ģ��
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
                IsNulErrorOccur = false;
            end
        end
    end
    
    % ============================================================================================================================
    % ����������ģ��
    % ============================================================================================================================
    
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
    
end
end

