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
% % ============================================================================================================================
% minEp = -10;
% maxEp = 10;
% FirstNumMin = 1;
% FirstNumMax = 1;
% % ============================================================================================================================
% % �����Ƕ�����������ά�ȵľ������
% DimA = randi([2,10],1,2);
% A_row = DimA(1);
% A_col = DimA(2);
% B_row = A_col;
% B_col = randi([2,10],1);
% % ============================================================================================================================
% % �����Ƕ��ڷ���Ĳ���
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
% ����Secure_2PMP��ȫ�����˷�Э��
function [Va, Vb, A, B, Re_Err_result, Ab_Err_result, MRE_result, MAPE_result, SAPE_result, F_norm, Theory_V, Real_V] =...
    S2PMP(A, B)
IsNulErrorOccur = false;

while(~IsNulErrorOccur)
    % ============================================================================================================================
    % �����������A��B�����������RA��RB����������A_hat��B_hat
    % ============================================================================================================================
    [Ra, A_hat] = RMG_2PMP(A);
    [Rb, B_hat] = RMG_2PMP(B);
    
    % ============================================================================================================================
    % �����������ra��rb�Լ�Bob��˽���������Vb
    [~, ra] = RMG_2PMP(Ra*Rb);
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
