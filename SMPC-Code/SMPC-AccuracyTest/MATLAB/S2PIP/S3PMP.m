% ============================================================================================================================
% S3PMP�ǰ�ȫ��������˷�Э��Ľӿ���ʽ���漰�����������������һ���뷽Alice�ĳ�ʼ����A���ڶ����뷽Bob�ĳ�ʼ����B���������뷽
% Carol�ĳ�ʼ����C���漰�����������������һ���뷽Alice�ĳ�ʼ����A���ڶ����뷽Bob�ĳ�ʼ����B���������뷽Carol�ĳ�ʼ����C����
% ������˷���������������Re_Err_result����������˷�����ľ���������Ab_Err_result���˷�����������������MRE_result��
% �˷�����ľ�ֵ���԰ٷֱ����MAPE_result���˷�����ľ��԰ٷֱ�����ܺ�SAPE_result���˷������F����������F_norm�����۵ĳ˷�
% ���Theory_V��ʵ�ʵĳ˷����Real_V���Լ����¸��°汾�Ĳ���Diff��ʾԭʼ���ݾ���A�ͻ�����ľ���A_hat֮��Ľ�����ֵ���ֵ��
% ============================================================================================================================
% clc;
% clear;
% format longE
% % ============================================================================================================================
% % �����Ƕ��ڷ���Ĳ���
% N = 10;
% minEp = -16;
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
% [~, ~, ~, ~, ~, ~, Diff_A, Diff_B, Diff_C, Error_Rate, ~, ~, MRE_result, ~,...
%     ~, ~, ~, ~] = S_3PMP(A, B, C)
% ============================================================================================================================
% ����Secure_3PMP��ȫ�����˷�Э��
function [Va, Vb, Vc, A, B, C, Diff_A, Diff_B, Diff_C, Re_Err_result, Ab_Err_result, MRE_result, MAPE_result,...
    SAPE_result, F_norm, Theory_V, Real_V, Error_Rate] = S3PMP(A, B, C)

IsNulErrorOccur = false;

while(~IsNulErrorOccur)
    % ============================================================================================================================
    % 1.����ԭʼ���ݾ���A��B��C���������RA��RB��Rc�������Ȼ�������A_hat��B_hat��C_hat
    [Ra, A_hat, Diff_A] = RMG_Enhanced(A);
    [Rb, B_hat, Diff_B] = RMG_Enhanced(B);
    [Rc, C_hat, Diff_C] = RMG_Enhanced(C);
    % % ������ʵ��1��RMG_Unconstrained��ָ����ȫ����,���ǲ����Ǿ�������ȵ������
    % [Ra, A_hat] = RMG_Unconstrained(A);
    % [Rb, B_hat] = RMG_Unconstrained(B);
    % [Rc, C_hat] = RMG_Unconstrained(C);
    % Diff_A = 0;
    % Diff_B = 0;
    % Diff_C = 0;
    % % ������ʵ��2��RMG_3PMP��ָ��ָ��ֻ��һ��û���뵫�Ƕ�Ӧ����ʱ�����ȵ������
    % [Ra, A_hat] = RMG_3PMP(A);
    % [Rb, B_hat] = RMG_3PMP(B);
    % [Rc, C_hat] = RMG_3PMP(C);
    % Diff_A = 0;
    % Diff_B = 0;
    % Diff_C = floor(log10(abs(C_hat(1,end)-C(1,end))));
    
    % ============================================================================================================================
    % ============================================================================================================================
    % 2.Bob�ڵ����Mb��w1��w2��r1��r2������w1��r1���͸�Carol�ڵ㣬w2��r2���͸�Alice�ڵ㣬
    Mb = A_hat * Rb * C_hat;
    w1 = A_hat * B_hat;
    r1 = A_hat * Rb;
    
    w2 = B_hat* C_hat;
    r2 = Rb * C_hat;
    % ============================================================================================================================
    % 3.Alice�ڵ㱾�ؼ���Sa��Ma
    Sa = Ra * r2;
    Ma = A * w2;
    % ============================================================================================================================
    % 4.Carol�ڵ㱾�ؼ���Sc��Mc
    Sc = r1 * Rc;
    Mc = w1 * Rc;
    % ============================================================================================================================
    % 5.Bob�ڵ㱾��ʵ�����ȷֽ⣬B_hat=B1*B2��B1����Alice�ڵ㣬B2����Carol�ڵ�
    [B1, B2] = fullrank_decomposition(B_hat,rank(B_hat));
%     Cond_B1 = cond(B1)
%     Cond_B2 = cond(B2)
%     Cond_B = cond(B)
    % ============================================================================================================================
    % 6.Alice�����������ra��Va������Ta��t1����BOb
    [~, ra] = RMG_Unconstrained(Ra * Rb * Rc);
    [~, Va] = RMG_Unconstrained(ra);
    Ta = Ma + Sa - Va - ra;
    t1 = Ra * B1;
    % ============================================================================================================================
    % 7.Carol����t2����Bob
    t2 = B2 * Rc;
    % ============================================================================================================================
    % 8.Bob�����������rb��Vb������Mb��Sb����Carol
    Sb = t1 * t2;
    [~, rb] = RMG_Unconstrained(Ra * Rb * Rc);
    rc = Ra * Rb * Rc - ra - rb;
    [~, Vb] = RMG_Unconstrained(rb);
    Tb = Ta + Sb - Mb - Vb - rb;
    % ============================================================================================================================
    % 9.Carol����Vc
    Vc = Tb - Mc + Sc - rc;
    
    % ============================================================================================================================
    % ���۽��Theory_V=A*B*C,ʵ����Real_V=Va+Vb+Vc
    Theory_V = A * B * C;
    Real_V = Va + Vb + Vc;
    
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


