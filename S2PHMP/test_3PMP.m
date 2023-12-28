% ============================================================================================================================
% Secure_3PMP��������ʵ�ְ�ȫ��������˷�Э�飬�漰����������������������ķ���ά��N�����뷽���Ԫ�ؽ�����СֵminEp�����뷽
% ���Ԫ�ؽ������ֵmaxEp�����뷽��Ԫ�ظ�λ����СֵFirstNumMin�����뷽��Ԫ�ظ�λ�����ֵFirstNumMax���漰�����������������һ
% ���뷽Alice�ĳ�ʼ����A���ڶ����뷽Bob�ĳ�ʼ����B���������뷽Carol�ĳ�ʼ����C����������˷���������������Re_Err_result��
% ��������˷�����ľ���������Ab_Err_result���˷�����������������MRE_result���˷�����ľ�ֵ���԰ٷֱ����MAPE_result��
% �˷�����ľ��԰ٷֱ�����ܺ�SAPE_result���˷������F����������F_norm�����۵ĳ˷����Theory_V��ʵ�ʵĳ˷����Real_V
% ============================================================================================================================
% clc;
% clear;
% format longE
% N = 5;
% minEp = -6;
% maxEp = 6;
% FirstNumMin = 1;
% FirstNumMax = 1;
% 
% [A, B, C, Re_Err_result, Ab_Err_result, MRE_result, MAPE_result, SAPE_result, F_norm, Theory_V, Real_V] =...
%     Secure3PMP(N, minEp, maxEp, FirstNumMin, FirstNumMax)

% % ============================================================================================================================
% clc;
% clear;
% format longE
% N = 5;
% minEp = -9;
% maxEp = 9;
% FirstNumMin = 1;
% FirstNumMax = 1;
% 
% A = RandInputMatrix(N,FirstNumMin,FirstNumMax,minEp,maxEp);
% B = RandInputMatrix(N,FirstNumMin,FirstNumMax,minEp,maxEp);
% C = RandInputMatrix(N,FirstNumMin,FirstNumMax,minEp,maxEp);
% [~, ~, ~, ~, ~, ~, ~, ~, ~, Theory_V, Real_V] =...
%     S_3PMP(A, B, C)
% % ============================================================================================================================
% % ����Secure_3PMP��ȫ�����˷�Э��
% function [A, B, C, Re_Err_result, Ab_Err_result, MRE_result, MAPE_result, SAPE_result, F_norm, Theory_V, Real_V] =...
%     S_3PMP(A, B, C)
% % ============================================================================================================================
% % 1.����ԭʼ���ݾ���A��B��C���������RA��RB��Rc�������Ȼ�������A_hat��B_hat��C_hat
% [Ra, A_hat] = RMG_Enhanced(A);
% [Rb, B_hat] = RMG_Enhanced(B);
% [Rc, C_hat] = RMG_Enhanced(C);

% ============================================================================================================================
% ����Secure_3PMP��ȫ�����˷�Э��
function [A, B, C, Re_Err_result, Ab_Err_result, MRE_result, MAPE_result, SAPE_result, F_norm, Theory_V, Real_V] =...
    test_3PMP(A,B,C,Ra,Rb,Rc)

% ============================================================================================================================
% 1.����ԭʼ���ݾ���A��B��C���������RA��RB��Rc�������Ȼ�������A_hat��B_hat��C_hat
A_hat = A+Ra;
B_hat = B+Rb;
C_hat = C+Rc;

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
% ============================================================================================================================
% 6.Alice�����������ra��Va������Ta��t1����BOb
% [ra, ~, Va] = RandMatrixGen3p(N, minEp, maxEp, FirstNumMin, FirstNumMax);
[~, ra] = RMG_Unconstrained(Ra * Rb * Rc);
[~, Va] = RMG_Unconstrained(ra);
Ta = Ma + Sa - Va - ra;
t1 = Ra * B1;
% ============================================================================================================================
% 7.Carol����t2����Bob
t2 = B2 * Rc;
% ============================================================================================================================
% 8.Bob�����������rb��Vb������Mb��Sb����Carol
% [rb, ~, Vb] = RandMatrixGen3p(N, minEp, maxEp, FirstNumMin, FirstNumMax);
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

end

