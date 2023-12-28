% clc;
% clear;
% % =========================================================================
% % ��λ����E_A���Ϊ�������л�����ɵľ���ErA1��ErA2
% % ���Է���ľ���ֽ�
% N=6;
% E = eye(N);
% Er = E;
% id_A = randperm(size(Er,2),2);
% Er(:,id_A(1)) = 0;
% Er(:,id_A(2)) = 0;
% P = randi(N,N);
% Q = randi(N,N);
% A = P\Er/Q;
% [U,V] = FRD(A,rank(A))
% % =========================================================================
% % ���ԷǷ���ľ���ֽ�
% m=3;n=6;N=2;
% % A = randi([1,10],m,N)*randi([1,10],N,N)*randi([1,10],N,n); % ����һ��mxn����ΪN�ķ������������
% A = rand(m,N)*rand(N,N)*rand(N,n); % ����һ��mxn����ΪN�ķ������������
% k = rank(A);
% [U,V] = FRD(A,k)
% % =========================================================================
% % ���Բ���
% Rank_U = rank(U)
% Rank_V = rank(V)
% Rank_A = rank(A)
% UV = U*V
% A
% RME = abs((U*V-A)./A)

% fullrank_decomposition()
function [Slice1, Slice2] = fullrank_decomposition(A,k)
m = size(A,1);
n = size(A,2);

[U,S,V] = svd(A, 'econ'); % ��������ֵ�ֽ�
% [U,S,V] = svd(A)
% ѡȡ��������ֵ
r = rank(A);
if r>k
    disp("Parameter k must be greater than " + r);
    return
else
    U = U(:,1:r);
    S = S(1:r,1:r);
    V = V(:,1:r);
    % ������B��C�������ȷֽ�Ľ����
    B = U*S; % �����B��һ��mxr�ľ���������A���пռ��һ���
    C = V';  % �����C��һ��rxn�ľ���������A���пռ��һ���
    Q = randn(k,k);
    Added = abs(k-r);
    NewB = [B zeros(m,Added)];
    NewC = [C ; zeros(Added,n)];
    if k ~= r
        Slice1 = NewB*Q;
        Slice2 = Q\NewC;
    else
        Slice1 = NewB;
        Slice2 = NewC;
    end
end
end





