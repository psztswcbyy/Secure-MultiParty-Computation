% clc;
% clear;
% % =========================================================================
% % 单位矩阵E_A拆分为两个行列互不相干的矩阵ErA1和ErA2
% % 测试方阵的矩阵分解
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
% % 测试非方阵的矩阵分解
% m=3;n=6;N=2;
% % A = randi([1,10],m,N)*randi([1,10],N,N)*randi([1,10],N,n); % 创建一个mxn且秩为N的非满秩随机矩阵
% A = rand(m,N)*rand(N,N)*rand(N,n); % 创建一个mxn且秩为N的非满秩随机矩阵
% k = rank(A);
% [U,V] = FRD(A,k)
% % =========================================================================
% % 测试参数
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

[U,S,V] = svd(A, 'econ'); % 进行奇异值分解
% [U,S,V] = svd(A)
% 选取非零奇异值
r = rank(A);
if r>k
    disp("Parameter k must be greater than " + r);
    return
else
    U = U(:,1:r);
    S = S(1:r,1:r);
    V = V(:,1:r);
    % 到这里B和C就是满秩分解的结果了
    B = U*S; % 这里的B是一个mxr的矩阵，其列是A的列空间的一组基
    C = V';  % 这里的C是一个rxn的矩阵，其行是A的行空间的一组基
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





