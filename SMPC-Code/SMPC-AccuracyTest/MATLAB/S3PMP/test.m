clc;
clear;
% N=5;
% r=N-1;
minEp = -8;
maxEp = 0;
% 
% Ori_M_hat = rand(r,r,'double') + randi([1,1],r,r);
% Exp_M_hat = 10.^(randi([minEp, maxEp],r,r));
% M_hat_base = Ori_M_hat.*Exp_M_hat;
% 
% Ori_Transformer_left = rand(N,r,'double') + randi([1,1],N,r);
% Exp_Transformer_left = 10.^(randi([0, 0],N,r));
% Transformer_left = Ori_Transformer_left.*Exp_Transformer_left;
% Ori_Transformer_right = rand(r,N,'double') + randi([1,1],r,N);
% Exp_Transformer_right = 10.^(randi([0, 0],r,N));
% Transformer_right = Ori_Transformer_right.*Exp_Transformer_right;
% M_hat = Transformer_left*M_hat_base*Transformer_right;
% 
% Rank_MHB = rank(M_hat_base)
% Rank_TL = rank(Transformer_left)
% Rank_TR = rank(Transformer_right)
% Rank_M_hat = rank(M_hat_base)

% P = randi([1,9],5,4);
% Q = randi([1,9],4,5);
% M = randi([1,9],4,4);
% M_hat = P*M*Q;
% rank(M_hat)
% rank(randi([1,1000],1)*M_hat+randi([1,1000],1))

n = 10^(maxEp+1);
m = 10^minEp;
n-m
delta = zeros(1,5);
index = randperm(5,4);
delta(index) = randi([-100,-1],1,4);
P = rand(5);
M = P*diag(delta)*inv(P);
a = min(min(M))
b = max(max(M))
if (max(delta) < (m*b-a*n)/(b-a))
    Standard = true;
else
    Standard = false;
end
e = delta(index(randi(4)))
k = (m-e)/a
M_hat = M*k+e*eye(5);
% =========================================================================================
% 该测试方法失败原因分析：
% =========================================================================================
% 动机是想通过对原始矩阵的数值区间，进行一个线性放缩，映射到指定的阶数范围区间内，但是实际
% 情况是这个线性映射中的比例关系可以保证在变换过程中，矩阵的秩不变，但是截距对应的平移关系
% 却无法保证秩的不变性，由于Y=k*X+σ*I中，σ就是截距，只有当σ等于矩阵X的负特征值时，才能
% 保证矩阵的秩不会变成满秩，该条件非常难以满足，因为大部分情况都是σ不等于特征值的，相加之
% 后往往得到的都是一个满秩矩阵，因此根据一个随机的非满秩矩阵试图将其均匀映射到一个精度区间
% 是一件非常困难的事情，该方案暂时搁置后续有空再研究






