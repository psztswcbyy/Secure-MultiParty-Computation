% =================================================================================
% RandBigMatrixGen2函数用于将原始数据矩阵转化为维度扩大的随机矩阵Big_M(其结
% 构为一个分块对角矩阵，其中对角线上的块为原始数据矩阵的拷贝，其余位置为0)涉
% 及到的输入参数包含原始矩阵的维度N，矩阵元素数值精度的阶数最大值maxEp和最小
% 值minEp、首位数字的最大值FirstNumMin和最小值FirstNumMax；输入数据元素的最大
% 有效数字位数S-bit；输出参数包含M即转化前的原始矩阵，RM用于混淆的随机矩阵，
% M_hat转化后的随机矩阵；
% 其最明显效果是可以实现一个计算结果正确性自检测
% 该方法后期需要完善和优化成一个同时满足非满秩和精度对齐的函数，思维导图中记
% 录具体方法
% =================================================================================
% 测试专用参数
% =================================================================================
% clc;clear;
% N=3;
% minEp=-5;
% maxEp=0;
% FirstNumMin=1;
% FirstNumMax=1;
% % 创建一个mxn且秩为r的非满秩随机矩阵
% % A = (randi([1,10],m,r)*rand(r,r)*randi([1,10],r,n)).*10.^(randi([minEp, maxEp],m,n));
%
% [A, Ra, A_hat] = RandBigMatrixGen(N,minEp,maxEp,FirstNumMin,FirstNumMax)

function [Big_M, RM, M_hat] = RandBigMatrixGen2(N, minEp, maxEp, FirstNumMin, FirstNumMax)
% =================================================================================
% 该模块负责：尾数矩阵Ori_M生成，阶数矩阵Exp_M生成，eg：3.14E-12，其中3.14
% 代表尾数项，E-12代表阶数项
% =================================================================================
% 原始矩阵M生成模块
% =================================================================================
% 直接生成15位有效数字作为小数尾数的思路
% Ori_M = (2*rand(N,N,'double')-1) + randi([FirstNumMin,FirstNumMax],N,N);
Ori_M = rand(N,N,'double') + randi([FirstNumMin,FirstNumMax],N,N);

Exp_M = 10.^(randi([minEp, maxEp],N,N));
M = Ori_M.*Exp_M;

% =================================================================================
% 扩展矩阵Big_M生成模块
% =================================================================================
% 直接生成15位有效数字作为小数尾数的思路
Big_M = zeros(2 * N);
Big_M(1:N, 1:N) = M;
Big_M(N+1:2*N, N+1:2*N) = M;

% =================================================================================
% 随机矩阵RM和M_hat生成模块
% =================================================================================
% 阶数对齐，Ali_M为M_hat的阶数对齐后的矩阵,OAM是阶数对齐矩阵OrderAlignMatrix简称
Ali_M = zeros(2*N,2*N,'double');
OAM = floor(log2(abs(M))) + randi([-1,1],N,N);
Ali_M(1:N,1:N) = OAM;
Ali_M(N+1:2*N,1:N) = OAM;
Ali_M(1:N,N+1:2*N) = OAM;
Ali_M(N+1:2*N,N+1:2*N) = OAM;
% Ori_M_hat为M_hat的原始尾数矩阵
Ori_M_hat = rand(2*N,2*N,'double') + randi([FirstNumMin,FirstNumMax],2*N,2*N);
% 生成混淆后的随机矩阵M_hat,由尾数矩阵和阶数矩阵构成
M_hat = Ori_M_hat.*(2.^(Ali_M));
% 生成随机矩阵RM
RM = M_hat - Big_M;

% =================================================================================
% % 无阶数对齐，直接生成随机矩阵时（对比实验方案）
% M_hat = rand(N,N,'double')+randi([FirstNumMin,FirstNumMax],N,N);
% RM = M_hat - M;
end