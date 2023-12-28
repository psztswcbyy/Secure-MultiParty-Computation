% =================================================================================
% RMG_2PMP函数用于复合多方计算对2PMP协议的调用；其中输出参数包含M
% 即转化前的原始矩阵，RM用于混淆的随机矩阵，M_hat转化后的随机矩阵；
% =================================================================================
clc;
clear;
format longE
N = 10;
minEp = -10;
maxEp = -minEp;
FirstNumMin = 1;
FirstNumMax = 1;

% A = RandInputMatrix(N,FirstNumMin,FirstNumMax,minEp,maxEp);
% 由于2PMP协议无非满秩的要求因此可以直接用RandMatrixGen2生成初始矩阵
A = RandMatrixGen2(N,minEp, maxEp, FirstNumMin, FirstNumMax);
[Ra, A_hat] = RFHMGen(A)

function [RM, M_hat] = RFHMGen(M)
% =================================================================================
% 随机矩阵RM和M_hat生成模块
% =================================================================================
% 阶数对齐，其中Ali_M为M_hat的阶数对齐后的矩阵
Check=false;
Num = 0;
while Check == false
    Num = Num + 1;
    [m, n] = size(M);
    % 由于原始矩阵中可能存在0元素，所以需要进行阶数矩阵特殊处理，以输入矩阵中非
    % 零元素最小值为参考标准，令阶数对其矩阵中的零元素对应阶数为差16阶的元素。
    M_new = sort(abs(M(:)));
    M_SecondMin = M_new(find(M_new>min(M_new),1));
    Ali_M = zeros(m,n);
    for i=1:m
        for j=1:n
            
            if M(i,j) ~= 0
                Ali_M(i,j) = floor(log2(abs(M(i,j)))) + randi([-1,1]);
            else
                Ali_M(i,j) = floor(log2(M_SecondMin)-54);
            end
        end
    end
    % Ori_M_hat为M_hat的原始尾数矩阵
    Ori_M_hat = rand(m,n,'double') + randi([1,1],m,n);
    % 生成混淆后的随机矩阵M_hat,由尾数矩阵和阶数矩阵构成
    M_hat = sign(M).* Ori_M_hat.* (2.^(Ali_M));
    if Num <= 10^4
        if cond(M_hat)>10^3
            Check = false;
        else
            fprintf('Find this matrix\n')
            RM = M_hat - M;
            break;
        end
    else
        fprintf('No matrix satisfy\n')
        RM = M_hat - M;
        break
    end
    % 生成随机矩阵RM
    % =================================================================================
    % % 无阶数对齐，直接生成随机矩阵时（对比实验方案）
    % M_hat = rand(N,N,'double')+ones(N);
    % RM = M_hat - M;
end
end