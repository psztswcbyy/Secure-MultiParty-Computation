% =================================================================================
% 【功能说明】
%  RMG_RMG_Unconstrained函数用于复合多方计算对3PMP协议的调用；需注意这里是不完全版
%  本，因为没有考虑混淆矩阵非满秩的条件，只是满足阶数对齐的情况，因此只适合
%  于半诚实的计算环境，其中输入参数M为转化前的原始矩阵，RM是用于混淆的随机矩
%  阵，M_hat是混淆后的随机矩阵；
% =================================================================================

function [RM, M_hat] = RMG_Unconstrained(M)
% =================================================================================
% 随机矩阵RM和M_hat生成模块
% =================================================================================
% 阶数对齐，其中Ali_M为M_hat的w阶数对齐后的矩阵
N = size(M);
DynamicRange = 0;
% 由于原始矩阵中可能存在0元素，所以需要进行阶数矩阵特殊处理，以输入矩阵中非
% 零元素最小值为参考标准，令阶数对其矩阵中的零元素对应阶数为差16阶的元素。
M_new = sort(abs(M(:)));
M_SecondMin = M_new(find(M_new>min(M_new),1));
Ali_M = zeros(N);
for i=1:N
    for j=1:N
        
        if M(i,j) ~= 0
            Ali_M(i,j) = floor(log2(abs(M(i,j)))) + randi([-DynamicRange,DynamicRange]);
        else
            Ali_M(i,j) = floor(log2(M_SecondMin)-52);
        end
    end
end
% Ori_M_hat为M_hat的原始尾数矩阵
Ori_M_hat = rand(N,'double') + randi([1,1],N);
% 生成混淆后的随机矩阵M_hat,由尾数矩阵和阶数矩阵构成
M_hat = sign(M).*Ori_M_hat.*(2.^(Ali_M));
% 生成随机矩阵RM
RM = M_hat - M;
end