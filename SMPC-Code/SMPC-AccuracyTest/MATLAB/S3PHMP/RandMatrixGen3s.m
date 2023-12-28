% =================================================================================
% 【函数功能说明】
%  本用例是对RandMatrixGen3s函数在3PMP协议非满秩但无阶数对齐环境下的测试，其
%  中输入包含原始矩阵M以及最后右上角的精度损失阈值EPS_Threshold，输出参数包
%  含：RM用于混淆的随机矩阵，M_hat转化后的随机矩阵；
% 【阶数对齐实验测试方案说明】
%  通过构建四个子模块，寻求最优的基向量组，这个基向量组通过特定的线性组合可以、
%  构成指定精度分布的非满秩矩阵，每个模块之间满足相应的关系：
%  1.首先由原始矩阵M的阶数分布生成对应需求的阶数对齐模块Core_M;
%  2.然后根据模块M1和模块M2之间存在的线性关系解析出及向量组（b1,b2...bn）的数
%    值解；
%  3.最后根据基底和模块M3共同生成模块M4对应的奇异值（该元素的精度此时并未实现
%    对齐），至此出去模块M4以外的三个模块均已实现阶数的对齐；
%  4.此时可以采用一个阈值作为循环涮选出符合要求的阶数对齐矩阵，由于只有一个元
%    素对应M4右上角的位置，因此很容易根据这个元素找到所有元素阶数均对齐的矩阵
% =================================================================================

function [M, RM, M_hat] = RandMatrixGen3s(N,minEp, maxEp, FirstNumMin, FirstNumMax)
% =================================================================================
% 生成原始输入矩阵M
% =================================================================================
Ori_M = 2*rand(N,'double')-1 ;%step1:生成正态[-1,1]区间分布；
Ori_M = sign(Ori_M).*randi([FirstNumMin, FirstNumMax], N, N) + Ori_M;%step2：符号函数对应加成，首位数为1；
Exp_M = 10.^(randi([minEp, maxEp],N,N));%生成原始阶数矩阵
M = Ori_M.*(Exp_M);

% =================================================================================
% 对输入矩阵做一个预处理，把其中的零元素先变成一个相对很小的非零数值
% =================================================================================
N = size(M,1);
M_Pre = M;
M_new = sort(abs(M(:)));
M_SecondMin = M_new(find(M_new>min(M_new),1));
if any(M_Pre(:) == 0)
    [row, col] = find(M_Pre == 0);
    for i = 1:numel(row)
        M_Pre(row(i), col(i)) = floor(log2(M_SecondMin)-52);
    end
end

% =================================================================================
% 阶数不对齐的情况
Ali_M = randi([minEp, maxEp],N,N);
Rand_M = sign(M_Pre).*(rand(N,'double')+randi([1,9],N));

% 目标随机混淆矩阵Core_M
Core_M = Rand_M.*10.^Ali_M;
% 模块1 <==>  Block_M1
Block_M1 = Core_M(1,1:N-1);
% 模块2 <==>  Block_M2\
Block_M2 = Core_M(2:end,1:N-1);
% 模块3 <==>  Block_M3
Block_M3 = Core_M(2:end,end);
% 基向量首位元素分布 <==>  Vector_X
[P,R,C] = equilibrate(Block_M2');
tol = 1e-11;
maxit = size(Block_M2, 1);
newBlock_M2 = R * P * Block_M2' * C;
newBlock_M1 = R * P * Block_M1';
[y,~,~,~,~] = gmres(newBlock_M2,newBlock_M1,[],tol,maxit);
Vector_X = C * y;
% Vector_X = Block_M2'\Block_M1';
% 模块4 <==>  Block_M4
Block_M4 = Block_M3'*Vector_X;
% 更新最终的M_hat
M_hat = [Block_M1 Block_M4; Block_M2 Block_M3];
% 生成随机RM
RM = M_hat-M;
end



