% =================================================================================
% 【函数功能说明】RMG_3PMP_Enhanced函数用于复合函数对3PMP函数的调用（此方法
% 以右上角元素的精度损失阶数阈值EPS_Threshold作为筛选标准，以力求做到元素阶
% 数上的完全对齐）；其中输出参数RM用于混淆的随机矩阵，M_hat转化后的随机矩阵；
% 【阶数对齐实验测试方案说明】
% 通过构建四个子模块，寻求最优的基向量组，这个基向量组通过特定的线性组合可以、
% 构成指定精度分布的非满秩矩阵，每个模块之间满足相应的关系：
% 1.首先由原始矩阵M的阶数分布生成对应需求的阶数对齐模块Core_M;
% 2.然后根据模块M1和模块M2之间存在的线性关系解析出及向量组（b1,b2...bn）的数
%   值解；
% 3.最后根据基底和模块M3共同生成模块M4对应的奇异值（该元素的精度此时并未实现
%   对齐），至此出去模块M4以外的三个模块均已实现阶数的对齐；
% 4.此时可以采用一个阈值作为循环涮选出符合要求的阶数对齐矩阵，由于只有一个元
%   素对应M4右上角的位置，因此很容易根据这个元素找到所有元素阶数均对齐的矩阵
% =================================================================================

function [RM, M_hat] = RMG_3PMP_Enhanced(M)
EPS_Threshold = 1;
ResultCheck = false;
NumCheck = 0;
M_hat_cache = [];
M_new = sort(abs(M(:)));
% =================================================================================
% 对输入矩阵做一个预处理，把其中的零元素先变成一个相对很小的非零数值
% =================================================================================
N = size(M,1);
M_Pre = M;
M_SecondMin = M_new(find(M_new>min(M_new),1));
if any(M_Pre(:) == 0)
    [row, col] = find(M_Pre == 0);
    for i = 1:numel(row)
        M_Pre(row(i), col(i)) = floor(log2(M_SecondMin)-52);
    end
end

while ~ResultCheck
    % 目标矩阵的阶数分布Ali_M和系数分布Rand_M
    Ali_M = floor(log10(abs(M_Pre)))+randi([-1,1],N);%阶数矩阵
    Rand_M = sign(M_Pre).*(rand(N,'double')+randi([1,9],N));
%     Rand_M = (rand(N,'double')+randi([1,9],N));
    
    % =================================================================================
    % 生成正负区间均有的随机系数分布矩阵
%     Rand_M = 2*rand(N,'double')-1 ;%step1:生成正态[-1,1]区间分布；
%     Rand_M = sign(Rand_M) + Rand_M;%step2：符号函数对应加成，首位数为1；
    % =================================================================================
    
    % 目标随机混淆矩阵Core_M
    Core_M = Rand_M.*10.^Ali_M;
    % 模块1 <==>  Block_M1
    Block_M1 = Core_M(1,1:N-1);
    % 模块2 <==>  Block_M2\
    Block_M2 = Core_M(2:end,1:N-1);
    % 模块3 <==>  Block_M3
    Block_M3 = Core_M(2:end,end);
    % 基向量首位元素分布 <==>  Vector_X
    Vector_X = Block_M2'\Block_M1';
    % 模块4 <==>  Block_M4
    Block_M4 = Block_M3'*Vector_X;
    % 更新最终的M_hat
    M_hat = [Block_M1 Block_M4; Block_M2 Block_M3];
    %生成随机RM
    RM = M_hat-M;
    % =================================================================================
    % 优化筛选模块（可选）
    % =================================================================================
    NumCheck = NumCheck + 1;
    if (abs(floor(log10(abs(M_hat(1,end))))-floor(log10(abs(M(1,end))))) > EPS_Threshold) || (M_hat(1,end)*M(1,end) < 0)
        fprintf('Warning! The element in M is %d but in M_hat is %d\n',M(1,end),M_hat(1,end));
        ResultCheck = false;
    else
        fprintf('Congratulation! The element in M is %d and in M_hat is %d\n',M(1,end),M_hat(1,end));
        ResultCheck = true;
    end
    M_hat_cache = cat(3,M_hat_cache,M_hat);
    % =================================================================================
    % 通过在有限次数内筛选出，模块4和标准阶数差距最小的情况作为最终输出的M_hat，
    % 此处决定了筛选的上届次数，这里是10000次
    % =================================================================================
    if NumCheck == 10^4
        Diff = abs(floor(log10(abs(M_hat_cache(1,end,:))))-floor(log10(abs(M(1,end)))));
        [OrderDiffValue,Index] = min(Diff);
        M_hat = M_hat_cache(:,:,Index);
        fprintf('Timeout! The element in M is %d and in M_hat is %d\n',M(1,end),M_hat(1,end));
        fprintf('The Minest Order Difference is %d\n',OrderDiffValue);
        break
    end
end
end


