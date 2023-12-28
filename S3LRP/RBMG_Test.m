% =================================================================================
% 【测试函数功能说明】
% RBMG_Test函数用于将维度N的原始数据矩阵转化为维度扩大的随机矩阵Big_M，涉及
% 到的输入参数包含原始矩阵的维度N，矩阵元素数值精度的阶数最大值maxEp和最小
% 值minEp、首位数字的最大值FirstNumMin和最小值FirstNumMax；输入数据元素的最
% 大有效数字位数S-bit；输出参数包含M即转化前的原始矩阵，RM用于混淆的随机矩阵，
% M_hat转化后的随机矩阵；
% 【阶数对齐实验测试方案说明】
% 该方法通过对N维原始矩阵M进行扩增，形成四个对角块均分布有原始矩阵的大矩阵
% Big_M，其余位置为0矩阵，然后生成随机矩阵RM,包含九个矩阵块，每个矩阵块都具
% 有和原始矩阵M相同的精度分布，随机筛选Big_M中除去四个对角块位置的中间部分矩
% 阵块中的任意两行作为线性关系，使得RM非满秩（该部分增加了一个行列式误差极限
% 得强约束）。通过RM+Big_M即可获得最终的非满秩精度对齐矩阵。
% =================================================================================
% 测试专用参数---用于对函数RBMG的测试
% =================================================================================
% clc;clear;
% N=20;
% minEp=-5;
% maxEp=0;
% FirstNumMin=1;
% FirstNumMax=1;
% 创建一个mxn且秩为r的非满秩随机矩阵
% A = (randi([1,10],m,r)*rand(r,r)*randi([1,10],r,n)).*10.^(randi([minEp, maxEp],m,n));
% 
% [A, Ra, A_hat] = RBMG_Test(N,minEp,maxEp,FirstNumMin,FirstNumMax);
% =================================================================================
function [Big_M, RM, M_hat] = RBMG_Test(N, minEp, maxEp, FirstNumMin, FirstNumMax)
% =================================================================================
% 该模块负责：尾数矩阵Ori_M生成，阶数矩阵Exp_M生成，eg：3.14E-12，其中3.14
% 代表尾数项，E-12代表阶数项
% =================================================================================
% 【原始矩阵M生成模块】
% =================================================================================
% 直接生成15位有效数字作为小数尾数的思路
% =================================================================================
% % 生成正负区间均有的随机系数分布矩阵
% %step1:生成正态[-1,1]区间分布；
% Ori_M = 2*rand(N,'double')-1 ;
% %step2：符号函数对应加成，首位数为1；
% Ori_M = sign(Ori_M).*randi([FirstNumMin,FirstNumMax],N,N) + Ori_M;
% =================================================================================
Ori_M = rand(N,N,'double') + randi([FirstNumMin,FirstNumMax],N,N);
Exp_M = 10.^(randi([minEp, maxEp],N,N));
M = Ori_M.*Exp_M;

% =================================================================================
% 【扩展矩阵Big_M生成模块】
% =================================================================================
% 直接生成15位有效数字作为小数尾数的思路,这里按照3x3的分块来设置，四个对角
% 有矩阵元素，其余位置为全0元素
Big_M = zeros(3 * N);
Big_M(1:N, 1:N) = M;
Big_M(1:N, 2*N+1:end) = M;
Big_M(2*N+1:end, 1:N) = M;
Big_M(2*N+1:end, 2*N+1:end) = M;
% =================================================================================
% 【随机矩阵RM和M_hat生成模块】
% =================================================================================
% 阶数对齐，Ali_M为Big_M阶数对齐精度分布矩阵,OAM是阶数对齐矩阵OrderAlignMatrix简称
OAM=zeros(size(M,1));
M_new = sort(abs(M(:)));
M_SecondMin = M_new(find(M_new>min(M_new),1));
%子矩阵块的阶数分布
for i=1:N
    for j=1:N
        if M(i,j) ~= 0
            OAM(i,j) = floor(log2(abs(M(i,j)))) + randi([-1,1]);
        else
            OAM(i,j) = floor(log2(M_SecondMin)-52);
        end
    end
end
%赋值大矩阵的阶数分布，对应九个块
Ali_M = zeros(3*N,3*N,'double');
Sign_M = zeros(3*N,3*N,'double');
for i=1:3
    for j=1:3
        row_start = (i-1)*N+1;
        col_start = (j-1)*N+1;
        Ali_M(row_start:row_start+N-1,col_start:col_start+N-1)=OAM;
        Sign_M(row_start:row_start+N-1,col_start:col_start+N-1)=sign(M);
    end
end

RM = Sign_M.*(rand(3*N,3*N,'double') + randi([FirstNumMin,FirstNumMax],3*N,3*N));
% 生成混淆后的随机矩阵M_hat,由尾数矩阵和阶数矩阵构成
RM = RM.*(2.^(Ali_M));
% =================================================================================
% 【非满秩矩阵优化模块】（无非满秩约束的可以去掉该模块）
% =================================================================================

RM_Test = RM;
RowGroup = N+1:1:2*N;
while true
    Index_row = randperm(numel(RowGroup),2);
    RM_Test(RowGroup(Index_row(1)),:)= RM(RowGroup(Index_row(2)),:)*rand(1);
%     floor(log10(abs(det(RM_Test))));
%     rank(RM_Test);
    while (floor(log10(abs(det(RM_Test)))) > -16)%E-16决定行列式趋于0的缩进程度
        %该语句是在上一步已经存在两行比例关系后，继续打乱选两行比例，可能会进一步降低矩阵的秩
        Index_row = randperm(numel(RowGroup),2);
        RM_Test(RowGroup(Index_row(1)),:)= RM_Test(RowGroup(Index_row(2)),:)*rand(1);
%         det(RM_Test)%打印测试用
%         floor(log10(abs(det(RM_Test))));
%         rank(RM_Test);
    end
    break
end

% 生成随机矩阵RM
M_hat = Big_M + RM_Test;
% det(RM_Test)
% det(M_hat)
%打印输出混淆矩阵的秩，验证是否非满秩
% fprintf('The rank of A_hat is:%d \n',rank(M_hat))

% =================================================================================
% 【无阶数对齐，直接生成随机矩阵时（对比实验方案）】
% =================================================================================
% M_hat = rand(N,N,'double')+randi([FirstNumMin,FirstNumMax],N,N);
% RM = M_hat - M;
end