% =================================================================================
% 【测试函数功能说明】
% RBMG_2PMP用于根据输入矩阵可以将原矩阵扩增到一个大矩阵中，并且保证阶数对齐
% 和非满秩。该函数适用于任意维度输入的矩阵情况。
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
% 
% DimA = randi([2,10],1,2);
% M = DimA(1);
% N = DimA(2);
% minEp=-5;
% maxEp=0;
% FirstNumMin=1;
% FirstNumMax=1;
% Ori_A = rand(M,N,'double') + randi([FirstNumMin,FirstNumMax],M,N);
% Exp_A = 10.^(randi([minEp, maxEp],M,N));
% A = Ori_A.*Exp_A
% 
% [Big_A, Ra, A_hat] = RBMG_Anything(A)
% =================================================================================
function [Big_M, RM, M_hat] = RBMG_Anything(M)
% =================================================================================
% 【扩展矩阵Big_M生成模块】
% =================================================================================
% 直接生成15位有效数字作为小数尾数的思路,这里按照3x3的分块来设置，四个对角
% 有矩阵元素，其余位置为全0元素
DimM = size(M,1);
DimN = size(M,2);
if DimM > DimN
    Dx = ceil((DimM + DimN)/2);
    Dy = Dx + 2*(DimM - DimN);
else
    Dy = ceil((DimM + DimN)/2);
    Dx = Dy - 2*(DimM - DimN);
end
NewM = 2*DimM + Dx;
NewN = 2*DimN + Dy;
Big_M = zeros(NewM, NewN);
Big_M(1:DimM, 1:DimN) = M;
Big_M(1:DimM, DimN+Dy+1:end) = M;
Big_M(DimM+Dx+1:end, 1:DimN) = M;
Big_M(DimM+Dx+1:end, DimN+Dy+1:end) = M;
% =================================================================================
% 【随机矩阵RM和M_hat生成模块】
% =================================================================================
% 阶数对齐，Ali_M为Big_M阶数对齐精度分布矩阵,OAM是阶数对齐矩阵OrderAlignMatrix简称
OAM=zeros(DimM, DimN);
M_new = sort(abs(M(:)));
M_SecondMin = M_new(find(M_new>min(M_new),1));
%子矩阵块的阶数分布
for i=1:DimM
    for j=1:DimN
        if M(i,j) ~= 0
            OAM(i,j) = floor(log2(abs(M(i,j)))) + randi([-1,1]);
        else
            OAM(i,j) = floor(log2(M_SecondMin)-52);
        end
    end
end
%赋值大矩阵的阶数分布，对应九个块

Ali_M = randi([min(min(OAM)),max(max(OAM))], NewM, NewN);
Sign_M = sign(2*randn(NewM, NewN, 'double') - 1);


Ali_M(1:DimM,1:DimN)=OAM;
Ali_M(1:DimM, DimN+Dy+1:end) = OAM;
Ali_M(DimM+Dx+1:end, 1:DimN) = OAM;
Ali_M(DimM+Dx+1:end, DimN+Dy+1:end) = OAM;
Sign_M(1:DimM,1:DimN)=sign(M);
Sign_M(1:DimM, DimN+Dy+1:end) = sign(M);
Sign_M(DimM+Dx+1:end, 1:DimN) = sign(M);
Sign_M(DimM+Dx+1:end, DimN+Dy+1:end) = sign(M);

RM = Sign_M.*(rand(NewM,NewN,'double') + randi([1,1],NewM,NewN));
% 生成混淆后的随机矩阵M_hat,由尾数矩阵和阶数矩阵构成
RM = RM.*(2.^(Ali_M));
% =================================================================================
% 【非满秩矩阵优化模块】（无非满秩约束的可以去掉该模块）
% =================================================================================

RM_Test = RM;
RowGroup = DimM+1:1:DimM+Dx;
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