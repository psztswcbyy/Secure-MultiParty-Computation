% =================================================================================
% 【函数功能说明】
%  RandMatrixGen3p函数对RandMatrixGen3p函数在3PMP协议任意秩但阶数对齐环境下
%  的测试，涉及到的输入参数包含待转换矩阵的维度N，矩阵元素数值精度的阶数最大
%  值maxEp和最小值minEp、首位数字的最大值FirstNumMin和最小值FirstNumMax；输
%  入数据元素的最大有效数字位数S-bit；输出参数包含M即转化前的原始矩阵，RM用
%  于混淆的随机矩阵，M_hat为转化后的随机矩阵；该函数仅用于3PMP内部测试精度使
%  用，不用作外部调用函数。
% =================================================================================

function [M, RM, M_hat] = RandMatrixGen3p(N, minEp, maxEp, FirstNumMin, FirstNumMax)
% =================================================================================
DynamicRange = 0;
% 该模块负责：尾数矩阵Ori_M生成，阶数矩阵Exp_M生成，eg：3.14E-12，其中3.14
% 代表尾数项，E-12代表阶数项
% =================================================================================
% % S作为有效位数实验的参数，生成最后一位有效位x_end的非零尾数值，介于[1,9]
% % 之间；生成S-1位有效位的x_front的数值，介于[1,10^(S-1)]之间；将前(s-1)位
% % 和第s位组合后移位加个位数元素数值，构成初始随机矩阵。
% x_end = randi(9,N);
% x_front = randi(10^(S-1),N);
% X = (x_front * 10 + x_end) * 10^(-S);
% Ori_M = X + randi([FirstNumMin,FirstNumMax],N,N);
% Exp_M = 10.^(randi([minEp, maxEp],N,N));
% M = Ori_M.*Exp_M;
% =================================================================================
% 原始矩阵M生成模块
% =================================================================================
% 原始矩阵为纯正数分布的情况，直接生成15位有效数字作为小数尾数的思路   
Ori_M = rand(N,N,'double') + randi([FirstNumMin,FirstNumMax],N,N);%step1:生成正态[0,1]区间分布；
Exp_M = 10.^(randi([minEp, maxEp],N,N));%生成原始阶数矩阵
M = Ori_M.*Exp_M;

% % 原始矩阵存在正负数分布的情况
% Ori_M = 2*rand(N,'double')-1 ;%step1:生成正态[-1,1]区间分布；
% Ori_M = sign(Ori_M).*randi([FirstNumMin, FirstNumMax], N, N) + Ori_M;%step2：符号函数对应加成，首位数为1；
% Exp_M = 10.^(randi([minEp, maxEp],N,N));%生成原始阶数矩阵
% M = Ori_M.*(Exp_M);

% 随机矩阵RM和M_hat生成模块，以2位对数，并且符号对齐，存在正负随机数分布
% =================================================================================
% 阶数对齐，其中Ali_M为M_hat的阶数对齐后的矩阵
Ali_M = floor(log2(M)) + randi([-DynamicRange,DynamicRange],N,N);
% Ori_M_hat为M_hat的原始尾数矩阵
Ori_M_hat = rand(N,N,'double') + randi([FirstNumMin,FirstNumMax],N,N);
% 生成混淆后的随机矩阵M_hat,由尾数矩阵和阶数矩阵构成
% M_hat = Ori_M_hat.*(2.^(Ali_M));
M_hat = sign(M).*Ori_M_hat.*(2.^(Ali_M));

% 生成随机矩阵RM
RM = M_hat - M;


end