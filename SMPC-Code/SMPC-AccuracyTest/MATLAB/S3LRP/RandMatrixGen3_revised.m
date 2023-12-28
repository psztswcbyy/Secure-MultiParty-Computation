% =========================================================================================
% 【函数功能说明】
%  RandMatrixGen3函数根据原始矩阵生成秩为r的随机混淆矩阵M_hat,涉及到的输入参
%  数包含待转换矩阵的维度N，待转换矩阵的秩r，矩阵元素数值精度的阶数最大值和 
%  最小值、首位数字的最大最小值；输出参数包含M即转化前的原始矩阵，RM用于混淆
%  的随机矩阵，M_hat转化后的随机矩阵；
%  该函数由于用作测试精度用，因此采用先生成混淆后矩阵然后针对该矩阵的精度分
%  布去设计原始矩阵，以控制原始矩阵和混淆后矩阵的精度配准问题（仅在精度试验
%  中可用）。
% 【使用注意事项】该函数仅用于3PMP内部测试精度使用，不用作外部调用函数，且该
%  函数是根据混淆矩阵M_hat反过来生成M得出RM，因此只在3PMP内部精度测试用，而
%  在实际调用中M是事先确定好的，因此需要另外采用RandMatrix3s；
% =========================================================================================


function [M, RM, M_hat] = RandMatrixGen3_revised(N,minEp,maxEp,FirstNumMin,FirstNumMax)
% =========================================================================================
% 生成最后一位有效位x_end的非零尾数值，介于[1,9]之间；生成前面部分S-1位有效
% 位的x_front的数值，介于[0,10^(S-1)-1]之间；将前(s-1)位和第s位组合后移位加
% 个位数元素数值，构成混淆后的随机矩阵。该部分用于测试有效位数对精度的影响。
% x_end = randi([1,9],r);
% x_front = randi([0,10^(S-1)-1],r);
% X = (x_front * 10 + x_end) * 10^(-S);
% Ori_M_hat = X + randi([FirstNumMin,FirstNumMax],r);
% Exp_M_hat = 10.^(randi([minEp, maxEp],r));
% M_hat_base = Ori_M_hat.*Exp_M_hat;
% =========================================================================================
% 1.生成用于混淆秩为r的满秩方阵M_hat_base
% =========================================================================================
% 通过构造M_hat=Transformer_left*M_hat_base*Transformer_right生成一个秩为r的
% 维度为N的非满秩方阵。
% =========================================================================================
r = N-1;
Ori_M_hat = rand(r,r,'double') + randi([FirstNumMin,FirstNumMax],r,r);
Exp_M_hat = 10.^(randi([minEp, maxEp],r,r));
M_hat_base = Ori_M_hat.*Exp_M_hat;

% =========================================================================================
% 2.生成用于混淆列满秩左侧变换阵M_hat_left，以及行满秩右侧变换阵M_hat_right
% =========================================================================================
Ori_Transformer_left = rand(N,r,'double') + randi([FirstNumMin,FirstNumMax],N,r);
Exp_Transformer_left = 10.^(randi([1, 1],N,r));
Transformer_left = Ori_Transformer_left.*Exp_Transformer_left;
Ori_Transformer_right = rand(r,N,'double') + randi([FirstNumMin,FirstNumMax],r,N);
Exp_Transformer_right = 10.^(randi([1, 1],r,N));
Transformer_right = Ori_Transformer_right.*Exp_Transformer_right;
M_hat = Transformer_left*M_hat_base*Transformer_right;

% Rank_MHB = rank(M_hat_base)
% Rank_TL = rank(Transformer_left)
% Rank_TR = rank(Transformer_right)
% Rank_M_hat = rank(M_hat_base)

% =========================================================================================
% 3.阶数对齐，其中Ali_M为原始矩阵M与M_hat对应的阶数对齐矩阵
% =========================================================================================
% 【需注意的是，这里先生成随机混淆后的结果矩阵M_hat然后根据M_hat的阶数分布去生成原始矩阵，
%  属于通过反向推理的方式保证初始矩阵M和混淆后矩阵M_hat阶数对齐的实验条件，严格意义上说，
%  应该是根据原始输入矩阵来生成混淆矩阵的，这里需要注意一下！】
Ali_M = floor(log2(abs(M_hat))) + randi([-1,1],N,N);
% Ori_M_hat为M_hat的原始尾数矩阵
Ori_M_hat = rand(N,N,'double') + randi([FirstNumMin,FirstNumMax],N,N);
% 生成混淆后的随机矩阵M_hat,由尾数矩阵和阶数矩阵构成
M = Ori_M_hat.*(2.^(Ali_M));

% =========================================================================================
% 4.生成随机矩阵RM
% =========================================================================================
RM = M_hat - M;

% =========================================================================================
% 无阶数对齐，直接生成随机矩阵时（对比实验方案）
% M_hat = rand(N,N,'double')+randi([FirstNumMin,FirstNumMax],N,N);
% RM = M_hat - M;

end