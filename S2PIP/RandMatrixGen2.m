% =================================================================================
% RandMatrixGen2函数用于2PMP协议独立测试中生成随机矩阵Ra以及混淆矩阵A_hat。
% 涉及到的输入参数包含待转换矩阵的维度N，矩阵元素数值精度的阶数最大值maxEp
% 和最小值minEp、首位数字的最大值FirstNumMin和最小值FirstNumMax；输入数据元
% 素的最大有效数字位数S-bit；输出参数包含M即转化前的原始矩阵，RM用于混淆的随
% 机矩阵，M_hat转化后的随机矩阵；
% =================================================================================

function [M, RM, M_hat] = RandMatrixGen2(N, minEp, maxEp, FirstNumMin, FirstNumMax)
% =================================================================================
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
% 直接生成15位有效数字作为小数尾数的思路   
Ori_M = rand(N,N,'double') + randi([FirstNumMin,FirstNumMax],N,N);
Exp_M = 10.^(randi([minEp, maxEp],N,N));
Sign_M = sign(2*rand(N)-1);
M = Sign_M.*Ori_M.*Exp_M;

% =================================================================================
% 随机矩阵RM和M_hat生成模块
% =================================================================================
% 阶数对齐，其中Ali_M为M_hat的阶数对齐后的矩阵
M_new = sort(abs(M(:)));
M_SecondMin = M_new(find(M_new>min(M_new),1));
Ali_M = zeros(N,N);
for i=1:N
    for j=1:N
        if M(i,j) ~= 0
            Ali_M(i,j) = floor(log2(abs(M(i,j)))) + randi([-1,1]);
        else
            Ali_M(i,j) = floor(log2(M_SecondMin)-54);
        end
    end
end

% Ali_M = floor(log2(abs(M))) + randi([-1,1],N,N);

% Ori_M_hat为M_hat的原始尾数矩阵
Ori_M_hat = rand(N,N,'double') + randi([FirstNumMin,FirstNumMax],N,N);
% 生成混淆后的随机矩阵M_hat,由尾数矩阵和阶数矩阵构成
% M_hat = Ori_M_hat.*(2.^(Ali_M));
M_hat = sign(M).*Ori_M_hat.*(2.^(Ali_M));

% 生成随机矩阵RM
RM = M_hat - M;
% =================================================================================
% % 无阶数对齐，直接生成随机矩阵时（对比实验方案）
% M_hat = rand(N,N,'double')+randi([FirstNumMin,FirstNumMax],N,N);
% RM = M_hat - M;

end