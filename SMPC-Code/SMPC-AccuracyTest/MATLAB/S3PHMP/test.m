clc;
clear;
% N=5;
% r=N-1;
minEp = -8;
maxEp = 0;
% 
% Ori_M_hat = rand(r,r,'double') + randi([1,1],r,r);
% Exp_M_hat = 10.^(randi([minEp, maxEp],r,r));
% M_hat_base = Ori_M_hat.*Exp_M_hat;
% 
% Ori_Transformer_left = rand(N,r,'double') + randi([1,1],N,r);
% Exp_Transformer_left = 10.^(randi([0, 0],N,r));
% Transformer_left = Ori_Transformer_left.*Exp_Transformer_left;
% Ori_Transformer_right = rand(r,N,'double') + randi([1,1],r,N);
% Exp_Transformer_right = 10.^(randi([0, 0],r,N));
% Transformer_right = Ori_Transformer_right.*Exp_Transformer_right;
% M_hat = Transformer_left*M_hat_base*Transformer_right;
% 
% Rank_MHB = rank(M_hat_base)
% Rank_TL = rank(Transformer_left)
% Rank_TR = rank(Transformer_right)
% Rank_M_hat = rank(M_hat_base)

% P = randi([1,9],5,4);
% Q = randi([1,9],4,5);
% M = randi([1,9],4,4);
% M_hat = P*M*Q;
% rank(M_hat)
% rank(randi([1,1000],1)*M_hat+randi([1,1000],1))

n = 10^(maxEp+1);
m = 10^minEp;
n-m
delta = zeros(1,5);
index = randperm(5,4);
delta(index) = randi([-100,-1],1,4);
P = rand(5);
M = P*diag(delta)*inv(P);
a = min(min(M))
b = max(max(M))
if (max(delta) < (m*b-a*n)/(b-a))
    Standard = true;
else
    Standard = false;
end
e = delta(index(randi(4)))
k = (m-e)/a
M_hat = M*k+e*eye(5);
% =========================================================================================
% �ò��Է���ʧ��ԭ�������
% =========================================================================================
% ��������ͨ����ԭʼ�������ֵ���䣬����һ�����Է�����ӳ�䵽ָ���Ľ�����Χ�����ڣ�����ʵ��
% ������������ӳ���еı�����ϵ���Ա�֤�ڱ任�����У�������Ȳ��䣬���ǽؾ��Ӧ��ƽ�ƹ�ϵ
% ȴ�޷���֤�ȵĲ����ԣ�����Y=k*X+��*I�У��Ҿ��ǽؾֻ࣬�е��ҵ��ھ���X�ĸ�����ֵʱ������
% ��֤������Ȳ��������ȣ��������ǳ��������㣬��Ϊ�󲿷�������ǦҲ���������ֵ�ģ����֮
% �������õ��Ķ���һ�����Ⱦ�����˸���һ������ķ����Ⱦ�����ͼ�������ӳ�䵽һ����������
% ��һ���ǳ����ѵ����飬�÷�����ʱ���ú����п����о�






