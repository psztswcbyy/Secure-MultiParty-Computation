% =================================================================================
% RandMatrixGen2�������ڽ�ԭʼ���ݾ���ת��Ϊ�ض����漴���������漰��������
% ����������ת�������ά��N������Ԫ����ֵ���ȵĽ������ֵmaxEp����СֵminEp��
% ��λ���ֵ����ֵFirstNumMin����СֵFirstNumMax����������Ԫ�ص������Ч����
% λ��S-bit�������������M��ת��ǰ��ԭʼ����RM���ڻ������������M_hatת
% ������������
% ������˵�����ú���������2PMP�ڲ����Ծ���ʹ�ã��������ⲿ���ú���
% =================================================================================

function [M, RM, M_hat] = RandMatrixGen2(N, minEp, maxEp, FirstNumMin, FirstNumMax)
% =================================================================================
% 
% ��ģ�鸺��β������Ori_M���ɣ���������Exp_M���ɣ�eg��3.14E-12������3.14
% ����β���E-12����������
% =================================================================================
% % S��Ϊ��Чλ��ʵ��Ĳ������������һλ��Чλx_end�ķ���β��ֵ������[1,9]
% % ֮�䣻����S-1λ��Чλ��x_front����ֵ������[1,10^(S-1)]֮�䣻��ǰ(s-1)λ
% % �͵�sλ��Ϻ���λ�Ӹ�λ��Ԫ����ֵ�����ɳ�ʼ�������
% x_end = randi(9,N);
% x_front = randi(10^(S-1),N);
% X = (x_front * 10 + x_end) * 10^(-S);
% Ori_M = X + randi([FirstNumMin,FirstNumMax],N,N);
% Exp_M = 10.^(randi([minEp, maxEp],N,N));
% M = Ori_M.*Exp_M;
% =================================================================================
% ԭʼ����M����ģ��
% =================================================================================
% ֱ������15λ��Ч������ΪС��β����˼·   
Ori_M = rand(N,N,'double') + randi([FirstNumMin,FirstNumMax],N,N);
Exp_M = 10.^(randi([minEp, maxEp],N,N));
M = Ori_M.*Exp_M;
% =================================================================================
% vpa��ʾ��Ч���ֵ�λ��
% vpa(M,16)
% =================================================================================
% �������RM��M_hat����ģ��
% =================================================================================
% �������룬����Ali_MΪM_hat�Ľ��������ľ���
Ali_M = floor(log2(M)) + randi([-1,1],N,N);
% Ori_M_hatΪM_hat��ԭʼβ������
Ori_M_hat = rand(N,N,'double') + randi([FirstNumMin,FirstNumMax],N,N);
% ���ɻ�������������M_hat,��β������ͽ������󹹳�
M_hat = Ori_M_hat.*(2.^(Ali_M));
% �����������RM
RM = M_hat - M;
% =================================================================================
% % �޽������룬ֱ�������������ʱ���Ա�ʵ�鷽����
% M_hat = rand(N,N,'double')+randi([FirstNumMin,FirstNumMax],N,N);
% RM = M_hat - M;

end