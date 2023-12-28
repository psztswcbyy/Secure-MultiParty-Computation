% =================================================================================
% �����Ժ�������˵����
% RBMG_2PMP���ڸ������������Խ�ԭ����������һ��������У����ұ�֤��������
% �ͷ����ȡ��ú�������������ά������ľ��������
% ����������ʵ����Է���˵����
% �÷���ͨ����Nάԭʼ����M�����������γ��ĸ��Խǿ���ֲ���ԭʼ����Ĵ����
% Big_M������λ��Ϊ0����Ȼ�������������RM,�����Ÿ�����飬ÿ������鶼��
% �к�ԭʼ����M��ͬ�ľ��ȷֲ������ɸѡBig_M�г�ȥ�ĸ��Խǿ�λ�õ��м䲿�־�
% ����е�����������Ϊ���Թ�ϵ��ʹ��RM�����ȣ��ò���������һ������ʽ����
% ��ǿԼ������ͨ��RM+Big_M���ɻ�����յķ����Ⱦ��ȶ������
% =================================================================================
% ����ר�ò���---���ڶԺ���RBMG�Ĳ���
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
% ����չ����Big_M����ģ�顿
% =================================================================================
% ֱ������15λ��Ч������ΪС��β����˼·,���ﰴ��3x3�ķֿ������ã��ĸ��Խ�
% �о���Ԫ�أ�����λ��Ϊȫ0Ԫ��
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
% ���������RM��M_hat����ģ�顿
% =================================================================================
% �������룬Ali_MΪBig_M�������뾫�ȷֲ�����,OAM�ǽ����������OrderAlignMatrix���
OAM=zeros(DimM, DimN);
M_new = sort(abs(M(:)));
M_SecondMin = M_new(find(M_new>min(M_new),1));
%�Ӿ����Ľ����ֲ�
for i=1:DimM
    for j=1:DimN
        if M(i,j) ~= 0
            OAM(i,j) = floor(log2(abs(M(i,j)))) + randi([-1,1]);
        else
            OAM(i,j) = floor(log2(M_SecondMin)-52);
        end
    end
end
%��ֵ�����Ľ����ֲ�����Ӧ�Ÿ���

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
% ���ɻ�������������M_hat,��β������ͽ������󹹳�
RM = RM.*(2.^(Ali_M));
% =================================================================================
% �������Ⱦ����Ż�ģ�顿���޷�����Լ���Ŀ���ȥ����ģ�飩
% =================================================================================

RM_Test = RM;
RowGroup = DimM+1:1:DimM+Dx;
while true
    Index_row = randperm(numel(RowGroup),2);
    RM_Test(RowGroup(Index_row(1)),:)= RM(RowGroup(Index_row(2)),:)*rand(1);
%     floor(log10(abs(det(RM_Test))));
%     rank(RM_Test);
    while (floor(log10(abs(det(RM_Test)))) > -16)%E-16��������ʽ����0�������̶�
        %�����������һ���Ѿ��������б�����ϵ�󣬼�������ѡ���б��������ܻ��һ�����;������
        Index_row = randperm(numel(RowGroup),2);
        RM_Test(RowGroup(Index_row(1)),:)= RM_Test(RowGroup(Index_row(2)),:)*rand(1);
%         det(RM_Test)%��ӡ������
%         floor(log10(abs(det(RM_Test))));
%         rank(RM_Test);
    end
    break
end

% �����������RM
M_hat = Big_M + RM_Test;
% det(RM_Test)
% det(M_hat)
%��ӡ�������������ȣ���֤�Ƿ������
% fprintf('The rank of A_hat is:%d \n',rank(M_hat))

% =================================================================================
% ���޽������룬ֱ�������������ʱ���Ա�ʵ�鷽������
% =================================================================================
% M_hat = rand(N,N,'double')+randi([FirstNumMin,FirstNumMax],N,N);
% RM = M_hat - M;
end