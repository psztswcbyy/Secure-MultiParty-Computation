% =================================================================================
% RMG_2PMP�������ڸ��϶෽�����2PMPЭ��ĵ��ã����������������M
% ��ת��ǰ��ԭʼ����RM���ڻ������������M_hatת������������
% =================================================================================
clc;
clear;
format longE
N = 10;
minEp = -10;
maxEp = -minEp;
FirstNumMin = 1;
FirstNumMax = 1;

% A = RandInputMatrix(N,FirstNumMin,FirstNumMax,minEp,maxEp);
% ����2PMPЭ���޷����ȵ�Ҫ����˿���ֱ����RandMatrixGen2���ɳ�ʼ����
A = RandMatrixGen2(N,minEp, maxEp, FirstNumMin, FirstNumMax);
[Ra, A_hat] = RFHMGen(A)

function [RM, M_hat] = RFHMGen(M)
% =================================================================================
% �������RM��M_hat����ģ��
% =================================================================================
% �������룬����Ali_MΪM_hat�Ľ��������ľ���
Check=false;
Num = 0;
while Check == false
    Num = Num + 1;
    [m, n] = size(M);
    % ����ԭʼ�����п��ܴ���0Ԫ�أ�������Ҫ���н����������⴦������������з�
    % ��Ԫ����СֵΪ�ο���׼���������������е���Ԫ�ض�Ӧ����Ϊ��16�׵�Ԫ�ء�
    M_new = sort(abs(M(:)));
    M_SecondMin = M_new(find(M_new>min(M_new),1));
    Ali_M = zeros(m,n);
    for i=1:m
        for j=1:n
            
            if M(i,j) ~= 0
                Ali_M(i,j) = floor(log2(abs(M(i,j)))) + randi([-1,1]);
            else
                Ali_M(i,j) = floor(log2(M_SecondMin)-54);
            end
        end
    end
    % Ori_M_hatΪM_hat��ԭʼβ������
    Ori_M_hat = rand(m,n,'double') + randi([1,1],m,n);
    % ���ɻ�������������M_hat,��β������ͽ������󹹳�
    M_hat = sign(M).* Ori_M_hat.* (2.^(Ali_M));
    if Num <= 10^4
        if cond(M_hat)>10^3
            Check = false;
        else
            fprintf('Find this matrix\n')
            RM = M_hat - M;
            break;
        end
    else
        fprintf('No matrix satisfy\n')
        RM = M_hat - M;
        break
    end
    % �����������RM
    % =================================================================================
    % % �޽������룬ֱ�������������ʱ���Ա�ʵ�鷽����
    % M_hat = rand(N,N,'double')+ones(N);
    % RM = M_hat - M;
end
end