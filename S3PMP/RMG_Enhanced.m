% =================================================================================
% ����������˵����RMG_Enhanced�������ڸ��Ϻ�����3PMP�����ĵ��ã��˷���
% �����Ͻ�Ԫ�صľ�����ʧ������ֵEPS_Threshold��Ϊɸѡ��׼������������Ԫ�ؽ�
% ���ϵ���ȫ���룩�������������RM���ڻ������������M_hatת������������
% ����������ʵ����Է���˵����
% ͨ�������ĸ���ģ�飬Ѱ�����ŵĻ������飬�����������ͨ���ض���������Ͽ��ԡ�
% ����ָ�����ȷֲ��ķ����Ⱦ���ÿ��ģ��֮��������Ӧ�Ĺ�ϵ��
% 1.������ԭʼ����M�Ľ����ֲ����ɶ�Ӧ����Ľ�������ģ��Core_M;
% 2.Ȼ�����ģ��M1��ģ��M2֮����ڵ����Թ�ϵ�������������飨b1,b2...bn������
%   ֵ�⣻
% 3.�����ݻ��׺�ģ��M3��ͬ����ģ��M4��Ӧ������ֵ����Ԫ�صľ��ȴ�ʱ��δʵ��
%   ���룩�����˳�ȥģ��M4���������ģ�����ʵ�ֽ����Ķ��룻
% 4.��ʱ���Բ���һ����ֵ��Ϊѭ����ѡ������Ҫ��Ľ��������������ֻ��һ��Ԫ
%   �ض�ӦM4���Ͻǵ�λ�ã���˺����׸������Ԫ���ҵ�����Ԫ�ؽ���������ľ���
% =================================================================================
% ���Դ���
% =================================================================================
% clc;
% clear;
% N=5;
% minEp = -10;
% maxEp = 10;
% FirstNumMin = 1;
% FirstNumMax = 1;
% Ori_M = rand(N,N,'double') + randi([FirstNumMin,FirstNumMax],N,N);%step1:������̬[1,2]����ֲ���
% Exp_M = 10.^(randi([minEp, maxEp],N,N));%����ԭʼ��������
% M = Ori_M.*Exp_M;
% %��������ʾ
% [RM, M_hat, Diff] = RMGEnhanced(M)
% M
% if rank(M_hat)  < N
%     fprintf('The rank of M_hat is:%d and the matrix M_hat is rank-deficient matrix\n',rank(M_hat))
% end
% 
% if max(RM+M-M_hat)<N*eps(class(M))*norm((M_hat), 'inf')
%     fprintf('The Test Programme is 0 bug, and the error is: \n')
%     disp(M_hat-M-RM)
% end

% =================================================================================
function [RM, M_hat,OrderDiffValue] = RMG_Enhanced(M)
EPS_Threshold = 0;
DynamicRange = 0;
ResultCheck = false;
NumCheck = 0;
M_hat_cache = [];
% =================================================================================
% �����������һ��Ԥ���������е���Ԫ���ȱ��һ����Ժ�С�ķ�����ֵ
% =================================================================================
N = size(M,1);
M_Pre = M;
M_new = sort(abs(M(:)));
M_SecondMin = M_new(find(M_new>min(M_new),1));
if any(M_Pre(:) == 0)
    [row, col] = find(M_Pre == 0);
    for i = 1:numel(row)
        M_Pre(row(i), col(i)) = floor(log2(M_SecondMin)-52);
    end
end

while ~ResultCheck
    NumCheck = NumCheck + 1;
    % Ŀ�����Ľ����ֲ�Ali_M��ϵ���ֲ�Rand_M
    Ali_M = floor(log2(abs(M_Pre)))+randi([-DynamicRange,DynamicRange],N);%��������
    Rand_M = sign(M_Pre).*(rand(N,'double')+randi([1,9],N));
%     Rand_M = (rand(N,'double')+randi([1,9],N));
    
    % =================================================================================
    % ��������������е����ϵ���ֲ�����
%     Rand_M = 2*rand(N,'double')-1 ;%step1:������̬[-1,1]����ֲ���
%     Rand_M = sign(Rand_M) + Rand_M;%step2�����ź�����Ӧ�ӳɣ���λ��Ϊ1��
    % =================================================================================
    
    % Ŀ�������������Core_M
    Core_M = Rand_M.*2.^Ali_M;
    % ģ��1 <==>  Block_M1
    Block_M1 = Core_M(1,1:N-1);
    % ģ��2 <==>  Block_M2\
    Block_M2 = Core_M(2:end,1:N-1);
    % ģ��3 <==>  Block_M3
    Block_M3 = Core_M(2:end,end);
    % ��������λԪ�طֲ� <==>  Vector_X
%     Vector_X = Block_M2'\Block_M1';
    % �������Ż�����
    [P,R,C] = equilibrate(Block_M2');
    tol = 1e-12;
    maxit = size(Block_M2, 1);
    newBlock_M2 = R * P * Block_M2' * C;
    newBlock_M1 = R * P * Block_M1';
    [y,~,~,~,~] = gmres(newBlock_M2,newBlock_M1,[],tol,maxit);
    Vector_X = C * y;
    
    % ģ��4 <==>  Block_M4
    Block_M4 = Block_M3'*Vector_X;
    % �������յ�M_hat
    M_hat = [Block_M1 Block_M4; Block_M2 Block_M3];
    % =================================================================================
    % �Ż�ɸѡģ�飨��ѡ��
    % =================================================================================
%     NumCheck = NumCheck + 1;
    if (abs(floor(log10(abs(M_hat(1,end))))-floor(log10(abs(M(1,end))))) > EPS_Threshold) || (M_hat(1,end)*M(1,end) < 0)
%         fprintf('Warning! The element in M is %d but in M_hat is %d\n',M(1,end),M_hat(1,end));
        ResultCheck = false;
    else
        fprintf('Congratulation! The element in M is %d and in M_hat is %d\n',M(1,end),M_hat(1,end));
        ResultCheck = true;
        %�������RM
        RM = M_hat-M;
        OrderDiffValue = abs(floor(log10(abs(M_hat(1,end))))-floor(log10(abs(M(1,end)))));
    end
   
    M_hat_cache = cat(3,M_hat_cache,M_hat);
    % =================================================================================
    % ͨ�������޴�����ɸѡ����ģ��4�ͱ�׼���������С�������Ϊ���������M_hat��
    % �˴�������ɸѡ���Ͻ������������10000��
    % =================================================================================
    if NumCheck == 10^2
        
        CacheM = M_hat_cache(1,end,:);
        
        if size(find(CacheM*M(1,end) > 0),1)~=0
            Index_z = find(CacheM*M(1,end) > 0);
        else
            fprintf('Please note that the matrix lacks complete exponential congruence!\n');
            Diff = abs(M_hat_cache(1,end,:)-M(1,end));
            [OrderDiffValue,Index] = min(Diff);
            M_hat = M_hat_cache(:,:,Index);
            RM = M_hat-M;
            fprintf('Timeout! The element in M is %d and in M_hat is %d\n',M(1,end),M_hat(1,end));
            fprintf('The Minest Numerical Difference is %d\n',OrderDiffValue);
            break
        end
        
        Diff = abs(floor(log10(abs(M_hat_cache(1,end,Index_z(:)))))-floor(log10(abs(M(1,end)))));
        [OrderDiffValue,Index_Diff] = min(Diff);
        
        M_hat = M_hat_cache(:,:,Index_z(Index_Diff));
        RM = M_hat-M;
        
        fprintf('Timeout! The element in M is %d and in M_hat is %d\n',M(1,end),M_hat(1,end));
        fprintf('The Minest Order Difference is %d\n',OrderDiffValue);
        break
    end
end
end


