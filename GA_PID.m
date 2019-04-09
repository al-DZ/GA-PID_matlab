%%清空环境
clc;
clear;
close all;
%% 参数设置
ChromosomeSize = 3;                                 %染色体长度，参数的个数，PID有三个参数
PopulationSize = 100;                               %种群规模
MaxIter = 100;                                      %最大迭代次数
MinFitness=0.01;                                    %最小适应值
% CrossRate=0.6;                                      %交叉概率
MutateRate=0.2;                                     %变异概率
ObjFun=@PSO_PID;                                    %适应值函数
NoChangeNo=5;                                              
UpLimit=30;                                              %上限
global Kp;
global Ki;
global Kd;
%% 初始化种群init.m
Population=rand(PopulationSize,ChromosomeSize);    %种群,预分配内存
for i=1:PopulationSize
    disp(['正在初始化种群',num2str(i)]);
    for j=1:ChromosomeSize
        Population(i,j)=rand*30;    %产生0-30之间的随机数用来初始化种群，30表示每个值不超过30
    end
end
clear i;
clear j;

%% 开始循环
%Iter =1;
PopulationFitness=zeros(PopulationSize,1);       %种群适应度值，预分配内存
BestFitness=zeros(MaxIter,1);                                   %初始化每一代的最佳适应度
AveFitness=zeros(MaxIter,1);
Elite=zeros(MaxIter,ChromosomeSize);                  %用于记录每一代的最优解

for Iter=1:MaxIter
    disp(['迭代次数：',num2str(Iter)]);         %显示迭代进度
    %% 适应值计算Fitness
    for i=1:PopulationSize
        Kp=Population(i,1);
        Ki=Population(i,2);
        Kd=Population(i,3);
        PopulationFitness(i,1) = fitness(Kp,Ki,Kd);
    end

    %% 适应值大小排序，并保存最佳个体和最佳适应度
    FitnessSum=sum(PopulationFitness);                 %种群累加适应度
    AveFitness(Iter,1)=FitnessSum/PopulationSize;           %种群平均适应度
    [PopulationFitness,Index]=sort(PopulationFitness);         %适应值从小到大排序
    BestFitness(Iter,1) = PopulationFitness(MaxIter,1);                %最佳适应度
    Elite(Iter,:) = Population(Index(MaxIter),:);                    %记录本代的精英
%     if(BestFitness(Iter,1)<MinFitness)                        %判断是否达到要求的适应值
%         break;
%     end
    
%     FitnessCumsum=cumsum(PopulationFitness);                   %累加适应度数组
    
    
   %复制适应值最大的不变的种群
    PopulationNew=zeros(PopulationSize,ChromosomeSize);
    for i=1:NoChangeNo
        PopulationNew(i,:)=Population(Index(MaxIter-i+1),:);
    end
    
    %轮盘赌法   选择selection要交叉的父母    
    for i=(NoChangeNo+1):2:PopulationSize
        %确定要交叉的父亲染色体序号
        idx1=0;
        idx2=0;
        m1=0;
        r1=rand*FitnessSum;
        for k=1:PopulationSize
            m1=m1+PopulationFitness(k);
            if r1<=m1
                idx1=Index(k);
                break;
            end
        end
         %确定要交叉的母亲染色体序号
        m2=0;
        r2=rand*FitnessSum;
        for k=1:PopulationSize
            m2=m2+PopulationFitness(k);
            if r2<=m2
                idx2=Index(k);
                break;
            end
        end
        acr_position = floor(ChromosomeSize*rand+1);                %要交叉的节点
        %交叉
        for j=1:acr_position
            temp = Population(idx1,j);
            Population(idx1,j) = Population(idx2,j);
            Population(idx2,j) = temp;
        end
        %将产生的两个子代添加到新的种群中
        PopulationNew(i,:)=Population(idx1,:);
        PopulationNew(i+1,:)=Population(idx2,:);
    end
    
    %变异  mutation
    for i=(NoChangeNo+1):PopulationSize
        for j=1:ChromosomeSize
            mut_rand = rand; %是否变异
            if mut_rand < MutateRate
                mut_pm = rand; %增加还是减少
                mut_num = rand*(1-AveFitness(Iter)/BestFitness(Iter))^2;
                if PopulationNew(i, j)>=UpLimit
                    PopulationNew(i, j)= PopulationNew(i, j)*(1-mut_num);
                elseif mut_pm<=0.5
                    PopulationNew(i, j)= PopulationNew(i, j)*(1-mut_num);
                else
                    PopulationNew(i, j)= PopulationNew(i, j)*(1+mut_num);
                end
                if PopulationNew(i, j)>=UpLimit
                    PopulationNew(i, j)=UpLimit;
                end
            end
        end
    end
	parfor i=1:PopulationSize
        for j=1:ChromosomeSize
            Population(i,j)=PopulationNew(i,j);
        end
	end
    clear i;
    clear j;
    clear first;
    clear last;
    clear idx;
    clear mid;
%     clear PopulationNew;
    
    %交叉操作  crossover
%     for i=1:PopulationSize
%         % rand<交叉概率，对两个个体的染色体串进行交叉操作
%         if(rand < CrossRate)
%             acr_position = floor(ChromosomeSize*rand+1);                %要交叉的节点
% %             if (cross_position == 0 || cross_position == 1)
% %                 continue;
% %             end
%             acr_chrom = floor((PopulationSize-1)*rand+1);                 %要交叉的染色体，floor取比它小的整数,acr_chrom取值在1-N    
%             for j=1:acr_position
%                 temp = Population(i,j);
%                 Population(i,j) = Population(acr_chrom,j);
%                 Population(acr_chrom,j) = temp;
%             end
%             
%         end
%     end
%     clear i;
% %     clear j;
%     clear temp;
%     clear acr_chrom;
    
    
    
    
    clear i;
    clear j;
    
    K_p(1,Iter)=Elite(Iter,1);
    K_i(1,Iter)=Elite(Iter,2);
    K_d(1,Iter)=Elite(Iter,3);
end

figure(1)
plot(BestFitness,'LineWidth',2);

figure(2)
plot(K_p)
hold on
plot(K_i,'k','LineWidth',3)
plot(K_d,'--r')






