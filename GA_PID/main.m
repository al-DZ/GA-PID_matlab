%%清空环境
clc;
clear;
close all;
%% 参数设置
ChromosomeSize = 1;                                 %染色体个数
ChromosomeLen=17;                                   %染色体长度   由最大值决定  
PopulationSize = 100;                               %种群规模
MaxIter = 50;                                      %最大迭代次数
% MinFitness=0.01;                                    %最小适应值
CrossRate=0.6;                                      %交叉概率
MutateRate=0.1;                                     %变异概率
% ObjFun=@PSO_PID;                                    %适应值函数
NoChangeNo=5;                                              
% UpLimit=30;                                              %上限
global Kp;
global Ki;
global Kd;
%% 初始化种群init.m
Population1=rand(PopulationSize,ChromosomeLen);    %种群,预分配内存   
Population2=rand(PopulationSize,ChromosomeLen);    %种群,预分配内存
Population3=rand(PopulationSize,ChromosomeLen);    %种群,预分配内存
for i=1:PopulationSize
    disp(['正在初始化种群',num2str(i)]);
    for j=1:ChromosomeLen
        Population1(i,j)=round(rand);   
    end 
end
for i=1:PopulationSize
    disp(['正在初始化种群',num2str(i)]);
    for j=1:ChromosomeLen
        Population2(i,j)=round(rand);   
    end 
end
for i=1:PopulationSize
    disp(['正在初始化种群',num2str(i)]);
    for j=1:ChromosomeLen
        Population3(i,j)=round(rand);   
    end 
end
clear i;
clear j;

%% 开始循环
PopulationFitness=zeros(PopulationSize,1);       %种群适应度值，预分配内存
BestFitness=zeros(MaxIter,1);                                   %初始化每一代的最佳适应度
AveFitness=zeros(MaxIter,1);                               %初始化每一代的平均适应度
K_p=zeros(MaxIter,1);                        %初始化 用于提高运算速度
K_i=zeros(MaxIter,1);                        %初始化 用于提高运算速度
K_d=zeros(MaxIter,1);                        %初始化 用于提高运算速度
Elite1=zeros(MaxIter,1);                  %用于记录每一代的最优解
Elite2=zeros(MaxIter,1);                  %用于记录每一代的最优解
Elite3=zeros(MaxIter,1);                  %用于记录每一代的最优解
for Iter=1:MaxIter
    disp(['迭代次数：',num2str(Iter)]);         %显示迭代进度
    %% 适应值计算Fitness
    for i=1:PopulationSize
%         PopulationFitness(i,1) = fitness(decode(Population1(i,:)));
        Kp=decode(Population1(i,:));
        Ki=decode(Population2(i,:));
        Kd=decode(Population3(i,:));
        PopulationFitness(i,1) = fitness(Kp,Ki,Kd);
    end
    %% 适应值大小排序，并保存最佳个体和最佳适应度
    FitnessSum=sum(PopulationFitness);                 %种群累加适应度
    AveFitness(Iter,1)=FitnessSum/PopulationSize;           %种群平均适应度
    [PopulationFitness,Index]=sort(PopulationFitness);         %适应值从小到大排序
    BestFitness(Iter,1) = PopulationFitness(PopulationSize,1);                %记录每一代的最佳适应度 
    Elite1(Iter,1) = decode(Population1(Index(PopulationSize),:));                    %记录本代的精英
    Elite2(Iter,1) = decode(Population2(Index(PopulationSize),:));                    %记录本代的精英
    Elite3(Iter,1) = decode(Population3(Index(PopulationSize),:));                    %记录本代的精英
    disp(['最佳适应度：',num2str(BestFitness(Iter,1))]);
    disp(['最佳个体：',num2str(Elite1(Iter,1)),' ',num2str(Elite2(Iter,1)),' ',num2str(Elite3(Iter,1))]);
     %% 复制适应值最大的不变的染色体
     PopulationNew1=zeros(PopulationSize,ChromosomeLen);                %初始化新的种群
     PopulationNew2=zeros(PopulationSize,ChromosomeLen);                %初始化新的种群
     PopulationNew3=zeros(PopulationSize,ChromosomeLen);                %初始化新的种群
    for i=1:NoChangeNo
        PopulationNew1(i,:)=Population1(Index(PopulationSize-i+1),:);
        PopulationNew2(i,:)=Population2(Index(PopulationSize-i+1),:);
        PopulationNew3(i,:)=Population3(Index(PopulationSize-i+1),:);
    end   
    %% 轮盘赌法   选择selection     
    for i=(NoChangeNo+1):2:PopulationSize
        [idx1,idx2] = selection(PopulationSize,FitnessSum,PopulationFitness,Index);        
    %%  父母交叉形成子代
        Rate=rand;
        if Rate<=CrossRate
            acr_position = floor(ChromosomeLen*rand+1);   %交叉节点
            [PopulationNew1(i,:),PopulationNew1(i+1,:)]=crossover(acr_position,Population1(idx1,:),Population1(idx2,:));
            [PopulationNew2(i,:),PopulationNew2(i+1,:)]=crossover(acr_position,Population2(idx1,:),Population2(idx2,:));
            [PopulationNew3(i,:),PopulationNew3(i+1,:)]=crossover(acr_position,Population3(idx1,:),Population3(idx2,:));
        end
    end
    %%  变异
    for i=(NoChangeNo+1):PopulationSize 
        PopulationNew1(i,:)=mutation(ChromosomeLen,MutateRate,PopulationNew1(i,:));
        PopulationNew2(i,:)=mutation(ChromosomeLen,MutateRate,PopulationNew2(i,:));
        PopulationNew3(i,:)=mutation(ChromosomeLen,MutateRate,PopulationNew3(i,:));
    end
    for i=1:PopulationSize
        Population1(i,:)=PopulationNew1(i,:);
        Population2(i,:)=PopulationNew2(i,:);
        Population3(i,:)=PopulationNew3(i,:);
    end
    
    K_p(Iter,1)=Elite1(Iter,1);
    K_i(Iter,1)=Elite2(Iter,1);
    K_d(Iter,1)=Elite3(Iter,1);
end
figure(1)
plot(BestFitness,'LineWidth',2);
title('最优个体适应值','fontsize',18);
xlabel('迭代次数');ylabel('适应值');
figure(2)
plot(K_p)
hold on
plot(K_i,'k','LineWidth',3)
plot(K_d,'--r')
title('pid参数优化曲线','fontsize',18);






