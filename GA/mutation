function y= mutation(n,MutateRate,x)
%n为染色体长度
%MutataRate为变异概率
%x为待变异染色体
rate=rand;    %随机概率判断是否变异
if rate<=MutateRate
    mutate_pos=floor(rand*n+1);    %产生一个变异位置
    x(mutate_pos)=floor(rand+1);  %bianyi 
end
y=x;
end

