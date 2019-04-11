function [idx1,idx2] = selection(PopulationSize,FitnessSum,PopulationFitness,Index)
%轮盘赌法，确定要交叉的父母编号
        %确定要交叉的父亲序号
        m1=0;
        r1=rand*FitnessSum;
        for k=1:PopulationSize
            m1=m1+PopulationFitness(k);
            if r1<=m1
                idx1=Index(k);
                break;
            end
        end
         %确定要交叉的母亲序号
        m2=0;
        r2=rand*FitnessSum;
        for k=1:PopulationSize
            m2=m2+PopulationFitness(k);
            if r2<=m2
                idx2=Index(k);
                break;
            end
        end
    
end

