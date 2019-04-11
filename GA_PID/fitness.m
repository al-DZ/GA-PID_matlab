function fit=fitness(Kp,Ki,Kd)

assignin('base','Kp',Kp);
assignin('base','Ki',Ki);
assignin('base','Kd',Kd);
[~,~,y_out]=sim('thetapid',[0,50]);
z=(y_out(end,1));   %%表示 取这个 矩阵的第一列最后一行的数据。
t0=yy.time;
y=yy.signals.values;
% res=yy;
% ymax=find(y>=1);
% tr=t0(ymax(1));%计算上升时间tr
[ym,tp]=max(y);
% tp=t0(tp);%计算峰值时间tp
Mp=ym-1;%计算超调量Mp
s0=length(t0);
while y(s0)>0.98&&y(s0)<1.02
    s0=s0-1;
end
ts=t0(s0);%计算调整时间ts
% error=sum(abs(yy1.signals.values-9));
w=[0.8 0.1 0.1];
zhibiao=[z ts Mp];
fit=1/(zhibiao*w');
