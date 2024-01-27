%The code of phase diagram of the  linear delay when $k=3$.
b=1;a=1/2;
f=@(t,x,Z)[a-b*x(1)-Z(1,1)+(x(1)^2)*x(2); b*x(1)-(x(1)^2)*x(2)];
ff=odeset;ff.RelTol=1e-6;
tau=@(t,x)[t-3*(x(1)-1/2)-0.62];
d=[0.6;1.8];
txx=ddesd(f,tau,d,[0,1000],ff);
xxx=0:0.001:1000;
yyy1=interp1(txx.x,txx.y(2,:),xxx,'pchip');
yyy2=interp1(txx.x,txx.y(1,:),xxx,'pchip');
o1=0
hh=size(yyy1)
zzz1=zeros(1,hh(2)-o1);
for i = 1:hh(2)-o1
         zzz1(1,i) = yyy1(i+o1);
end
zzz2=zeros(1,hh(2)-o1);
for i = 1:hh(2)-o1
         zzz2(1,i) = yyy2(i+o1);
end
plot(zzz2(1,:),zzz1(1,:))
xlabel('x(t)');
%xlabel('t');
ylabel('y(t)');
%legend('x(t)','y(t)')