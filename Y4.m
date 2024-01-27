%The code of ¡°${\tau }_{2}(x^{*})$ determining criticality for different $k$ for the exponential delay¡±. 
clear all;
a=1/2;
b=1;
w0=((-a^4-b^2+2*a^2*b+1+((-a^4-b^2+2*a^2*b+1)^2+4*a^4)^(1/2))/2)^(1/2);
cwr=(b*(w0)^2)/(a^4+(w0)^2);
swr=((w0)^3+(a^2-b)*a^2*w0)/(a^4+(w0)^2);
c2wr=2*(cwr)^2-1;
s2wr=2*(cwr)*(swr);
t=(c2wr)^2+(s2wr)^2;
a1=2*(w0)-s2wr;
a2=c2wr-b;
M1=(swr-w0)/a^2;
M2=(cwr-b)/a^2;
n=0;
r0=(asin(swr)+2*n*pi)/w0;
r1=0.62;
kvals=linspace(4.5,5,150);
for i=length(kvals):-1:1
k=kvals(i);
f1=k*r0;
f2=k*k*r0;
A=[a1,a2,0,0,-a^2,0;a2,-a1,0,-a^2,0,0;0,0,0,0,0,-a^2;
0,1,0,2*(w0),a^2,0;1,0,0,a^2,-2*(w0),0;0,0,1,0,0,a^2];
b=[a*M2+(1/2)*1.305+(1/2)*w0*swr*f1;a*M1-(1/2)*w0*cwr*f1;
a*M2+(1/2)*1.305+(1/2)*w0*swr*f1;-a*M2-(1/2)*1.305;-a*M1;-a*M2-(1/2)*1.305];
R=A\b;
e1=norm(A*R-b);
A1=[2*M1*M2-r0*(swr),-swr*(w0);(M1)^2-(M2)^2-1+r0*(cwr),w0*(cwr)];
b11=f1*w0*(R(1)*c2wr+R(2)*s2wr)-(1/2)*w0*f1*R(1)*cwr+(3/8)*w0*f2*swr+w0*f1*R(3)*swr+(3/8)*((w0*f1)^2)*cwr+1.305*R(2)+2*1.305*R(3)+1.305*R(1)*M1+a*R(2)*M2+2*a*R(3)*M2+a*R(5)+2*a*R(6)+(3/4)*M2+a*R(1)*M1-a*R(2)*(M1)^2+2*a*R(3)*(M1)^2+R(1)*M1*M2+(1/4)*(M1)^2+a*R(4)*M1-1.305*R(2)*M2-1.305*2*R(3)*M2-a*R(1)*M1*M2-a*R(2)*(M2)^2-2*a*R(3)*(M2)^2-(3/4)*(M2)^2-a*R(5)*M2-2*a*M2*R(6)+(1/2)*(w0)*f1*R(2)*swr;
b12=f1*w0*(R(1)*s2wr-R(2)*c2wr)+(1/2)*w0*f1*R(2)*cwr-(1/8)*w0*f2*cwr-w0*f1*R(3)*cwr+(1/2)*w0*f1*R(1)*swr+(1/8)*((w0*f1)^2)*swr+1.305*R(1)+a*R(1)*M2-a*R(2)*M1+2*a*R(3)*M1+(1/4)*M1+a*R(4)-1.305*R(1)*M2+a*R(2)*M1*M2-2*a*M1*M2*R(3)-a*R(1)*(M2)^2-(1/4)*M1*M2-a*R(4)*M2-1.305*R(2)*M1-2*1.305*R(3)*M1-a*R(1)*(M1)^2-a*R(2)*M1*M2-2*a*R(3)*M2*M1-(3/4)*M1*M2-a*R(5)*M1-2*a*R(6)*M1;
b1=[b11;b12];
b1=[b11;b12];
X=A1\b1;
e2=norm(A1*X-b1);
w2(i)=vpa(X(1),4);
t2(i)=vpa(X(2),4) ;
end
plot(kvals,-t2);
grid on;legend();hold off
positive_indices = t2 > 0;
negative_indices = t2 < 0;
figure;
plot(kvals(positive_indices), t2(positive_indices), 'k-', 'LineWidth', 2);
hold on;
plot(kvals(negative_indices), t2(negative_indices), 'k--', 'LineWidth', 2);
grid on;
legend('supercritical', 'subcritical');
hold off;
xlabel('k');
ylabel('{\tau }_{2}(x^{*})');
title('{\tau }_{2}(x^{*}) determining criticality for different k for the linear delay')