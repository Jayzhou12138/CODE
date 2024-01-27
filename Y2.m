%The code of ¡°Lyapunov coefficients determining criticality for different $k$¡±.
clear;format compact
parnames={'a','b','k','beta'}; % parameter names
cind=[parnames;num2cell(1:length(parnames))]; % convert parameter names to struct
ip=struct(cind{:});
f=@(x,a,b)[a-b.*x(1,1,:)-x(1,2,:)+x(1,1,:).^2.*x(2,1,:);...
b.*x(1,1,:)-x(1,1,:).^2.*x(2,1,:)]; %r.h.s.
tauexp=@(x,k,beta)beta.*exp(k.*(x(1,1,:)-1/2));% exponential delay
taulin=@(x,k,beta)k.*(x(1,1,:)-1/2)+beta; % linear delay
f=@(tau)set_funcs('sys_rhs',@(x,p)f(x,p(ip.a),p(ip.b)),...
'sys_tau',@(it,x,p)tau(x,p(ip.k),p(ip.beta)),...
'sys_ntau',@()1,'x_vectorized',true);
problems={f(taulin),f(tauexp)};
names={'linear delay','exponential delay'}; % define r.h.s for both delay types
kvals=linspace(2,5,10); % range of values for k
par([ip.a,ip.b,ip.beta])=[0.5,1,0.61]; 
for i=length(problems):-1:1
for k=length(kvals):-1:1
par(ip.k)=kvals(k);
[eqbr,suc]=SetupStst(problems{i},'parameter',par,'x',[0.5;2]);
[hbr(i,k),suc(i,k)]=SetupHopf(problems{i},eqbr,1,'contpar',ip.beta); 
L1(i,k)=HopfLyapunovCoefficients(problems{i},hbr(i,k)); 
fprintf('problem: %s, k=%g, L1=%g\n',names{i},kvals(i),L1(i,k));
end
end
%%
figure(1);clf;hold on
for i=1:length(names)
plot(kvals,L1(i,:),'DisplayName',names{i});
end
title('Lyapunov coefficients determining criticality for different k')
xlabel('k');ylabel('L1');grid on;legend();hold off