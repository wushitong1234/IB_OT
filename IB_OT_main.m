%% Gaussian source
M=100;N=100;K=100;snr=1;
%r0=1/N*ones(N,1);
r0=rand(N,1);r0=r0+1;r0=r0/sum(r0);
h=20/K;y=((-10):h:(10-h))';x=((-10):h:(10-h))';
q=h*1/sqrt(2*pi)*exp(-y.^2/2);q(q<10^(-60))=10^(-60);
s1=diag(q)*(1/sqrt(2*pi*1/snr)*exp(-snr*(y*ones(1,M)-ones(K,1)*x').^2/2))*h;
p=sum(s1)';s=s1*diag(1./p);s(s<10^(-60))=10^(-60);
p=p/sum(p);q=s*p;
w0=rand(M,N)+1;w0=w0*diag(1./sum(w0));z0=s*w0*diag(r0);
I=0.1;

%% Bernoulli source
M=2;N=2;K=2;
r0=rand(N,1)+1;r0=r0/sum(r0);
u=0.3;
p=[0.5;0.5];I=log(2)-(-u*log(u)-(1-u)*log(1-u));e=0.15;
s=[e,1-e;1-e,e];v=(u-e)/(1-2*e);q=s*p;
w0=[0.99,0.01;0.01,0.99];
z0=[0.01,0.49;0.49,0.01];
% w0=rand(M,N);w0=w0*diag(1./sum(w0));z0=s*w0*diag(r0);


%% The GAS
time_tol_AS=0;number=1;
for i1=1:number
time_AS=tic;
g=1/N*ones(N,1);f=1/M*ones(M,1);R=zeros(M,N);
r=r0;
zeta=1;
w=w0;z=z0;
lambda=-zeta*ones(K,N);G=0;
T=I+sum(q.*log(q));
R1=-s'*lambda;R2=s'*log(z);
value=0;delta=1;

num=0;
while (abs(delta)>10^(-6))
    num=num+1;
    R1=-s'*lambda;R2=s'*log(z);
    R=R1+zeta*(R2-ones(M,1)*max(R2));
    g=g-log(sum(exp(R+f*ones(1,N)+ones(M,1)*g')))';
    f=log(p)+f-log(exp(R+f*ones(1,N)+ones(M,1)*g')*r);
    r(r<10^(-60))=10^(-60);
    c=T+sum(r.*log(r));     
    i=0;
    G=sum(sum((s*(exp(R1+zeta*(R2-ones(M,1)*max(R2))+f*ones(1,N)+ones(M,1)*g'))*diag(r)).*log(s*w*diag(r))))-c;
    while ((i<20)&&(abs(G)>10^(-12)))
        i=i+1;
        zeta=zeta-(G)/(sum(sum((s*(exp(R1+zeta*(R2-ones(M,1)*max(R2))+f*ones(1,N)+ones(M,1)*g').*(R2-ones(M,1)*max(R2)))*diag(r)).*(log(s*w*diag(r))))));
        G=sum(sum((s*(exp(R1+zeta*(R2-ones(M,1)*max(R2))+f*ones(1,N)+ones(M,1)*g'))*diag(r)).*log(s*w*diag(r))))-c;
    end
    w=exp(R+f*ones(1,N)+ones(M,1)*g');
    
    lambda=-zeta*ones(K,N);
    
    z=s*w*diag(r);
    z(z<10^(-60))=10^(-60);
    
    ww0=w;ww0(ww0<10^(-30))=1;
    a=sum((s*w).*log(z))'-sum(ww0.*log(ww0))'/zeta-sum(diag(-f-0.5)*w)'/zeta-sum(w*diag(-g+zeta*max(R2)'-0.5))'/zeta+(-g+zeta*max(R2)'-0.5)/zeta-sum((s*w).*lambda)'/zeta-1;
    r=exp(a-max(a)); 
    r=r/sum(r);
    
    ww=w*diag(r);ww(ww<10^(-30))=1;
    value=sum(sum(ww.*log(ww)))-sum(r.*log(r))-sum(p.*log(p));
    %delta=value-(-0.5*log(((1+snr)*exp(-2*I)-1)/snr));
    delta=value-(log(2)+v*log(v)+(1-v)*log(1-v)); 
end
 
time_AS=toc(time_AS)
time_tol_AS=time_tol_AS+time_AS;
end

time_avg_AS=time_tol_AS/number
