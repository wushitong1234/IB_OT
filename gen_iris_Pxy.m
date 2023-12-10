x = load('iris_data.txt');
x1 = x(:,1); x1_a = min(x1); x1_b = max(x1);
x2 = x(:,2); x2_a = min(x2); x2_b = max(x2);
x3 = x(:,3); x3_a = min(x3); x3_b = max(x3);
x4 = x(:,4); x4_a = min(x4); x4_b = max(x4);
%4-8; 
%2-4.5; 
%1-7; 
%0-2.5; 

step = 1;
n1 = ceil((8-4)/step); 
n2 = ceil((4.5-2)/step); 
n3 = ceil((7-1)/step); 
n4 = ceil((2.5-0)/step);
len_x = n1*n2*n3*n4; 
len_y = 3;
p_xy = zeros(len_x,len_y);
for i=1:length(x1)
    t1 = floor((x(i,1)-4)/step)+1;
    t2 = floor((x(i,2)-2)/step)+1;
    t3 = floor((x(i,3)-1)/step)+1;
    t4 = floor((x(i,4)-0.1)/step)+1;
    nx = t1+(t2-1)*n1+(t3-1)*n1*n2+(t4-1)*n1*n2*n3; 
    ny = x(i,5);
    p_xy(nx,ny) = p_xy(nx,ny)+1;
end
p_xy = p_xy/(length(x1));
% p_xy=p_xy/sum(sum(p_xy));
P_xy=p_xy(sum(p_xy')~=0,:);
q=sum(P_xy)';
p=sum(P_xy')';
s=(P_xy)'*diag(1./p);
s(s<10^(-60))=10^(-60);

M=length(p);
N=29;
K=3;
r0=rand(N,1)+1;
r0=r0/sum(r0);
w0=rand(M,N)+1;
w01=diag(1./sum(w0'));
w0=diag(1./sum(w0'))*w0;z0=s*w0*diag(r0);










