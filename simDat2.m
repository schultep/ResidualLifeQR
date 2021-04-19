
function [QRobj,b1] = simDat2(N,theta,rndEf)
Z=[rand([N,1]),binornd(1,.5*ones([N,1]))];
b0=normrnd(rndEf(1),rndEf(2),[N,1]);
b1=normrnd(rndEf(3),rndEf(4),[N,1]);
delta=zeros([N,1]);
const=exp(theta(2)+Z*theta(4:end)'+theta(1)*b1);

C=5+exprnd(10,[N,1]);
if theta(3)==0
    Y=-log(rand([N,1]))./const;
else
    Y=log(1-theta(3)*log(rand([N,1]))./const)./theta(3);
end

ind=~isreal(Y);
Y(ind)=2*max(C);
ind=C>=Y;
delta(ind)=1;
Y=min([Y,C],[],2);
T=[]; X=[];id=[];

if sum(Y<=0)
    disp('Y is negative')
    Y(Y<=0)
    Y(~isreal(Y))
end

if sum(~isreal(Y))
    disp('Y is not real')
    Y(Y<=0)
    Y(~isreal(Y))
end
r=poissrnd(1.5*Y);
if sum(isnan(r))
    disp('nan')
    disp('min r_ max r')
    min(r)
    max(r)
end
if sum(isinf(r))
    disp('inf')
    disp('min r_ max r')
    min(r)
    max(r)
end
for i=1:N
   t=Y(i)*[0;sort(rand([r(i),1]))]; 
   T=[T;t];
   X=[X;b0(i)+t* b1(i) + normrnd(0,rndEf(5),[r(i)+1,1])];
   id=[id;i*ones([r(i)+1,1])];
end

QRobj=QRdata(id,X,T,Y(id),delta(id),Z(id,:),0,theta,b1);
    
end
