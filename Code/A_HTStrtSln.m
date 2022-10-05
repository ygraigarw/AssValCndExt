function p0=HTStrtSln(X,Y,b0)
%HT Starting solution estimates a, mu | b=0 and Z=0 then s | on the rest
% X main margin
% Y associated margin

sX=X;
P=0;

sXpb=X.^b0;
P=blkdiag(P,0);

Q=[sX,sXpb];
p=(Q'*Q+P)\Q'*Y;

a0=min(p(1),1-1e-6);
if any(a0(:)<-0.9) || any(a0(:)>0.9)
    a0(a0<-0.9)=-0.9;
    a0(a0>0.9)=0.9;
end

m0=p(2,:);

%stationary in b TODO think about nono stationary b0 case here
R=(Y-Q*p)./(X.^b0);

s0=var(R);

p0=[a0;b0;m0;s0];

end