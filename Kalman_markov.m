function [Pr1,y1,P1,Pr2,y2,P2,P1_,y1_]=Kalman_markov(r,Pr1,y1,P1,Pr2,y2,P2,CSADt,sigv,sige,p11,p22)

Pr={Pr1,Pr2};
P={P1,P2};
y={y1,y2};
%滤波计算
y_=cellmat(2,2);
P_=cellmat(2,2);
ep_=cellmat(2,2);
P1_=cellmat(2,2);%%P(t|t-1)(i,j)
y1_=cellmat(2,2);%y(t|t-1)(i,j)
for i=1:2;
    for j=1:2;
        [yij,pij,epij,y1ij,P1ij]=Kalman(r,y{i},P{i},j,CSADt,sigv,sige);
        y_{i,j}=yij;
        P_{i,j}=pij;
        ep_{i,j}=epij;
        y1_{i,j}=y1ij;
        P1_{i,j}=P1ij;
    end
end

%计算后验概率分布
Prt_=zeros(2,2);
gauss=@(x,mu,P)1/sqrt(2*pi)/sqrt(det(P))*exp(-x(:)'*(P\x(:))/2);
ptrans=[p11,1-p11;1-p22,p22];
[~,~,H,~]=getmodel(sigv,sige,r,1);
sumPrt=0;
for i=1:2;
    for j=1:2;
        Vj=sigv(j);pij=ptrans(i,j);
        Prt_(i,j)=gauss(ep_{i,j},0,H*P_{i,j}*H'+Vj)*pij*Pr{i};
        sumPrt=sumPrt+Prt_(i,j);
    end
end
for i=1:2;
    for j=1:2;
        Prt_(i,j)=Prt_(i,j)/sumPrt;
    end
end

%输出
Pr_=zeros(2,1);
ynew_=cellmat(2,1);
Pnew_=cellmat(2,1);
for j=1:2;
    Pr_(j)=sum(Prt_(:,j));
    ynew_{j}=(y_{1,j}*Prt_(1,j)+y_{2,j}*Prt_(2,j))/(Prt_(1,j)+Prt_(2,j));
    ep1=y_{1,j}-ynew_{j};li1=P_{1,j}+ep1(:)*ep1(:)';
    ep2=y_{2,j}-ynew_{j};li2=P_{2,j}+ep2(:)*ep2(:)';
    Pnew_{j}=(li1*Prt_(1,j)+li2*Prt_(2,j))/(Prt_(1,j)+Prt_(2,j));
end
Pr1=Pr_(1);Pr2=Pr_(2);
y1=ynew_{1};y2=ynew_{2};
P1=Pnew_{1};P2=Pnew_{2};


