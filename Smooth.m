function [Pr1,y1,P1,Pr2,y2,P2,yt]=Smooth(Pr1,y1,P1,Pr2,y2,P2,Prx,Px,P1x,yx,y1x,p11,p22)
%Prx为t时刻考虑前t个观测的先验概率密度分布=P(S(t)=i|t);
%Pr为考虑T个观测的后验概率分布=P(S(t+1)=i|T)
%Prx=Pr(S(t)=i|t),P1x=P(t+1|t)(i,j),Px=P(t|t)(i),y1x=y(t+1|t)(i,j),yx=y(t|t)(i)

[~,~,~,F]=getmodel(0,0,1,1);

PrT=zeros(2,2);
Pr=[Pr1,Pr2];
ptrans=[p11,1-p11;1-p22,p22];
for i=1:2;
    for j=1:2;
        PrT(i,j)=Pr(j)*ptrans(i,j)*Prx(i)/(Prx(1)*ptrans(1,j)+Prx(2)*ptrans(2,j));
    end
end

Prnew=zeros(2,1);
for i=1:2;
    Prnew(i)=PrT(i,1)+PrT(i,2);
end
Pr1=Prnew(1);Pr2=Prnew(2);

yT1=cellmat(2,2);
yT={y1,y2};
PT1=cellmat(2,2);
PT={P1,P2};
for i=1:2;
    for j=1:2
        pij=Px{i}*F'/P1x{i,j};
        yT1{i,j}=yx{i}+pij*(yT{j}-y1x{i,j});
        PT1{i,j}=Px{i}+pij*(PT{j}-P1x{i,j})*pij';
    end
end

yTnew=cellmat(2,1);
for i=1:2;
    yTnew{i}=yT1{i,1}*PrT(i,1)/(PrT(i,1)+PrT(i,2))...
        +yT1{i,2}*PrT(i,2)/(PrT(i,1)+PrT(i,2));
end
y1=yTnew{1};y2=yTnew{2};
yt=yTnew{1}*(PrT(i,1)+PrT(i,2))+yTnew{2}*(PrT(i,1)+PrT(i,2));

PTnew=cellmat(2,1);
for i=1:2;
    j=1;epij=yT1{i,j}-yTnew{i};
    PTnew{i}=(PT1{i,j}+epij(:)*epij(:)')*PrT(i,j)/(PrT(i,1)+PrT(i,2));
    j=2;epij=yT1{i,j}-yTnew{i};
    PTnew{i}=PTnew{i}+(PT1{i,j}+epij(:)*epij(:)')*PrT(i,j)/(PrT(i,1)+PrT(i,2));
end
P1=PTnew{1};P2=PTnew{2};

        
