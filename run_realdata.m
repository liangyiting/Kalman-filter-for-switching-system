function []=erere()
global zx_ratio;
zx_ratio=0.05;%在画置信区间的时候只波动 zx_ratio*方差
close all
%股市CASD数据模拟所需参数
%线性-马尔可夫切换-状态空间模型参数
parts={'All Stocks'};
if 1
    data=xlsread('模拟数据1.xlsx',1,'D4:E1219');%读取真实数据
    rr=data(:,1);%Rm回报率时间序列
    zz=data(:,2);%CSADt时间序列
%     rr=xlsread('111.xlsx',1,'M2:M169');
%     zz=xlsread('111.xlsx',1,'W2:W169');
else
    data=load('simulatedata.txt');%读取模拟数据
    rr=data(:,1);%Rm回报率时间序列
    zz=data(:,2);%CSADt时间序列
    
end

T=numel(rr)/1;
rr=rr(1:T,:);
zz=zz(1:T,:);
tt=1:T;
tt=(1:T)/365+2004;%
figure(1);
a=2;b=2;
subplot(a,b,1);plot(tt,zz),title(strcat('CASD of--',' ',parts{1}));
names={'\alpha_1','\alpha_2','herding coefficient (\alpha_3)'};
disp('样本数目');disp(T)
%%
y0=[0 0 0.1]';
maxiterations=12;%优化最大迭代次数
[sigv,sige,p,y0]=estimate(zz,rr,maxiterations);
%%
if 0
    %卡尔曼滤波器测试，假设St已知
    Pt0=0*eye(3);
    [logL,yyk,yykup,yyklow]=objfun_kalman(ss,zz,rr,sigv,sige,y0,Pt0);
    figure(1);for i=1:3;subplot(a,b,1+i);hold on;plot(tt,yyk(:,i),'r');title(names{i});if i==1;legend('实际值','滤波值');end;end;
    figure(4);
    hold on;
    maxyy=max(yykup(:,i))+0.2;minyy=min(yyklow(:,i))-0.2;
    i=3;fill([tt(:);flipud(tt(:))],[yyklow(:,i);flipud(yykup(:,i))],[0.6 0.6 0.6]);hold on;plot(tt,yyk(:,i),'b--','linewidth',1.0);
    %plot(tt,yy(:,i),'b--','linewidth',1.0);
    %plot(tt,(cpr1<0.5)-0.5,'b-');
    %figure;li=ss(:)-1;fill([tt(1);tt(:);tt(end)],[minyy;li*maxyy+(1-li)*minyy;minyy],[1 1 0.7]);
    legend('St','置信区间(95%)','herding coefficient估计值');title(parts{1});
    box on;
    if 0;axis([-inf,inf,minyy+0.1,maxyy-0.1]);end
    hold off;
end
%%
%卡尔曼滤波-马尔可夫切换
%y0=[0 0 0.1]';
P1=11*eye(3);P2=11*eye(3);Pr1=0.9;
[logL,cpr1,cP1x,cy1x,cPx,cyx,cPr,yyk1,Pr1,y1,P1,Pr2,y2,P2,yyk1up,yyk1low]=objfun_kalman_markov(zz,rr,p,sigv,sige,y0,Pr1,P1,P2);
m=60;
figure(1);for i=1:3;subplot(a,b,1+i);plot(tt(m:end),yyk1(m:end,i),'r');title(names{i});if i==1;legend('实际值','滤波值');end;end;
figure(4);
maxyy=max(yyk1up(m:end,i))*1;minyy=min(yyk1low(m:end,i))*1;
li=(cpr1(:)<0.5);fill([tt(1);tt(:);tt(end)],[100;li*(-20)+(1-li)*20;20],[0.8 0.8 0.8]);hold on;
i=3;fill([tt(m:end)';flipud(tt(m:end)')],[yyk1low(m:end,i);flipud(yyk1up(m:end,i))],[0.6 0.6 0.6]);
hold on;plot(tt(m:end),yyk1(m:end,i),'b--','linewidth',1.0);
legend('St','置信区间(95%)','herding coefficient估计值');title(parts{1});
axis([-inf,inf,minyy-0.1*abs(maxyy-minyy),maxyy+0.1*abs(maxyy-minyy)]);
box on;

%%
%平滑
p11=p(1);p22=p(2);
yyk2=zeros(T,3);
yyk2(T,:)=yyk1(end,:);
yyk2up=yyk2;
yyk2low=yyk2;
cpr1=zeros(T,1);
cpr2=zeros(T,2);
cer=zeros(T,1);
zz2=zeros(T,1);
r=rr(T);H=[1,abs(r),r^2];yt=yyk1(T,:);zt=H*yt(:);zz2(T)=zt;
for t=T:-1:2;
    Prx=cPr{t-1};%这里对应的是t-1时刻！！！,=Pr(S(t-1)=i|t-1);
    Px=cPx{t-1};%这里对应的是t-1时刻！！！,=P(S(t-1)=i|t-1)
    P1x=cP1x{t};%这里对应的是t时刻！！！，=P(S(t-1)=i,S(t)=j|t-1)
    yx=cyx{t-1};%这里对应的也是t-1时刻！！！，=y(t-1)(S(t-1)=i|t-1);
    y1x=cy1x{t};%这里对应的是t时刻！！！，=y(t)(S(t-1)=i,S(t)=j|t-1);
    [Pr1,y1,P1,Pr2,y2,P2,yt,Pt]=Smooth(Pr1,y1,P1,Pr2,y2,P2,Prx,Px,P1x,yx,y1x,p11,p22);
    yyk2(t-1,:)=yt;
    r=rr(t-1);H=[1,abs(r),r^2];zt=H*yt(:);
    zz2(t-1,:)=zt;
    cer(t-1,:)=zz(t-1)-zt;
    yyk2up(t-1,:)=yt(:)+sqrt(diag(Pt))*1.96;
    yyk2low(t-1,:)=yt(:)-sqrt(diag(Pt))*1.96;
    cpr1(t-1,:)=Pr1;
    cpr2(t-1,:)=Pr2;
end
figure(1);for i=1:3;subplot(a,b,1+i);hold on;plot(tt(m:end),yyk2(m:end,i),'k');title(names{i});if i==1;legend('滤波值','平滑值');end;end
disp('平均绝对误差')
figure(2);
a=4;b=1;
subplot(a,b,2);plot(tt,cpr1);hold on;plot(tt,cpr2);legend('Pr(状态St=0)','Pr(状态St=1)');hold off;
subplot(a,b,4);plot(tt,cer);legend('equation残差');hold off;
subplot(a,b,3);hold on;plot(tt,[zz,zz2]);legend('CSADt','CSAD_t-平滑值预测');
subplot(a,b,1);plot(tt,rr);legend('回报率(Rmt)')
herding=yyk2(:,3);
disp('save(herding.txt)');
save('herding.txt','herding','-ascii');
%%
function [Pr1,y1,P1,Pr2,y2,P2,yt,Pt]=Smooth(Pr1,y1,P1,Pr2,y2,P2,Prx,Px,P1x,yx,y1x,p11,p22)
%Prx为t时刻考虑前t个观测的先验概率密度分布=P(S(t)=i|t);
%Pr为考虑T个观测的后验概率分布=P(S(t+1)=i|T)
%Prx=Pr(S(t)=i|t),P1x=P(t+1|t)(i,j),Px=P(t|t)(i),y1x=y(t+1|t)(i,j),yx=y(t|t)(i)
[~,~,~,F]=getmodel(0,0,1,1);
PrT=zeros(2,2);%Pr(S(t)=i,S(t+1)=j|T)
Pr=[Pr1,Pr2];
ptrans=[p11,1-p11;1-p22,p22];
for i=1:2;
    for j=1:2;
        PrT(i,j)=Pr(j)*ptrans(i,j)*Prx(i)/(Prx(1)*ptrans(1,j)+Prx(2)*ptrans(2,j));
    end
end

Prnew=zeros(2,1);%Pr(S(t)=i|T)
for i=1:2;
    Prnew(i)=PrT(i,1)+PrT(i,2);
end
Pr1=Prnew(1);Pr2=Prnew(2);

yT1=cellmat(2,2);%y(t)(S(t)=i,S(t+1)=j|T)
yT={y1,y2};%y(t)(S(t+1)=j|T)
PT1=cellmat(2,2);
PT={P1,P2};
for i=1:2;
    for j=1:2
        if 1
            pij=Px{i}*(F)'/P1x{i,j};%貌似错了
            if cond(P1x{i,j})>1e20;
                1;
            end
            yT1{i,j}=yx{i}+pij*(yT{j}-y1x{i,j});
            PT1{i,j}=Px{i}+pij*(PT{j}-P1x{i,j})*pij';
        else
            yT1{i,j}=y1x{i,j};
            PT1{i,j}=P1x{i,j};
        end
    end
end

yTnew=cellmat(2,1);
for i=1:2;
    yTnew{i}=yT1{i,1}*PrT(i,1)/(PrT(i,1)+PrT(i,2))...
        +yT1{i,2}*PrT(i,2)/(PrT(i,1)+PrT(i,2));
end
y1=yTnew{1};y2=yTnew{2};
yt=yTnew{1}*(PrT(1,1)+PrT(1,2))+yTnew{2}*(PrT(2,1)+PrT(2,2));

PTnew=cellmat(2,1);
for i=1:2;
    j=1;epij=yT1{i,j}-yTnew{i};
    PTnew{i}=(PT1{i,j}+epij(:)*epij(:)')*PrT(i,j)/(PrT(i,1)+PrT(i,2));
    j=2;epij=yT1{i,j}-yTnew{i};
    PTnew{i}=PTnew{i}+(PT1{i,j}+epij(:)*epij(:)')*PrT(i,j)/(PrT(i,1)+PrT(i,2));
end
P1=PTnew{1};P2=PTnew{2};
Pt=P1*Pr1+P2*Pr2;

%%
function [yt,Pt,ept,y1,P1,dL]=Kalman(r,yt0,Pt0,St,CSADt,sigv,sige)
[W,V,H,F]=getmodel(sigv,sige,r,St);
y1=F*yt0;
P1=F*Pt0*F'+W;
z1=H*y1;
ep1=CSADt-z1;
yt=y1+P1*H'*((H*P1*H'+V)\ep1);
Pt=(eye(3)-P1*H'*((H*P1*H'+V)\H))*P1;
ept=ep1;
%gauss=@(x,mu,P)(2*pi)^(-1/2)*(det(P))^(-1/2)*exp(-x(:)'*(P\x(:))/2);
%dL=log(gauss(ept,0,H*P1*H'+V));
x=ept;P=H*P1*H'+V;
dL=-1/2*log(2*pi)-1/2*log(det(P))-x(:)'*(P\x(:))/2;
%%
function [Pr1,y1,P1,Pr2,y2,P2,P1_,y1_,dL,yt,Pt]=Kalman_markov(r,Pr1,y1,P1,Pr2,y2,P2,CSADt,sigv,sige,p11,p22)
Pr={Pr1,Pr2};%Pr=Pr(S(t-1)=i|t-1),
P={P1,P2};%P=P(S(t-1)=i|t-1);
y={y1,y2};%y=y(S(t-1)=i|t-1);
%滤波计算
y_=cellmat(2,2);%y_=y(S(t-1)=i,S(t)=j|t);
P_=cellmat(2,2);%P_=P(S(t-1)=i,S(t)=j|t);
ep_=cellmat(2,2);%ep_=(CASD-z(t))(S(t-1)=i,S(t)=j|t-1);
P1_=cellmat(2,2);%%P(t)(S(t-1)=i,S(t)=j|t-1)
y1_=cellmat(2,2);%y(t)(S(t-1)=i,S(t)=j|t-1)
dL_=zeros(2,2);
for i=1:2;
    for j=1:2;
        [yij,pij,epij,y1ij,P1ij,dlij]=Kalman(r,y{i},P{i},j,CSADt,sigv,sige);
        y_{i,j}=yij;
        P_{i,j}=pij;
        ep_{i,j}=epij;
        y1_{i,j}=y1ij;
        P1_{i,j}=P1ij;
        dL_(i,j)=dlij;
    end
end

%计算后验概率分布
Prt_=zeros(2,2);%Pr(S(t-1)=i,S(t)=j|t)
gauss=@(x,mu,P)(2*pi)^(-1/2)*(det(P))^(-1/2)*exp(-x(:)'*(P\x(:))/2);
ptrans=[p11,1-p11;1-p22,p22];
[~,~,H,~]=getmodel(sigv,sige,r,1);
sumPrt=0;
for i=1:2;
    for j=1:2;
        Vj=sigv(j);pij=ptrans(i,j);
        if 0
            Prt_(i,j)=gauss(ep_{i,j},0,H*P1_{i,j}*H'+Vj)*pij*Pr{i};%Pr(S(t-1)=i,S(t)=j,z(t)=CSAD(t)|t-1)
        else
            Prt_(i,j)=exp(dL_(i,j))*pij*Pr{i};
        end
        sumPrt=sumPrt+Prt_(i,j);
    end
end
for i=1:2;
    for j=1:2;
        Prt_(i,j)=Prt_(i,j)/sumPrt;%Pr(S(t-1)=i,S(t)=j|t)
    end
end

%输出
Pr_=zeros(2,1);%Pr(S(t)=j|t)
ynew_=cellmat(2,1);%y(t)(S(t)=j|t)
Pnew_=cellmat(2,1);%P(t)(S(t)=j|t)
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
dL=log(sumPrt);
yt=y1*Pr1+y2*Pr2;
Pt=P1*Pr1+P2*Pr2;


function [W,V,H,F]=getmodel(sigv,sige,r,St)
W=diag(sigv.^2);
V=sige(St)^2;
H=[1,abs(r),r^2];
F=eye(3);

%%
function [logL,cpr1,cP1x,cy1x,cPx,cyx,cPr,yyk1,Pr1,y1,P1,Pr2,y2,P2,yyk1up,yyk1low]=objfun_kalman_markov(zz,rr,p,sigv,sige,y0,Pr1,P1,P2)
global zx_ratio;
T=size(zz,1);
p11=p(1);p22=p(2);
y1=y0;y2=y0;
Pr2=1-Pr1;
yyk1=zeros(T,3);
cpr1=zeros(T,1);
cP1x=cellmat(T,1);
cy1x=cellmat(T,1);
cPx=cellmat(T,1);
cyx=cellmat(T,1);
cPr=cellmat(T,1);
logL=0;%极大似然函数
yyk1up=yyk1;
yyk1low=yyk1;
for t=1:T;
    %St=ss(t);
    CSADt=zz(t,:);r=rr(t);
    [Pr1,y1,P1,Pr2,y2,P2,P1x,y1x,dL,yt,Pt]=Kalman_markov(r,Pr1,y1,P1,Pr2,y2,P2,CSADt,sigv,sige,p11,p22);
    logL=logL+dL;
    yyk1(t,:)=yt;
    cpr1(t)=Pr1;
    cP1x{t}=P1x;
    cy1x{t}=y1x;
    cPx{t}={P1,P2};
    cyx{t}={y1,y2};
    cPr{t}=[Pr1,Pr2];
    yyk1up(t,:)=yt(:)+sqrt(diag(Pt))*zx_ratio;
    yyk1low(t,:)=yt(:)-sqrt(diag(Pt))*zx_ratio;
end


function [logL,yyk,yykup,yyklow]=objfun_kalman(ss,zz,rr,sigv,sige,y0,Pt0)
global zx_ratio;
T=numel(ss);
yt0=y0;
yyk=zeros(T,3);
logL=0;
yykup=yyk;
yyklow=yyk;
for t=1:T;
    St=ss(t);CSADt=zz(t,:);r=rr(t);
    %St=2;
    [yt,Pt,~,~,~,dL]=Kalman(r,yt0,Pt0,St,CSADt,sigv,sige);
    yt0=yt;Pt0=Pt;
    yyk(t,:)=yt;
    logL=logL+dL;
    yykup(t,:)=yt(:)+sqrt(diag(Pt))*zx_ratio;
    yyklow(t,:)=yt(:)-sqrt(diag(Pt))*zx_ratio;
    
end

function [ss,zz,rr,yy,tt]=simulate(p,sigv,sige,rmin,rmax,meanr,stdr,T);
%F=eye(3);
y0=[0 0 0.1]';S0=2;
rr=min(rmax,max(rmin,meanr+randn(T,1)*stdr));%产生模拟的股市回报率数据
yy=zeros(T,3);yt=y0;St=S0;
zz=zeros(T,1);
ss=zeros(T,1);
for t=1:T;
    r=rr(t);
    %H=[1,abs(r),r^2];
    [~,~,H,F]=getmodel(sigv,sige,r,St);
    wt=sigv.*randn(size(sigv));
    et=sige.*rand(size(sige));
    yt1=F*yt+wt;
    pt=p(St);
    li=rand();
    if 1;
        if li>pt; St=3-St;end
    else
        St=1;
    end
    ss(t)=St;
    zt1=H*yt1+et(St);
    yy(t,:)=yt1;
    zz(t,:)=zt1;
    yt=yt1;
end
tt=1:T;


function [sigv,sige,p,y0]=estimate(zz,rr,maxiteration,ss,sigv,sige,p,y0);
disp('开始极大似然估计，估计参数sigv,sige,p');
if nargin<3;
    maxiteration=10;%优化程序的最大迭代次数
    ss=ones(size(rr));
    sigv=zeros(3,1);sige=zeros(2,1);p=zeros(2,1);y0=zeros(3,1);%没什么意义，可以随便定。
end
%%
if 0
    %已知切换状态St，基于一般卡尔曼滤波方法进行参数估计 并比较不同的初始Pt0设置对参数估计结果的影响
    if 0
        %Pt0最好不要都是0，否则会出现对奇异矩阵求逆的情况
        y0=[0 0 0.1]';Pt0=0.1*eye(3);
        fun=@(sig)-objfun_kalman(ss,zz,rr,sig(1:3),sig(4:5),y0,Pt0);
        pmin=zeros(5,1);pmax=ones(5,1);ptrue=[sigv;sige];
        p0=pmin+rand(size(pmin)).*(pmax-pmin);
        %argmin fun s.t. sige(1)-sige(2)<0  sigv>0,sige>0;
        opt=optimset('display','iter');
        popt=fmincon(fun,p0,[0 0 0 1 -1],0,[],[],pmin,pmax,[],opt);
        [popt(:),ptrue],
    end
    if 0
        %在不同的初始Pt下估计 好像也没有多少效果
        y0=[0 0 0.1]';Pt0=10*eye(3);
        fun=@(sig)-objfun_kalman(ss,zz,rr,sig(1:3),sig(4:5),y0,Pt0);
        pmin=zeros(5,1);pmax=ones(5,1);ptrue=[sigv;sige];
        p0=pmin+rand(size(pmin)).*(pmax-pmin);
        %argmin fun s.t. sige(1)-sige(2)<0  sigv>0,sige>0;
        opt=optimset('display','iter');
        popt=fmincon(fun,p0,[0 0 0 1 -1],0,[],[],pmin,pmax,[],opt);
        [popt(:),ptrue],
    end
    if 1
        %同时对初始点进行估计  试验表明这样做好像不会有改善
        fun=@(sig)-objfun_kalman(ss,zz,rr,sig(1:3),sig(4:5),y0,diag(sig(6:8)));
        pmin=zeros(8,1);pmax=ones(8,1);ptrue=[sigv;sige;[0;0;0]];
        p0=pmin+rand(size(pmin)).*(pmax-pmin);
        %argmin fun s.t. sige(1)-sige(2)<0  sigv>0,sige>0;
        opt=optimset('display','iter');
        popt1=fmincon(fun,p0,[0 0 0 1 -1 0 0 0],0,[],[],pmin,pmax,[],opt);
        [popt1(:),ptrue],
    end
end
%%
if 1
    %切换状态St未知，基于卡尔曼-马尔可夫切换滤波进行参数估计
    if 0
        %1,在切换概率p已知的情况下估计
        y0=[0 0 0.1]';P1=0*eye(3);P2=0*eye(3);Pr1=1;
        fun=@(pk)-objfun_kalman_markov(zz,rr,p,pk(1:3),pk(4:5),y0,Pr1,P1,P2);
        pmin=zeros(5,1);pmax=ones(5,1);ptrue=[sigv;sige];
        p0=pmin+rand(size(pmin)).*(pmax-pmin);
        %argmin fun s.t. sige(1)-sige(2)<0  sigv>0,sige>0;
        opt=optimset('display','iter');
        popt=fmincon(fun,p0,[0 0 0 1 -1],0,[],[],pmin,pmax,[],opt);
        [popt(:),ptrue],
    end
    
    if 1
        %2，在切换概率p未知的条件下估计
        y0=[0 0 0.1]';P1=1*eye(3);P2=1*eye(3);Pr1=1;
        fun=@(pk)-objfun_kalman_markov(zz,rr,pk(6:7),pk(1:3),pk(4:5),y0,Pr1,P1,P2);
        pmin=zeros(7,1);pmax=ones(7,1);%ptrue=[sigv;sige;p;y0];
        pmax(1)=0.1;pmax(2:3)=0.01;
        pmax(4)=0.5;pmax(5)=1.4;
        pmax(6:7)=1;pmin(6)=0.9;pmin(7)=0.75;
        p0=pmin+rand(size(pmin)).*(pmax-pmin);
        %argmin fun s.t. sige(1)-sige(2)<0  sigv>0,sige>0;
        opt=optimset('display','iter','maxiter',maxiteration);
        popt=fmincon(fun,p0,[0 0 0 1 -1 0 0],0,[],[],pmin,pmax,[],opt);
        disp('参数估计值');
        disp(popt);
        sigv=popt(1:3);sige=popt(4:5);p=popt(6:7);
        disp('alpha波动协方差：sigv=');disp(sigv);
        disp('CSAD输出误差协方差：sige=');disp(sige);
        disp('切换概率：p=[p11,p22]=');disp(p);
        disp(strcat('最优似然对数值=',num2str(-fun(popt))));
    end
    if 0
        %2，在切换概率p未知,初始的y0未知的条件下估计
        y0=[0 0 0.1]';P1=1*eye(3);P2=1*eye(3);Pr1=1;
        fun=@(pk)-objfun_kalman_markov(zz,rr,pk(6:7),pk(1:3),pk(4:5),pk(8:10),Pr1,P1,P2);
        pmin=zeros(10,1);pmax=ones(10,1);%ptrue=[sigv;sige;p;y0];
        pmax(6:7)=1;
        p0=pmin+rand(size(pmin)).*(pmax-pmin);
        %argmin fun s.t. sige(1)-sige(2)<0  sigv>0,sige>0;
        opt=optimset('display','iter','maxiter',maxiteration);
        popt=fmincon(fun,p0,[0 0 0 1 -1 0 0 0 0 0],0,[],[],pmin,pmax,[],opt);
        disp('参数估计值');
        disp(popt);
        sigv=popt(1:3);sige=popt(4:5);p=popt(6:7);y0=popt(8:10);
        disp('sigv');disp(sigv);
        disp('sige');disp(sige);
        disp('p');disp(p);
        disp('y0');disp(y0);
        disp(strcat('最优似然对数值=',num2str(-fun(popt))));
    end
end
