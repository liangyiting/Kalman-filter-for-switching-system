function []=erere()
close all
%股市CASD数据模拟所需参数
stdr=[0.76 0.87 1 0.76 0.78 0.61 0.77 1.71]'*1e-0;
meanr=[1.61 1.58 1.68 1.34 1.62 1.11 1.54 1.74]*1e-0;
rmin=[0.11 0.04 0.31 0.05 0.06 0.05 0.11 0.01]*1e-0;
rmax=[6.72 8.10 17.5 9.18 7.36 5.47 6.77 44.72]*1e-0;
%线性-马尔可夫切换-状态空间模型参数
para=[0.07681 0.00753 0.00583 0.19042 0.81964 0.95986 0.86559;%all
    0.13808 0.00013 0       0.81443 2.46036 0.98570 0.77807;%cons
    0.07150  0.00342  0.00176 0.32287  0.88134 0.94506 0.771 ;%finacials
    0.07864 0.00719 0.00455 0.26603 0.82205 0.95576 0.87924;%industry
    0.08631 0.01004 0.00274 0.2148 0.90433 0.96425 0.86470;%small
    0.059 0.00957 0.006 0.20281 0.60805 0.95546 0.87943;%big
    0.07878 0.00561 0.00644 0.17066 0.755 0.96119 0.89147;%concentional
    0.13295 0 0.00489 0.34593 1.137 0.960 0.93088];%islamic
parts={'All Stocks','Consumer Cyclical','Finacials','Industry','Small','Big','Concentional','Islamic'};
%%
%模拟产生数据
k=1;
if 1
    pk=para(k,:)';
else
    pk=[0.0027
        0.0100
        0.0023
        0.0256
        0.2496
        0.9999
        0.9981];
end
sigv=pk(1:3);sige=pk(4:5);p=pk(6:7);p=[0.99,0.99]';
T=365*2;%模拟采样个数
[ss,zz,rr,yy,tt]=simulate(p,sigv,sige,rmin(k),rmax(k),meanr(k),stdr(k),T);
%xlswrite('simulatedata.xlsx',[rr,zz],1);
li=[rr,zz];save('simulatedata.txt','li','-ascii');%保存数据
figure(10);
a=2;b=2;
subplot(a,b,1);plot(tt,zz),title(strcat('CASD of--',' ',parts{k}));
names={'\alpha_1','\alpha_2','herding coefficient (\alpha_3)'};
for i=1:3;subplot(a,b,1+i);plot(tt,yy(:,i));title(names{i});end;

function [W,V,H,F]=getmodel(sigv,sige,r,St)
W=diag(sigv.^2);
V=sige(St)^2;
H=[1,abs(r),r^2];
F=eye(3);
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


    