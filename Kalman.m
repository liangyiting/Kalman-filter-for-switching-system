
function [yt,Pt,ept,y1,P1]=Kalman(r,yt0,Pt0,St,CSADt,sigv,sige)
[W,V,H,F]=getmodel(sigv,sige,r,St);
y1=F*yt0;
P1=F*Pt0*F'+W;
z1=H*y1;
ep1=CSADt-z1;
yt=y1+P1*H'*((H*P1*H'+V)\ep1);
Pt=(eye(3)-P1*H'*((H*P1*H'+V)\H))*P1;
ept=ep1;


