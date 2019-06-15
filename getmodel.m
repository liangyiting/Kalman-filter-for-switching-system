function [W,V,H,F]=getmodel(sigv,sige,r,St)
W=diag(sigv.^2);
V=sige(St)^2;
H=[1,abs(r),r^2];
F=eye(3);