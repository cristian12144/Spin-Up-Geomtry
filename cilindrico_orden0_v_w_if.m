function [v,w]=cilindrico_orden0_v_w_if(nde,h,radio,B,C,D,E,F,lz0,omTau)


A=zeros(2*nde,2*nde);
CC=zeros(2*nde,1);
A(1,1)=1;
% A(1,1)=-3; A(1,2)=4; A(1,3)=-1;                               % CONDICIÓN DE FRONTERA EN R=1, v'(0) = 0.
A(nde+1,nde+1)=-3; A(nde+1,nde+2)=4; A(nde+1,nde+3)=-1;


A(nde,nde)=1;                                                 % CONDICIÓN DE FRONTERA EN R=1,  v(1) = 0.
A(2*nde,2*nde)=1;

for i=2:nde-1
    
    A(i,i-1)=B-C/radio(i);
    A(i,i)=D/radio(i)^2-2*B;
    A(i,i+1)=B+C/radio(i);
    A(i,nde+i-1)=-E;
    A(i,nde+i+1)=E;
    
    A(nde+i,i-1)=-1/h;
    A(nde+i,i)=2/radio(i);
    A(nde+i,i+1)=1/h;
    A(nde+i,nde+i-1)=F/h^2-F/(2*h*radio(i));
    A(nde+i,nde+i)=-(4+2*F/h^2);
    A(nde+i,nde+i+1)=F/(2*h*radio(i))+F/h^2;
    CC(nde+i)=-lz0/omTau;
            
end

A=sparse(A);

res=A\CC;
v=zeros(nde,1); w=zeros(nde,1);
for i=1:nde
    
    v(i)=res(i);
    w(i)=res(nde+i);
    
end


