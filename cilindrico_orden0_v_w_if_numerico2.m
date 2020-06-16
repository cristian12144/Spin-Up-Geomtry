function [v,w]=cilindrico_orden0_v_w_if_numerico2(nder,her,r,B,C,D,E,F,lzm,omTau,thetaGrafica)

for i=1:nder
    
    A=zeros(2*nder,2*nder);
    CC=zeros(2*nder,1);
    A(1,1)=1;                                                       % CONDICIÓN DE FRONTERA EN R=0, v (0) = 0.
%     A(1,1)=-3; A(1,2)=4; A(1,3)=-1;                               % CONDICIÓN DE FRONTERA EN R=0, v'(0) = 0.
    A(nder+1,nder+1)=-3; A(nder+1,nder+2)=4; A(nder+1,nder+3)=-1;   % CONDICIÓN DE FRONTERA EN R=0, w'(0) = 0.
    
    
    A(nder,nder)=1;                                                 % CONDICIÓN DE FRONTERA EN R=1,  v(1) = 0.
    A(2*nder,2*nder)=1;                                             % CONDICIÓN DE FRONTERA EN R=1,  w(1) = 0.
    
end

for i=2:nder-1
    
    A(i,i-1)=B-C/r(i);
    A(i,i)=D/r(i)^2-2*B;
    A(i,i+1)=B+C/r(i);
    A(i,nder+i-1)=-E;
    A(i,nder+i+1)=E;
    
    A(nder+i,i-1)=-1/her;
    A(nder+i,i)=2/r(i);
    A(nder+i,i+1)=1/her;
    A(nder+i,nder+i-1)=F/her^2-F/(2*her*r(i));
    A(nder+i,nder+i)=-(4+2*F/her^2);
    A(nder+i,nder+i+1)=F/(2*her*r(i))+F/her^2;
    
    if i<=round(nder*30/100)
          CC(nder+i)=-lzm(round(nder*50/100)+1)/omTau; %CC(nder+i)=-r(i)*abs(lzm(round(nder*30/100))/omTau);  % HAGO ESTA CONDICIÓN, YA QUE LOS PUNTOS CERCANOS AL CENTRO i<=3 PRODUCEN ERROR EN LOS RESULTADOS (NO SÉ XQ)
    else
        CC(nder+i)=-lzm(i)/omTau;
    end
    
end

A=sparse(A);

res=A\CC;

v=zeros(nder,1); w=zeros(nder,1);

for i=1:nder
    
    v(i)=res(i);
    w(i)=res(nder+i);
    
end


