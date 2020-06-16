function [Hr, Ht]=CampoMagnetico_CamposAltos(nder,ndet,ndt,tf,Mr,Mt,ji)

her=1/(nder-1);                % Paso en coordenada radial
het=2*pi/ndet;                 % Paso en coordenada angular
htt=tf/(ndt-1);                % Paso en coordenada temporal

r=zeros(nder,1);  teta=zeros(ndet,1); t=zeros(ndt,1);  % Inicialización de variables

for k=1:nder
    r(k)=her*(k-1); % Se establece el dominio de los puntos del espacio en la coordenada radial r=0 a r=1.
end

for k=1:ndet
    teta(k)=het*(k-1); % Se establece el dominio de los puntos del espacio en la coordenada tetha tetha=0 a tetha=2pi.
end

for k=1:ndt
    t(k)=htt*(k-1); % Se establece el dominio de los puntos del espacio en la coordenada tiempo t=0 a t=tf.
end


Hr=zeros(nder,ndet,ndt); Ht=zeros(nder,ndet,ndt); %Hkn=zeros(nder,ndet,ndt);
Hrk=zeros(nder*ndet,ndt); Htk=zeros(nder*ndet,ndt); Hkk=zeros(nder*ndet,ndt);

a=zeros(nder,1); b=zeros(nder,1); Z=zeros(nder,1); Zp=zeros(nder,1);  W=zeros(nder,1);

for i=1:nder
    
    a(i)=1+r(i)/her;
    b(i)=-r(i)/her;
    Z(i)=r(i)/her;
    Zp(i)=1+Z(i);
    
    W(i)=ji*r(i)/her;
    
end

c=1/het; cp=ji/het;

%%

A=zeros(2*nder*ndet,2*nder*ndet); C=zeros(2*nder*ndet,ndt);



for i=1:nder-1
    for j=1:ndet
        
        A((i-1)*ndet+j,i*ndet+j)=-c;
        
        if j-1==0
            
            A((i-1)*ndet+j,(i+1)*ndet+j-1)=c;
            
        else
            
            A((i-1)*ndet+j,i*ndet+j-1)=c; %OK
            
        end
        
        A((i-1)*ndet+j,nder*ndet+(i-1)*ndet+j)=-Z(i+1);
        A((i-1)*ndet+j,nder*ndet+(i-0)*ndet+j)=Zp(i+1); %OK
        
        
        
        
        
        
        A(nder*ndet+(i-1)*ndet+j,(i-1)*ndet+j)=-Z(i+1); % ACÁ COMETÍIII UN ERROR QUE ME COSTÓ 2 DÍAS DE RETRASO. TAN SOLO OLVIDÉ COLOCAR (i-1) OJO !!!
        A(nder*ndet+(i-1)*ndet+j,(i-0)*ndet+j)=Zp(i+1); %OK
        A(nder*ndet+(i-1)*ndet+j,nder*ndet+(i-0)*ndet+j)=c;
        
        if j-1==0
            
            A(nder*ndet+(i-1)*ndet+j,nder*ndet+(i+1)*ndet+j-1)=-c;
            
            
        else
            
            A(nder*ndet+(i-1)*ndet+j,nder*ndet+(i+0)*ndet+j-1)=-c; %OK
            
            
        end
        
        
        
        
    end
    
end


for i=1:ndet
    


    A(ndet*(nder-1)+i,i)=              1; % CONDICION DE FRONTERA EN EL CENTRO DEL CILINDRO EN DONDE LA COMPONENTE NORMAL DEL FLUX B DEBE SER CONTINUO Y LA UNICA FORMA ES QUE SEA CERO!
    A(2*nder*ndet+1-i,2*nder*ndet+1-i)=1;
    
end

A=sparse(A); %OK
% A, pause
%% EN ESTA PARTE SE HACEN LOS VALORES DE C QUE DEPENDEN DEL TIEMPO.

for k=1:ndt
    for i=1:nder-1
        for j=1:ndet
    
            if j==1


                C(nder*ndet+(i-1)*ndet+j,k)=-ji*Zp(i+1)*Mr(i+1,j,k)+W(i+1)*Mr(i,j,k)-cp*Mt(i+1,j,k)+cp*Mt(i+1,ndet,k);

            else


                C(nder*ndet+(i-1)*ndet+j,k)=-ji*Zp(i+1)*Mr(i+1,j,k)+W(i+1)*Mr(i,j,k)-cp*Mt(i+1,j,k)+cp*Mt(i+1,j-1,k);

            end
    
        end
    end

    
    for j=1:ndet
        
        C(ndet*(nder-1)+j,k)=-ji*Mr(1,j,k);              % CONDICIÓN DE FRONTERA EN R=0!
%       C(ndet*(nder-1)+j,k)=-(ur-1)*Mr(1,j,k);                                                                                                                                                      %         C(1*nder*ndet+1-j,k)=+sin(t(k)-teta(ndet+1-j));  % CONDICIÓN DE FRONTERA SENGÚN LO ESTABLECIDO EN LA TESIS DE RINALDI P-129 -ji*Mr(nder,ndet+1-j,k);%
        C(2*nder*ndet+1-j,k)=-cos(t(k)-teta(ndet+1-j));  % CONDICIÓN DE FRONTERA APLICANDO LA ECUACIÓN PARA QUE EXISTA CAMPO ROTATIVO
        
    end
    
%     det(A),pause
    res=A\C(:,k);
    
    
    
    for i=1:nder*ndet
        
        Hrk(i,k)=res(i);
        Htk(i,k)=res(nder*ndet+i);
        Hkk(i,k)=sqrt(Hrk(i,k)^2+Htk(i,k)^2);
        
    end
    
    
    Hr_=reshape(Hrk(:,k),ndet,nder); % DEBIDO A LA FORMA EN QUE MATLAB REORDENA EL VECTOR Hrk DE UN VECTOR COLUMNA A UNA MATRIZ (ES POR ESO, OJO AHÍ)
    Ht_=reshape(Htk(:,k),ndet,nder); % POR ESO ES QUE EN ESTAS TRES LINEAS ESTÁ ndet,nder Y NO AL CONTRARIO.
%     Hk_=reshape(Hkk(:,k),ndet,nder);
    
    
   
    Hr(:,:,k)=Hr_';                 % Para que quede como debe cada pareja (i,j) sea "i" el radio y "j" el ángulo.
    Ht(:,:,k)=Ht_';                 % Para que quede como debe cada pareja (i,j) sea "i" el radio y "j" el ángulo.
%     Hkn(:,:,k)=Hk_';                 % Para que quede como debe cada pareja (i,j) sea "i" el radio y "j" el ángulo.
    
end                               % Aquí finaliza la el cálculo de las componentes de los campos H para un tiempo t(k).

% Hr=Hrn;
% Ht=Htn;
    
  