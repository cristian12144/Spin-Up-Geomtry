
function SpinUpFlowDescriptor(nder,ndet,ndt,tf,maxIterT, varargin) 

tic
close all
clc
datetime('now')

    %% DEFINICIÓN DE LAS VARIABLES DEL SISTEMA
    
if nargin > 5 
    BmT = varargin{1};
else
    BmT=[3.4 4.5 5.6 6.8 7.9 9.0 10.1 11.3]; % Valores de densidad de campo magnetico en mT.
end

if nargin > 6 
    ff = varargin{2};
else
    ff={}; % Valores de densidad de campo magnetico en mT.
end

if ~isfield(ff,'name'),  ff.name = 'WBF1'; end                % Ferrofluid name
if ~isfield(ff,'eta'),   ff.eta=1.03e-3;   end                % Sheear viscosity
if ~isfield(ff,'eta0'),  ff.eta0=1.02e-3;  end                % Viscosidad del liquido portador de las nanopartículas
if ~isfield(ff,'ji'),    ff.ji=0.106;      end                % susceptibilidad magnética inicial
if ~isfield(ff,'phi'),   ff.phi=2.13e-3;   end                % Fracción volumétrica
if ~isfield(ff,'tau'),   ff.tau=1.67e-5;   end                % Tiempo de relajación (Browniano) [s]
if ~isfield(ff,'kappa'), ff.kappa=3.3;     end                % Constante
if ~isfield(ff,'om'),    ff.om=150;        end                % Frecuencia de rotación del campo magnético [Hz]
if ~isfield(ff,'R0'),    ff.R0=24.7e-3;    end                % Radio del cilindro REPORTADO POR TORRES-DIAZ.
if ~isfield(ff,'Md'),    ff.Md=425e3;      end                % Magnetización del dominio magnético de las nanopartículas
if ~isfield(ff,'T'),     ff.T=294;         end                % Temperatura absoluta del sistema
if ~isfield(ff,'d'),     ff.d=14.3e-9;     end                % Diámetro promedio de nanopartículas magnéticas

% filenm = sprintf( 'simulationData_%s.mat',ff.name);   
% save( filenm , 'ff' );

u0=4*pi*10^-7;                       % Permeabilidad del aire o vacío
omTau=2*pi*ff.om*ff.tau;             % Frecuencia adimensional
zita=1.5*ff.phi*ff.eta0;             % Vortex viscosity
eta_e=ff.eta+zita;                   % Constante = eta + zita

Kb=1.38064852e-23;

    

her=1/(nder-1);                % Paso en coordenada radial
het=2*pi/ndet;                 % Paso en coordenada angular
htt=tf/(ndt-1);                % Paso en coordenada temporal

lz0=omTau/(1+omTau^2);         % Torque análitico de orden cero.

r=zeros(nder,1);  teta=zeros(ndet,1); t=zeros(ndt,1);  % Inicialización de variables
HrP=zeros(nder,ndt);HtP=zeros(nder,ndt);MrP=zeros(nder,ndt);MtP=zeros(nder,ndt);lzP=zeros(nder,ndt); cont=zeros(nder,1);
g1=zeros(ndt,1);g2=zeros(ndt,1);g3=zeros(ndt,1);g4=zeros(ndt,1);g5=zeros(ndt,1);


for k=1:nder
    r(k)=her*(k-1);%+(R1/R0); % Se establece el dominio de los puntos del espacio en la coordenada radial r=0 a r=1.
end

for k=1:ndet
    teta(k)=het*(k-1); % Se establece el dominio de los puntos del espacio en la coordenada tetha tetha=0 a tetha=2pi.
end

for k=1:ndt
    t(k)=htt*(k-1); % Se establece el dominio de los puntos del espacio en la coordenada tiempo t=0 a t=tf.
end

Mr=zeros(nder,ndet,ndt); Mt=zeros(nder,ndet,ndt); lzm=zeros(nder,1); mr_teory=zeros(ndt,1); mt_teory=zeros(ndt,1); lz=zeros(nder,ndet,ndt);


B=eta_e/(ff.eta*her^2);           % Constante de la discretización de las ecuaciones diferenciales hidrodinámicas
Cx=eta_e/(ff.eta*2*her);          % Constante de la discretización de las ecuaciones diferenciales hidrodinámicas
D=-eta_e/ff.eta;                  % Constante de la discretización de las ecuaciones diferenciales hidrodinámicas
E=-zita/(her*ff.eta);             % Constante de la discretización de las ecuaciones diferenciales hidrodinámicas
F=4*ff.eta/(eta_e*ff.kappa^2);       % Constante de la discretización de las ecuaciones diferenciales hidrodinámicas


%% FUNCIONES PARA LA VALIDACIÓN DE LOS RESULTADOS NUMÉRICOS
[v_an, w_an]=cilindrico_analitica_bessel(nder,ff.eta,ff.eta0,ff.phi,ff.tau,ff.kappa,ff.om); % SOLUCIÓN ANALÍTICA, A TRAVÉS DE FUNCIONES DE BESSEL MODIFICADAS
[v1,w1]=cilindrico_orden0_v_w_if(nder,her,r,B,Cx,D,E,F,lz0,omTau);        % HALLA PERFILES CON EL TORQUE TEÓRICO lz0- EN FUNCIÓN DE R.

%% VALOR INICIAL DEL TORQUE PROMEDIO

lzm=zeros(nder,1);

thetaGrafica=round(ndet/2); radioGrafica=round(nder/2);                     % AQUÍ SE ESTABLECE QUE PARTE SE QUIERE MOSTRAR EN LAS GRÁFICAS (SÓLO ES UNA DE LAS 2).



alf=zeros(1,length(BmT));

vf=zeros(nder,length(BmT));wf=zeros(nder,length(BmT)); 






for campo=1:length(BmT)
    
    K=BmT(campo)*1e-3/u0;                                                   % VALOR DE CAMPO MAGNETICO DE LA ITETACIÓN 'campo'
    
    alf(campo)=pi/6*(u0*ff.Md*ff.d^3*K)/(Kb*ff.T);
    epsi=u0*ff.ji*K^2*ff.tau/zita;
    
    for xx=1:1

        fprintf('RESULTADOS OBTENIDOS PARA WBF1: \n\n');
        fprintf('kappa     = %.2f   [dimensionless]\n',ff.kappa);
        fprintf('f         = %.0f    [Hz] \n',ff.om);
        fprintf('omTau     = %.3f  [dimensionless] \n',omTau);
        fprintf('B         = %.1f    [mT] \n',K*u0/1e-3);
        fprintf('her       = %.4f [dimensionless] \n',her);
        fprintf('het       = %.4f [dimensionless] \n',het);
        fprintf('htt       = %.4f [dimensionless] \n',htt);
        fprintf('alpha     = %.2f   [dimensionless] \n',alf(campo));
        fprintf('epsi      = %.3f  [dimensionless] \n',epsi);%,pause
        fprintf('lz0       = %.4f [dimensionless] \n\n\n',lz0);%pause
        
    end                                                           % DISPLAY DE LOS VALORES CON LOS QUE SE ESTÁ HACIENDO LA SIMULACIÓN
    

    for m=1:maxIterT
        
        [v,w]=cilindrico_orden0_v_w_if_numerico2(nder,her,r,B,Cx,D,E,F,lzm,omTau,thetaGrafica);
        
        
        %% CAMPO MAGNETICO
           
        [Hr, Ht]=CampoMagnetico_CamposAltos(nder,ndet,ndt,tf,Mr,Mt,ff.ji);     % CÁLCULO DEL CAMPO MAGNÉTICO

        %% CÁLCULOS DE LA MAGNETIZACIÓN Y TORQUE MAGNÉTICO INSTANTÁNEO
        
        [Mr, Mt, lz]=NR_(nder,ndet,ndt,r,Mr,Mt,v,w,epsi,omTau,Hr,Ht,htt,het);
        
        %%
    


            for i=1:nder
                for k=1:ndt
                    HrP(i,k)=Hr(i,thetaGrafica,k);
                    HtP(i,k)=Ht(i,thetaGrafica,k);
                    MrP(i,k)=Mr(i,thetaGrafica,k);
                    MtP(i,k)=Mt(i,thetaGrafica,k);
                    lzP(i,k)=lz(i,thetaGrafica,k);
                end
            end                                                     % PREPARACIÓN DE LOS VECOTRES PARA EL PLOTEO DEL CAMPO MAGNÉTICO Y LA MAGNETIZACIÓN Y VECTOR PARA EL CÁLCULO DEL TORQUE PROMEDIO DE LA ITERACIÓN m 


            for k=1:ndt

                g1(k)=MrP(radioGrafica,k);
                g2(k)=MtP(radioGrafica,k);
                g3(k)=HrP(radioGrafica,k);
                g4(k)=HtP(radioGrafica,k);
                g5(k)=lzP(radioGrafica,k);

            end                                                      

                                 
        
        
        %% CÁLCULO DEL TORQUE PROMEDIO
        
        lzm_=lzm;                                                                                % AQUÍ SE GUARDA EL VALOR DEL TORQUE PROMEDIO CALCULADO EN LA ITERACIÓN ANTERIOR.
        
        lzm=valorMedio(nder,ndt,htt,lzP,t);                                                      % CÁLCULO DEL TORQUE PROMEDIO DE ITERACIÓN ACTUAL 
         
        
                                                                                                 %% NUEVO CÁLCULO DE LOS PERFILES CON NUEVO TORQUE PROMEDIO
       
        [v,w]=cilindrico_orden0_v_w_if_numerico2(nder,her,r,B,Cx,D,E,F,lzm,omTau,thetaGrafica);  % HALLA NUEVOS PERFILES CON NUEVO TORQUE NUMÉRICO OBTENIDO
        
                                                                                                 %% EVALUACIÓN DE LA CONVERGENCIA EN EL TORQUE PROMEDIO

            for i=round(nder*.30):round(nder)

                cont(i)=abs((lzm(i)-lzm_(i))/lzm(i))*100;                                           

            end                                                 % DIFERENCIA DEL TORQUE PROMEDIO ACTUAL CON EL ANTERIOR


            if max(cont)<3.0                                                                      % CRITERIO DE CONVERGENCIA DEL TORQUE PROMEDIO 3% (T_actual - T_anterior).

                break
                
            end

     end
    
%         fprintf('Iteraciones requeridas         = %i \n\n',m);                                  % REPORTA EL NÚMERO DE ITERACIONES REALIZADAS POR EL PROGRAMA.
    
    
                                                                                                %% SACA LOS PERFILES CON DIMENSIONES [mm/s] Y [rad/s]
    
    for xx=1:1
    
        vd=      v.*u0*ff.ji*K^2*omTau*ff.R0/zita*1000;                                               % radio de 1 pulgada=0.0254 metros  --  se multiplica por 1000 para dar en [mm/s]
        v1d=    v1.*u0*ff.ji*K^2*omTau*ff.R0/zita*1000;
        v_and=v_an.*u0*ff.ji*K^2*omTau*ff.R0/zita*1000;
        
        wd=      w.*u0*ff.ji*K^2*omTau/(zita);                                                     % se divide en 2*pi para que dé en Hertz.
        w1d=    w1.*u0*ff.ji*K^2*omTau/(zita);
        w_and=w_an.*u0*ff.ji*K^2*omTau/(zita);
        
        rd=r.*ff.R0*1000; 
        
    end                                                                               % VALORES DE LOS PERFILES CON DIMENSIONES [mm/s] y [rad/s]
    
    vf(:,campo)=vd; wf(:,campo)=wd;                                                   % PERFIL PARA EL VALOR DE CAMPO MAGNÉTICO DEL VECTOR BmT(campo)

  
  
    %%
    
    hoja4=figure(4);                                                               % GRAFICA LOS PERFILES DE VELOCIDAD LINEAL PARA CADA VALOR DE CAMPO MAGNÉTICO
    subplot(2,1,1)
    h=plot(rd,vf(:,campo),'o:','DisplayName',[num2str(BmT(campo)),' mT']); set(h,'LineWidth',1.5,'MarkerSize',2)%'Color',colorVec(campo,:));
    hold on, grid on, title(['WBF1 -- ','\kappa = ',num2str(ff.kappa),' --  f = ',num2str(ff.om),' Hz']);% -- \nabla (xM+H) -- \Omega\epsilonv\nablaM'
    xlabel('r [mm]','FontWeight','bold');ylabel('v [mm/s]','FontWeight','bold');xlim([0 ff.R0*1000]);
    legend('Location','bestoutside');

    subplot(2,1,2)                                                          % GRAFICA LOS PERFILES DE VELOCIDAD ANGULAR PARA CADA VALOR DE CAMPO MAGNÉTICO
    h=plot(rd,wf(:,campo),'o:','DisplayName',[num2str(BmT(campo)),' mT']); 
    hold on, grid on, set(h,'LineWidth',1.5,'MarkerSize',2);
    xlabel('r [mm]','FontWeight','bold');ylabel('\omega [rad/s]','FontWeight','bold');xlim([0 ff.R0*1000]);
    legend('Location','bestoutside');
    
    
    set(hoja4,'Units','Inches');
    pos=get(hoja4,'position');
    set(hoja4,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)])
    print(hoja4,'cilindro diferentes frecuencias Torres-Diaz','-dpdf','-r0')

%%
    lzm=zeros(nder,1);                                                      % REINICIA EL VALOR DEL TORQUE PROMEDIO A CERO PARA PROCEDER CON EL SIGUIENTE VALOR DE CAMPO MAGNÉTICO 
    
    
end                                                                         % FINALIZA LA ITERACIÓN PARA UN DETERMINADO VALOR DE CAMPO MAGNÉTICO Y COMIENZA CON OTRO.

 

% for i=1:nder
%     r(i)=r(i)*R0*1000;
% end
% 
% vMRShVal=vf(:,1);      wMRShVal=wf(:,1);
% vMRShConv3=vf(:,2);  wMRShConv3=wf(:,2);
% vMRShConv4=vf(:,3);  wMRShConv4=wf(:,3);
% vMRShConv5=vf(:,4);  wMRShConv5=wf(:,4);
% vMRShConv6=vf(:,5);  wMRShConv6=wf(:,5);
% vMRShConv7=vf(:,6);  wMRShConv7=wf(:,6);
% vMRShConv9=vf(:,7);  wMRShConv9=wf(:,7);
% vMRShConv10=vf(:,8); wMRShConv10=wf(:,8);
% vMRShConv11=vf(:,9); wMRShConv11=wf(:,9);
% 
% 
% %%
% 
% 
% hoja30=figure(30);                                                               % PERFILES DE VELOCIDAD LINEAL ANALÍTICO, NUMÉRICO SÓLO HIDRODINÁMICO, NUMÉRICO TOTAL. ESTO LO HAGO SÓLO PARA LA VALIDACIÓN A CAMPOS BAJOS.
% subplot(2,1,1)
% h=plot(rd,vMRShVal,':ok',rd,v_and,'k',rd,v1d,':pk'); set(h(1),'linewidth',2,'MarkerSize',2);set(h(2),'linewidth',2);set(h(3),'linewidth',2,'MarkerSize',2); 
% hold on, grid on, %title(['WBF1 -- B = ',num2str(BmT(campo)),' mT -- ','\alpha = ',num2str(round(alf,2)),'-- \kappa = ',num2str(kappa),' --  f = ',num2str(om),' Hz -- MRSh -- \nabla (xM+H) -- \Omega\epsilonv\nablaM']);
% xlabel('r [mm]','FontWeight','bold');ylabel('v [mm/s]','FontWeight','bold');
% legend('v_\theta MRSh','v_\theta analitic','v_\theta lz0','Location','south'), xlim([0 R0*1000]);
% 
% subplot(2,1,2)                                                          % PERFILES DE VELOCIDAD ANGULAR ANALÍTICO, NUMÉRICO SÓLO HIDRODINÁMICO, NUMÉRICO TOTAL.
% h=plot(rd,wMRShVal,':ok',rd,w_and,'k',rd,w1d,':pk'); set(h(1),'linewidth',2,'MarkerSize',2);set(h(2),'linewidth',2);set(h(3),'linewidth',2,'MarkerSize',2);
% hold on, grid on, set(h(1),'linewidth',2);set(h(2),'linewidth',2);set(h(3),'linewidth',2);%  title('Spin');
% xlabel('r [mm]','FontWeight','bold');ylabel('\omega [rad/s]','FontWeight','bold'); xlim([0 R0*1000]);
% legend('\omega_z MRSh','\omega_z analitic','\omega_z lz0');pause(.1)
% 
% set(hoja30,'Units','Inches');
% pos=get(hoja30,'position');
% set(hoja30,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)])
% print(hoja30,'cilindroValidacionNegro','-dpdf','-r0')
%     
%     
% 
% hoja50=figure(50);  
% figure(50)
% subplot(2,1,1)
% h=plot(rd,vMRShConv3,':ok',rd,vMRShConv4,':+k',rd,vMRShConv5,':*k',rd,vMRShConv6,':xk',rd,vMRShConv7,':sk',rd,vMRShConv9,':dk',rd,vMRShConv10,':pk',rd,vMRShConv11,':hk');  hold on, grid on, %title(['WBF1 -- ','\kappa = ',num2str(kappa),' --  f = ',num2str(om),' Hz -- MRSh -- \nabla (xM+H) -- \Omega\epsilonv\nablaM']);
% for a=1:length(BmT)
%     set(h(a),'LineWidth',1.5,'MarkerSize',2.5);
% end
% xlabel('r [mm]','FontWeight','bold');ylabel('v [mm/s]','FontWeight','bold');xlim([0 R0*1000]);
% legend('3.4 mT','4.5 mT','5.6 mT','6.8 mT','7.9 mT','9.0 mT','10.1 mT','11.3 mT','Location','bestoutside'),
% 
% subplot(2,1,2)
% h=plot(rd,wMRShConv3,':ok',rd,wMRShConv4,':+k',rd,wMRShConv5,':*k',rd,wMRShConv6,':xk',rd,wMRShConv7,':sk',rd,wMRShConv9,':dk',rd,wMRShConv10,':pk',rd,wMRShConv11,':hk'); hold on, grid on,% title('Spin');
% for a=1:length(BmT)
%     set(h(a),'LineWidth',1.5,'MarkerSize',2.5);
% end
% xlabel('r [mm]','FontWeight','bold');ylabel('\omega [rad/s]','FontWeight','bold');xlim([0 R0*1000]);
% legend('3.4 mT','4.5 mT','5.6 mT','6.8 mT','7.9 mT','9.0 mT','10.1 mT','11.3 mT','Location','bestoutside'),
% 
% 
% set(hoja50,'Units','Inches');
% pos=get(hoja50,'position');
% set(hoja50,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)])
% print(hoja50,'CilindroDiferentesCamposNegro','-dpdf','-r0')



% hoja51=figure(51);  
% figure(51)
% h=plot(rd,vMRShConv3,':^c',rd,vMRShConv4,':>m',rd,vMRShConv5,':<b',rd,vMRShConv6,':pg',rd,vMRShConv7,':h',rd,vMRShConv9,':dr',rd,vMRShConv10,':py',rd,vMRShConv11,':h');  hold on, grid on, %title(['WBF1 -- ','\kappa = ',num2str(kappa),' --  f = ',num2str(om),' Hz -- MRSh -- \nabla (xM+H) -- \Omega\epsilonv\nablaM']);
% for a=1:length(BmT)
%     set(h(a),'LineWidth',2,'MarkerSize',4);
% end
% xlabel('r [mm]','FontWeight','bold');ylabel('v [mm/s]','FontWeight','bold');xlim([0 R0*1000]);ylim([0 9]);
% legend('3.4 mT','4.5 mT','5.6 mT','6.8 mT','7.9 mT','9.0 mT','10.1 mT','11.3 mT','Location','best'),
% 
% set(hoja51,'Units','Inches');
% pos=get(hoja51,'position');
% set(hoja51,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)])
% print(hoja51,'CilindroDiferentesCamposV','-dpdf','-r0')
% 
% 
% hoja52=figure(52);  
% figure(52)
% h=plot(rd,wMRShConv3,':^c',rd,wMRShConv4,':>m',rd,wMRShConv5,':<b',rd,wMRShConv6,':pg',rd,wMRShConv7,':h',rd,wMRShConv9,':dr',rd,wMRShConv10,':py',rd,wMRShConv11,':h');  hold on, grid on, %title(['WBF1 -- ','\kappa = ',num2str(kappa),' --  f = ',num2str(om),' Hz -- MRSh -- \nabla (xM+H) -- \Omega\epsilonv\nablaM']);
% for a=1:length(BmT)
%     set(h(a),'LineWidth',2,'MarkerSize',4);
% end
% xlabel('r [mm]','FontWeight','bold');ylabel('\omega [rad/s]','FontWeight','bold');xlim([0 R0*1000]);
% legend('3.4 mT','4.5 mT','5.6 mT','6.8 mT','7.9 mT','9.0 mT','10.1 mT','11.3 mT','Location','best'),
% 
% set(hoja52,'Units','Inches');
% pos=get(hoja52,'position');
% set(hoja52,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)])
% print(hoja52,'CilindroDiferentesCamposW','-dpdf','-r0')
% 

toc
end

