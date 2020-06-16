function [v_an, w_an]=cilindrico_analitica_bessel(nder,eta,eta0,phi,tau,kappa,om) 

omTau=2*pi*om*tau;             % Frecuencia adimensional
zita=1.5*phi*eta0;             % Vortex viscosity
eta_e=eta+zita;                % Constante = eta + zita


c1=-2*eta_e*kappa*besseli(0,kappa)/(4*(1+omTau^2)*(2*zita*besseli(1,kappa)-eta_e*kappa*besseli(0,kappa)))-eta_e/(2*eta*(1+omTau^2));

c2=eta_e*kappa/(4*(1+omTau^2)*(2*zita*besseli(1,kappa)-eta_e*kappa*besseli(0,kappa)));


her=1/(nder-1);                % Paso en coordenada radial   

r=zeros(1,nder);

for k=1:nder
    r(k)=her*(k-1);            % Se establece el dominio de los puntos del espacio en la coordenada radial r=0 a r=1.
end

w_an=zeros(1,nder); v_an=zeros(1,nder);

for i=1:nder
    
    w_an(i)=c1/2+c2*besseli(0,kappa*r(i))+eta_e/(4*eta*(1+omTau^2));
    
    v_an(i)=r(i)/2*c1+2*zita/(eta_e*kappa)*c2*besseli(1,kappa*r(i))+zita*r(i)/(4*eta*(1+omTau^2));
    
end




