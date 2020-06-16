function G=Grc(r,Mr,Mt,Mtk_,Mtj_,v,w,epsi,omTau,Hr,Ht,htt,het)
             
norma=sqrt(Hr^2+Ht^2);                          % Norma del vector H

alpha=sqrt(3/2*epsi*norma^2);                   % Funcion Alfa

if alpha<3.189615789606358e-08                  % Esto se hace debido a la limitación computacional al calcular << coth >>
    
    alpha=3.189615789606358e-08;
    
end

L=coth(alpha)-1/alpha;                         % Función Langevin << L >> L(alfa)

Phi=(coth(alpha))/alpha-1/alpha^2;             % Función Phi = L/alpha

Bper=2*L/(alpha-L);                            % Parámetro de la ecuación MRSh (Componente perpendicular de B)

delta=1e-3;

Ldelta=coth(alpha+delta)-1/(alpha+delta);      % L(alfa+delta)

dlnL=log(Ldelta)-log(L);                       % ln[L(alfa+delta)]-ln[L(alfa)]

Bpar=alpha/delta*dlnL;                         % Parámetro de la ecuación MRSh (Componente perpendicular de B) = Alfa/delta*{ln[L(Alfa+delta)]-ln[L(Alfa)]}


G=omTau*(Mt-Mtk_)/htt+omTau*epsi*v*(1/(r*het)*(Mt-Mtj_)+Mr/r)-omTau*epsi*w*Mr-Ht*(Hr*Mr+Ht*Mt)/norma^2*(1/Bper-1/Bpar)-3*Ht*Phi/Bpar+Mt/Bper;       % Componente HORIZONTAL de MRSh.



% G=omTau*(Mt-Mtk_)/htt+omTau*epsi*v*(1/(r*het)*(Mt-Mtj_)+Mr/r)-omTau*epsi*w*Mr-Ht*(Hr*Mr+Ht*Mt)/norma^2*(1/Bper-1/Bpar)-3*Ht*Phi/Bpar+Mt/Bper;       % Componente HORIZONTAL de MRSh.

