function lzm=valorMedio(nder,ndt,htt,lz,t)

lzm=zeros(nder,1); 

for i=1:nder
    lzm(i)=mean(lz(i,:));
%     for k=1:ndt
%     
%         g(k)=lz(i,k);
%     
%     end
%     figure(7)
%     plot(t,g),grid on, pause(.1),pause
    
    pp=0; qq=0;

    for k=round(ndt/10):ndt-1

        if lz(i,k)>lz(i,k+1) && pp==0 && qq==0                             % Primer punto del ciclo. Punto A. Punto de máximo.

            be=k-1; v_(1)=lz(i,k); pp=1; max_=lz(i,k); a=t(k);

        elseif lz(i,k)>lz(i,k+1) && pp==1 && qq==0                         % Recopila los puntos de la zona B.

            v_(k-be)=lz(i,k); 

        elseif lz(i,k)<lz(i,k+1) && pp==1                                  % Recopila los puntos de la zona C.

            v_(k-be)=lz(i,k); qq=1; 

        elseif lz(i,k)>lz(i,k+1) && pp==1 && qq==1                         % Se ha vuelto a encontrar con un máximo de la señal.

            if abs(max(v_)-max_)<1e-5                                      % Se compara con el máximo anterior para observar estabilidad de la señal.

                v_(k-be)=lz(i,k); nd=k-be; t_=htt*(nd-1); b=t(k);                  % Último punto del ciclo adecuado para el cálculo del torque promedio lym.
                lzm1(i)=trapecio1D_puntual(nd,t_,htt,v_);  [i,lzm1(i), t(k),mean(lz(i,:)),1];clear v_; %v_, [a b lzm1(i)], pause %ja=1;                % Función que se encarga de calcular el torque promedio para cada r(i).
                break

            else
                max_=max(v_); be=k-1; clear v_; v_(1)=lz(i,k); pp=1; qq=0; a=t(k); %fr=1,pause % Si no se ha estabilizado la señal, se vuelve a realizar el proceso.

            end

        elseif abs(lz(i,k)-lz(i,k+1))<1e-18
            
            lzm1(i)=lz(i,k);[i,lzm1(i), t(k),2];   % jas=k,pause
            break

        end

    end
end