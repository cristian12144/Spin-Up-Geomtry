function [Mr, Mt, lz]=NR_(nder,ndet,ndt,r,Mr,Mt,v,w,epsi,omTau,Hr,Ht,htt,het)

lz=zeros(nder,ndet,ndt);

delta=1e-7;

for k=2:ndt

    for i=2:nder % para evitar la división por cero en r(1) que es igual a 0.

        for j=1:ndet

            r0=0; t0=0; r0d=r0+delta; t0d=t0+delta;

            if j==1

                for l=1:1000

                    f=Frc(r(i),r0,t0,Mr(i,j,k-1),Mr(i,ndet,k),v(i),w(i),epsi,omTau,Hr(i,j,k),Ht(i,j,k),htt,het);
                    g=Grc(r(i),r0,t0,Mt(i,j,k-1),Mt(i,ndet,k),v(i),w(i),epsi,omTau,Hr(i,j,k),Ht(i,j,k),htt,het);

                    frd=Frc(r(i),r0d,t0,Mr(i,j,k-1),Mr(i,ndet,k),v(i),w(i),epsi,omTau,Hr(i,j,k),Ht(i,j,k),htt,het);
                    ftd=Frc(r(i),r0,t0d,Mr(i,j,k-1),Mr(i,ndet,k),v(i),w(i),epsi,omTau,Hr(i,j,k),Ht(i,j,k),htt,het);

                    grd=Grc(r(i),r0d,t0,Mt(i,j,k-1),Mt(i,ndet,k),v(i),w(i),epsi,omTau,Hr(i,j,k),Ht(i,j,k),htt,het);
                    gtd=Grc(r(i),r0,t0d,Mt(i,j,k-1),Mt(i,ndet,k),v(i),w(i),epsi,omTau,Hr(i,j,k),Ht(i,j,k),htt,het);

                    dfr=(frd-f)/delta; dft=(ftd-f)/delta; dgr=(grd-g)/delta; dgt=(gtd-g)/delta;

                    J=[dfr, dft; dgr, dgt];det(J);

                    if isnan(det(J))

                        break

                    end

                     dM=J^-1*[-f;-g]; 

                    if max(abs(dM))<1e-11

                        Mr(i,j,k)=r0; Mt(i,j,k)=t0; lz(i,j,k)=Mr(i,j,k)*Ht(i,j,k)-Mt(i,j,k)*Hr(i,j,k);

                        break

                    end

                    if l==1000

                        Mr(i,j,k)=r0; Mt(i,j,k)=t0; lz(i,j,k)=Mr(i,j,k)*Ht(i,j,k)-Mt(i,j,k)*Hr(i,j,k);

                    end

                    r0=r0+dM(1); t0=t0+dM(2); r0d=r0+delta; t0d=t0+delta;

                end

            else


                for l=1:1000

                    f=Frc(r(i),r0,t0,Mr(i,j,k-1),Mr(i,j-1,k),v(i),w(i),epsi,omTau,Hr(i,j,k),Ht(i,j,k),htt,het);
                    g=Grc(r(i),r0,t0,Mt(i,j,k-1),Mt(i,j-1,k),v(i),w(i),epsi,omTau,Hr(i,j,k),Ht(i,j,k),htt,het);


                    frd=Frc(r(i),r0d,t0,Mr(i,j,k-1),Mr(i,j-1,k),v(i),w(i),epsi,omTau,Hr(i,j,k),Ht(i,j,k),htt,het);
                    ftd=Frc(r(i),r0,t0d,Mr(i,j,k-1),Mr(i,j-1,k),v(i),w(i),epsi,omTau,Hr(i,j,k),Ht(i,j,k),htt,het);

                    grd=Grc(r(i),r0d,t0,Mt(i,j,k-1),Mt(i,j-1,k),v(i),w(i),epsi,omTau,Hr(i,j,k),Ht(i,j,k),htt,het);
                    gtd=Grc(r(i),r0,t0d,Mt(i,j,k-1),Mt(i,j-1,k),v(i),w(i),epsi,omTau,Hr(i,j,k),Ht(i,j,k),htt,het);

                    dfr=(frd-f)/delta; dft=(ftd-f)/delta; dgr=(grd-g)/delta; dgt=(gtd-g)/delta;

                    J=[dfr, dft; dgr, dgt]; det(J); 

                    if isnan(det(J))
                        break


                    end

                    dM=J^-1*[-f;-g];

                    if max(abs(dM))<1e-11

                        Mr(i,j,k)=r0; Mt(i,j,k)=t0; lz(i,j,k)=Mr(i,j,k)*Ht(i,j,k)-Mt(i,j,k)*Hr(i,j,k);

                        break

                    end

                    if l==1000

                       Mr(i,j,k)=r0; Mt(i,j,k)=t0; lz(i,j,k)=Mr(i,j,k)*Ht(i,j,k)-Mt(i,j,k)*Hr(i,j,k);

                    end

                    r0=r0+dM(1); t0=t0+dM(2); r0d=r0+delta; t0d=t0+delta;

                end

            end

        end

    end

end                                                      % CÁLCULO DE LA MAGNETIZACIÓN Y TORQUE INSTANTÁNEO.