function F = ZitaRHS(t, ZITA, Re, h, hq, U, V, G, Nx, Ny)
% function  F  = ZitaRHS(~,ZITA)
% global V U h hq Re G Nx Ny conv

% RHS della equazione per Zita per il codice TunnelPsiZita

global conv
F = zeros(Nx,Ny);
for i = 2:Nx-1
    for j = 2:Ny-1
% Convective term
        if G(i,j)>0
            switch conv
                case 'div'
                    F(i,j) = -((U(i+1,j)*ZITA(i+1,j)-U(i-1,j)*ZITA(i-1,j))/(2*h)+...
                               (V(i,j+1)*ZITA(i,j+1)-V(i,j-1)*ZITA(i,j-1))/(2*h) );
                case 'adv'
                    F(i,j) = -( U(i,j)*(ZITA(i+1,j)-ZITA(i-1,j))/(2*h)+...
                                V(i,j)*(ZITA(i,j+1)-ZITA(i,j-1))/(2*h) );
                case 'skew'
                    F(i,j) = -0.5*((U(i+1,j)*ZITA(i+1,j)-U(i-1,j)*ZITA(i-1,j))/(2*h)+...
                                   (V(i,j+1)*ZITA(i,j+1)-V(i,j-1)*ZITA(i,j-1))/(2*h) );
                    F(i,j) = F(i,j)-0.5*(U(i,j)*(ZITA(i+1,j)-ZITA(i-1,j))/(2*h)+...
                                         V(i,j)*(ZITA(i,j+1)-ZITA(i,j-1))/(2*h) );
            end
% Diffusive term
            F(i,j) = F(i,j) + (1/Re)*(ZITA(i+1,j)+ZITA(i-1,j)+ZITA(i,j+1)+...
                                      ZITA(i,j-1) - 4*ZITA(i,j))/hq;
        end
    end
end
end

