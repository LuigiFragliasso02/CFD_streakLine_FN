function [U, V] = psi2uv(PSI, U, V, G, h, Nx, Ny, uRef)
    i = 2:Nx-1;   j = 2:Ny-1;
    
    U(i, j) = (PSI(i, j+1) - PSI(i, j-1)) / (2 * h);
    V(i, j) = -(PSI(i+1, j) - PSI(i-1, j)) / (2 * h);
    
    U(G < 0) = 0;
    V(G < 0) = 0;
    
    U(1, :) = uRef;
    U(end, :) = uRef;
end