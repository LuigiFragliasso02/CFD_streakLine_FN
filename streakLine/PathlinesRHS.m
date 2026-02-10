function dxdt = PathlinesRHS(t, xP, x, y, U, V, Dx, Dy, Nx, Ny)
    
    %localizzazione deglla cella in cui si trova la particella 
    i = floor(xP(1) / Dx) + 1;
    j = floor(xP(2) / Dy) + 1;
    
    % Ensure indices are within bounds
    i = max(1, min(i, Nx-1));
    j = max(1, min(j, Ny-1));

    %coordinate del nodo in basso a sinistra deella cella in cui si trova
    %la particella 
    xNode = [x(i); y(j)];
    
    %adimensionalizzazione della posizione della particella q nella cella
    csix = (xP(1) - xNode(1)) / Dx;
    csiy = (xP(2) - xNode(2)) / Dy;
    
    if i+1 <= Nx && j+1 <= Ny && i > 0 && j > 0
        % Bilinear interpolation
        u = [1-csix, csix] * U([i, i+1], [j, j+1]) * [1-csiy; csiy];
        v = [1-csix, csix] * V([i, i+1], [j, j+1]) * [1-csiy; csiy];
    else
        u = 0.05; % Minimal flow backup
        v = 0;
    end
    
    dxdt = [u; v];
end

%interp2(Xm, Ym, U_transp, xx, yy, 'linear', 0)  tutto cio si riassume in
%questa funzione