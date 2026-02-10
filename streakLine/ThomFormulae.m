function ZITA = ThomFormulae(ZITA, PSI, G, hq)
    [Nx, Ny] = size(ZITA);
    
    % South Wall
    i = 2:Nx-1;     j = 1;
    ZITA(i, j) = 2 * (PSI(i, j+1) - PSI(i, j)) / hq;   %u_wall=0
    
    % North Wall
    i = 2:Nx-1;       j = Ny;
    ZITA(i, j) = 2 * (PSI(i, j-1) - PSI(i, j)) / hq;
    
    % Obstacle Walls
    for i = 2 : Nx-1
        for j = 2 : Ny-1
            if G(i, j) == -1 && G(i+1, j) > 0 % East neighbor fluid
                ZITA(i, j) = 2 * (PSI(i+1, j) - PSI(i, j)) / hq;
            elseif G(i, j) == -1 && G(i-1, j) > 0 % West neighbor fluid
                ZITA(i, j) = 2 * (PSI(i-1, j) - PSI(i, j)) / hq;
            elseif G(i, j) == -1 && G(i, j-1) > 0 % South neighbor fluid
                ZITA(i, j) = 2 * (PSI(i, j-1) - PSI(i, j)) / hq;
            elseif G(i, j) == -1 && G(i, j+1) > 0 % North neighbor fluid
                ZITA(i, j) = 2 * (PSI(i, j+1) - PSI(i, j)) / hq;
            end
        end
    end
end