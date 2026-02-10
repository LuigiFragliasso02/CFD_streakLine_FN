function [PSI, zita] = zita2psi(zita, ZITA, PSI, Lap, G, Nx, Ny, hq, PsiS, PsiN, PsiW, PsiE, PsiCyl)

    % Inizializzazione: Mappatura della vorticità (ZITA 2D) nel vettore termine noto (zita 1D)
    for i = 2:Nx-1
        for j = 2:Ny-1
            if G(i,j) > 0
                % Mettiamo il valore di vorticità nel punto k corrispondente
                zita(G(i,j)) = ZITA(i,j);  
            end
        end
    end

    % Modifica del termine noto per le Condizioni al Contorno (BCs)
    % Ciclo su tutti i nodi interni per controllare i vicini
    for i = 2:Nx-1
        for j = 2:Ny-1
            if G(i,j) > 0       % Se siamo in un nodo fluido
                k = G(i,j);     % Indice lineare nel sistema risolutivo
                
                % --- Controllo Pareti Esterne (G == 0) ---
                
                % Parete Sud (vicino j-1 è bordo 0)
                if G(i,j-1) == 0; zita(k) = zita(k) - PsiS(i)/hq; end 
                % Parete Nord (vicino j+1 è bordo 0)
                if G(i,j+1) == 0; zita(k) = zita(k) - PsiN(i)/hq; end 
                % Parete Ovest (vicino i-1 è bordo 0)
                if G(i-1,j) == 0;  zita(k) = zita(k) - PsiW(j)/hq; end 
                % Parete Est (vicino i+1 è bordo 0)
                if G(i+1,j) == 0;  zita(k) = zita(k) - PsiE(j)/hq; end 

                % Controllo Pareti Ostacolo (G == -1)
                % Se un vicino è ostacolo, sottraiamo il valore costante PsiCyl
                
                if G(i,j+1) == -1;  zita(k) = zita(k) - PsiCyl/hq;     end    
                if G(i,j-1) == -1;  zita(k) = zita(k) - PsiCyl/hq;     end    
                if G(i+1,j) == -1;  zita(k) = zita(k) - PsiCyl/hq;     end 
                if G(i-1,j) == -1;  zita(k) = zita(k) - PsiCyl/hq;     end 
            end
        end
    end
    
    % Risoluzione dell'equazione di Poisson (Sistema Lineare) con metodo diretto
    % Lap è la matrice del Laplaciano discreto, zita è il termine noto modificato
    psi = Lap \ zita;

    % Riversamento del risultato nel vettore 2D PSI
    PSI(G > 0) = psi;
end