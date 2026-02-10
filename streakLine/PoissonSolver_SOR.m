
function PHI = PoissonSolver_SOR(Zita, PSI_old, G, PsiS, PsiN, PsiW, PsiE, PsiCyl)
% Risolutore iterativo SOR per l'equazione di Poisson nabla^2 psi = zita
% Gestisce le BCs di Dirichlet non omogenee e l'ostacolo interno.

    global Nx Ny hq om tol

    PHI = PSI_old; % Guess iniziale = soluzione al passo precedente
    iter = 0;
    normres = 1;
    max_iter = 2000; % Limite sicurezza per evitare loop infiniti

    while normres > tol && iter < max_iter
        iter = iter + 1;
        max_res = 0; % Per calcolare la norma infinito del residuo

        for i = 2:Nx-1
            for j = 2:Ny-1
                if G(i,j) > 0 % Calcola solo sui nodi fluidi

                    %  Recupero valori dai vicini (Gestione BCs), in particolare se il vicino 
                    % è un punto del bordo allora prenditi il valore noto di PSI altrimenti 
                    % prendiamo il valore in quel puunto di PSI. 

                    % EST (i+1)
                    if G(i+1,j) > 0;      valE = PHI(i+1,j);   % Fluido
                    elseif G(i+1,j) == 0; valE = PsiE(j);      % Bordo Est
                    else;                 valE = PsiCyl;       % Ostacolo
                    end

                    % OVEST (i-1)
                    if G(i-1,j) > 0;      valW = PHI(i-1,j);   % Fluido
                    elseif G(i-1,j) == 0; valW = PsiW(j);      % Bordo Ovest
                    else;                 valW = PsiCyl;       % Ostacolo
                    end

                    % NORD (j+1)
                    if G(i,j+1) > 0;      valN = PHI(i,j+1);   % Fluido
                    elseif G(i,j+1) == 0; valN = PsiN(i);      % Bordo Nord
                    else;                 valN = PsiCyl;       % Ostacolo
                    end

                    % SUD (j-1)
                    if G(i,j-1) > 0;      valS = PHI(i,j-1);   % Fluido
                    elseif G(i,j-1) == 0; valS = PsiS(i);      % Bordo Sud
                    else;                 valS = PsiCyl;       % Ostacolo
                    end

                    % Formula di Aggiornamento SOR 
                    % Laplaciano discreto: (E + W + N + S - 4P)/h^2 = Zita
                    % P_new = (1-om)*P_old + om * 0.25 * (E + W + N + S - h^2*Zita)

                    RHS_Term = hq * Zita(i,j);
                    Phi_Star = 0.25 * (valE + valW + valN + valS - RHS_Term);

                    % Calcolo residuo locale (solo ogni 10 iterazioni per velocità)
                    if mod(iter, 10) == 0
                        res_loc = abs(Phi_Star - PHI(i,j)); % Semplificazione del residuo
                        if res_loc > max_res; max_res = res_loc; end
                    end
                    %perchè calcoliamo così il residuo?
                    % In ogni singolo istante, guardi i tuoi vicini e calcoli: "Se i miei vicini restassero fermi così come sono ora, quale valore dovrei avere io per rispettare la regola?"
                    % Quel valore calcolato è Phi_Star. È un target locale momentaneo. Non è la soluzione finale assoluta, è solo il valore che ti metterebbe "in pace" con i tuoi vicini adesso.
                    % 2. Perché il bersaglio si muove?
                    % Il problema è che, appena tu ti sposti verso il tuo Phi_Star, anche i tuoi vicini stanno facendo lo stesso calcolo e si spostano verso il loro Phi_Star. 
                    % Cambiando loro, cambia la loro media, e quindi cambia di nuovo il tuo obiettivo.
                    % Iterazione 1: I vicini sono a 0. Il tuo target è X. Ti muovi verso X.
                    % Iterazione 2: I vicini si sono mossi. La loro media è cambiata. Il tuo target ora è Y. Ti muovi verso Y.
                    % Convergenza: Dopo tanti passaggi, nessuno si muove più significativamente. Il "target" (Phi_Star) coincide con dove sei (PHI). Hai raggiunto l'equilibrio.


                    % Aggiornamento
                    PHI(i,j) = (1 - om) * PHI(i,j) + om * Phi_Star;

                end
            end
        end
        
        %il prof calcola il residuo in questo modo cioè semplciemente la
        %differenza tra lap psi - zita e se è sotto una certa sogni allora fermati
        % if mod(iter,10) == 0
        % i = 2:Nx-1;    j = 2:Ny-1;
        % Res(i-1,j-1) = (PHI(i+1,j) + PHI(i-1,j) + PHI(i,j-1) + PHI(i,j+1) - 4*PHI(i,j) )/hq - Zita(i,j);   %residuo calcolato come la differenza tra lap psi - zita = 0
        % normres = max(abs(Res(:)));

        % Aggiorna residuo globale
        if mod(iter, 10) == 0
            normres = max_res;
        end
    end
end
