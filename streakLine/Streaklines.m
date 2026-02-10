% Codice per la simulazione delle equazioni di Navier-Stokes 2D nella
% formulazione vorticita'-funzione di corrente (Psi-Zita).
% Le equazioni sono integrate in una regione non connessa interna ad un
% rettangolo con un 'ostacolo' nel dominio.
%
%        ------------------------------------------------------
%        ->                                                   ->
%        ->                                                   ->
%        ->              -----                                ->
%        ->            -       -                              ->
%        ->            -       -                              ->
%        ->              -----                                ->
%        ->                                                   ->
%        ->                                                   ->
%        ------------------------------------------------------
%
%
% La collocazione delle variabili sul mesh è di tipo 'collocated', ovvero
% tutte le variabili sono collocate nei nodi di griglia, che cadono anche
% direttamente sul bordo.
%
% La indicizzazione delle variabili segue una notazione standard nella
% quale la variabile x e' associata all'indice i, mentre la variabile y e'
% associata all'indice j
%
%
%                                 Psi(i,j+1)
%                                     o
%                                     |
%                                     |
%                                     |
%                  Psi(i-1,j) o-------o-------o Psi(i+1,j)
%                                  Psi(i,j)
%                                     |
%                                     |
%                                     o
%                                 Psi(i,j-1)
close all; clear; clc;

global V U h hq Re G Nx Ny conv om tol
% ---------------------------

method = "diretto";    % Scegli il metodo per la discretizzazione spaziale: "diretto" (delsq) oppure "iterativo_sor/iterativo_sor_cpp" (SOR)
conv   = 'skew';       % Discretizzazione spaziale Schema convettivo
time_method = "AB3";   % metodo temporale: scegli: "RK4", "AB2", "AB3"
om     = 1.8;          % Parametro di rilassamento per SOR (1 < om < 2)
tol    = 1e-5;         % Tolleranza per il metodo iterativo
% ---------------------------

% Domain geometry
Ly = 1;
AR = 3;
Lx = AR * Ly;
Nhy = 79;   
Nhx = AR * Nhy;
Ny = Nhy + 1;
Nx = Nhx + 1;
y = linspace(0, Ly, Ny);
x = linspace(0, Lx, Nx);
h = x(2) - x(1);
hq = h*h;
[X, Y] = meshgrid(x, y);

% Boundary conditions su psi
Q = 1;
PsiN = Q * ones(1, Nx)';
PsiS = zeros(1, Nx)';
PsiW = linspace(0, Q, Ny)';
PsiE = PsiW;
PsiCyl = Q / 2;

% Physical parameters
T = 5;      %%%% era 30s
Re = 500;   %%%%% era 500
Cou = 0.2;  
%NB: il codice puo avere una stabilità non proprio legata all'instabilitá del metodo 
% ma instabilità legata alle BC che sono approssimate male, ci sono delle situazioni di singolarità 
% perchè all'istante inziale è come se ci fosse un vortice infinito poichè
% c'é a parete velocità 0 e il punto successivo sono diversi da 0 che creano dei picchi. 
% Cio puo essere regolarizzato con la viscosità ma all'inizio possono creare vortici locali 
% (cosa buona potrebbe essere mettere un profilo parabolico in ingresso). 
% per questo ora mettiamo un courant <<1. 
beta = 0.49;
uRef = Q / Ly;
Dt = min(Cou * h / uRef, beta * Re * hq);
Nt = ceil(T / Dt) + 1;
t = linspace(0, T, Nt);
Dtt = t(2) - t(1);

% Obstacle geometry (Cylinder)
Center = [Lx/6, Ly/2];
Radius = Ly / 8;

% % Obstacle geometry (Airfoil NACA 2418)
% Chord = 0.8 * Ly;        % Lunghezza della corda (es. 80% dell'altezza canale)
% LE = [Lx/5, Ly/2];       % Posizione del bordo d'attacco (Leading Edge)
% Alpha = 0;              % Angolo d'attacco in gradi 

% Domain topology and discrete Laplacian operator
k = 0;
G = zeros(Nx, Ny);
for j = 2 : Ny-1
    for i = 2 : Nx-1
        if norm([x(i), y(j)] - Center) <= Radius  %=if (x(i)-center(1))^2 + (y(j) - center(2))^2 < radius^2
            G(i, j) = 0;
        else
            k = k + 1;
            G(i, j) = k;
        end
    end
end

% PREPARAZIONE SOLUTORE
Lap = []; 
if method == "diretto"
    disp('Inizializzazione metodo DIRETTO (Matrice Laplaciana)...');
    Lap = -delsq(G) / hq;
else
    disp('Inizializzazione metodo ITERATIVO (SOR)...');
    % Non serve calcolare Lap per il SOR
end

% Set G=-1 inside obstacle for identification
for i = 2 : Nx-1
    for j = 2 : Ny-1
        if G(i, j) == 0
            G(i, j) = -1;
        end
    end
end

% Domain and obstacle visualization
figurePosition = [0, 50, 1640, 1640/AR];
f = figure('Position', figurePosition, 'Units', 'pixels');
D = zeros(Nx, Ny);
D(G == -1) = 1; 
D(G == 0) = 1;  
D = D';
spy(D);
xlabel('\it i'); ylabel('\it j'); 
set(gca, 'YDir', 'normal', 'FontSize', 18); 
drawnow;

% Array initialization
Ninc = k;  %numero di incognite 
zita = zeros(Ninc, 1);
PSI = zeros(Nx, Ny);
ZITA = PSI;
U = PSI;
V = PSI;

% Boundary conditions into the arrays
U(1, :) = uRef;         U(end, :) = uRef;    %velocita' costante parete O-E
PSI(1, :) = PsiW;       PSI(end, :) = PsiE;
PSI(:, 1) = PsiS;       PSI(:, end) = PsiN;  %parete SUD psi=0, NORD psi=Q
PSI(G == -1) = PsiCyl;  %=Q/2

% Pathlines initial condition
xP{1} = [0; Ly/2 + 0.100]; xP{2} = [0; Ly/2 + 0.125]; xP{3} = [0; Ly/2 + 0.150];
xP{4} = [0; Ly/2 - 0.100]; xP{5} = [0; Ly/2 - 0.125]; xP{6} = [0; Ly/2 - 0.150];

% 6 Streaklines initialization
for s=1:6
    xS{s} = zeros(2, Nt);
    xS{s}(:, 1) = xP{s};                % inizializziamo la posizione all't=0
    xS{s}(:, :) = repmat(xP{s}, 1, Nt); %permette di copiare la prima colonna inizializzata in tutte le altre
end

% Graphics preparation
clf(f);     %pulisce figura
Uquiv = zeros(Nx, Ny); Vquiv = zeros(Nx, Ny);  % matrice che conterrà U per il grafico a frecce
iq = 1 : 3 : Nx;   % si selezione un punto per la freccia ogni 3 lungo la direzione x altrimenti troppe frecce
jq = [1 : 6 : floor(Ny/3), floor(Ny/3+3) : 3 : floor(2*Ny/3), floor(2*Ny/3+6) : 6 : Ny]; %campionamento qui variabile 
%per focalizzare l'attenzione dove serve, cioè fascia bassa si prende un punto ogni 6, fascia centrale ogni 3 e
%fascia alta ogni 6

% mappa per i colori di vorticità - bianco zita=0, blu per zita<0, rosso per zita>0
Map = [linspace(75/255, 1, 500)', linspace(90/255, 1, 500)', linspace(241/255, 1, 500)'; ... %sfumatura dal bianco al blu
       linspace(1, 180/255, 500)', linspace(1, 10/255, 500)', linspace(1, 32/255, 500)'];    %sfumatura dal rosso al blu

% Variabili per analisi
enstrophy = NaN(Nt, 1);        %per valutare la conservatività dell'estrofia
time_per_step = zeros(Nt, 1);  %per valutare la velocità del metodo


%%modifica per AB
RHS_hist = cell(3,1); 
% RHS_hist{1} = f^{n-2}
% RHS_hist{2} = f^{n-1}
% RHS_hist{3} = f^{n}


% Time integration 
disp('Avvio integrazione temporale...');
for it = 1 : Nt
    time = t(it);
    
    ZITA(G == -1) = 0; 
    
    % Vorticità al bordo (Thom)
    ZITA = ThomFormulae(ZITA, PSI, G, hq);    %cc sulla zita
    

    % Integrazione Eq. Vorticità (RK4  - AB2 - AB3)

    % === Calcolo RHS corrente (f^n) ===
    RHSn = ZitaRHS(time, ZITA, Re, h, hq, U, V, G, Nx, Ny);
    
    switch time_method
    
        case "RK4"
            ZITA = RK4(time, ZITA, @(t,y) ZitaRHS(t,y,Re,h,hq,U,V,G,Nx,Ny), Dt);
    
        case "AB2"
            if it <= 2
                % Startup con RK4
                ZITA = RK4(time, ZITA, @(t,y) ZitaRHS(t,y,Re,h,hq,U,V,G,Nx,Ny), Dt);
    
                % Riempio storico progressivamente
                RHS_hist{it} = RHSn;
            else
                % AB2: usa f^n e f^{n-1}
                ZITA = AB2step(ZITA, RHSn, RHS_hist{2}, Dt);
    
                % Shift storico
                RHS_hist{1} = RHS_hist{2};
                RHS_hist{2} = RHSn;
            end
    
        case "AB3"
            if it <= 3
                % Startup con RK4
                ZITA = RK4(time, ZITA, @(t,y) ZitaRHS(t,y,Re,h,hq,U,V,G,Nx,Ny), Dt);
    
                RHS_hist{it} = RHSn;
            else
                % AB3: usa f^n, f^{n-1}, f^{n-2}
                ZITA = AB3step(ZITA, RHSn, RHS_hist{3}, RHS_hist{2}, Dt);
    
                % Shift storico
                RHS_hist{1} = RHS_hist{2};
                RHS_hist{2} = RHS_hist{3};
                RHS_hist{3} = RHSn;
            end
    
        otherwise
            error("Metodo temporale non riconosciuto");
    end


    % Risoluzione Ellittica per PSI (Poisson)
    tic; % Start timer risolutore
    
    if method == "diretto"
        % Metodo Diretto (si usa delsq classico)
        [PSI, zita] = zita2psi(zita, ZITA, PSI, Lap, G, Nx, Ny, hq, PsiS, PsiN, PsiW, PsiE, PsiCyl);
    elseif method == "iterativo_sor"
        % Metodo Iterativo (SOR)
        % Nota: Non restituisce 'zita' vettore, ma aggiorna direttamente PSI matrice
        PSI = PoissonSolver_SOR(ZITA, PSI, G, PsiS, PsiN, PsiW, PsiE, PsiCyl);
    else
        PSI = PoissonSolver_SOR_CPP(ZITA, PSI, G, PsiS, PsiN, PsiW, PsiE, PsiCyl, hq, om, tol);
    end
    
    time_per_step(it) = toc; % Stop timer
    
    % Calcolo Enstrofia Globale (uguale per entrambi i casi)
    % Integrale di zita^2 sul dominio fluido
    enstrophy(it) = 0.5 * sum(ZITA(G > 0).^2) * hq;
    
    % Calcolo Velocità
    [U, V] = psi2uv(PSI, U, V, G, h, Nx, Ny, uRef);
    
    %  Aggiornamento Streaklines
    for k = 1 : it
        for s=1:6
            xS{s}(:, k) = RK4(time, xS{s}(:, k), @(t, y_) PathlinesRHS(t, y_, x, y, U, V, h, h, Nx, Ny), Dt);
        end
    end
    

    % Plotting ogni 10 step
    if mod(it, 10) == 0
        ZITAg = ZITA;
        ZITAg(G == -1) = NaN;
        
        subplot(2,1,1); % Plot fluido sopra
        pcolor(x, y, ZITAg');
        colormap(Map);
        shading interp;
        
        vals = ZITA(G > 0); vals = vals(~isnan(vals));
        if ~isempty(vals)
            c_min = 0.7 * min(vals); c_max = 0.7 * max(vals);
            if c_min < c_max; clim([c_min, c_max]); end
        end
        colorbar; hold on;
        
        Uquiv(iq, jq) = U(iq, jq);
        Vquiv(iq, jq) = V(iq, jq);
        quiver(X(jq, iq), Y(jq, iq), Uquiv(iq, jq)', Vquiv(iq, jq)', 3, 'Color', 'k');

        % %per profilo naca
        % contour(X, Y, G', [0 0], 'k', 'LineWidth', 2);
        
        % Disegno ostacolo generico
        contour(X, Y, G', [0 0], 'k', 'LineWidth', 2);
        
        % Plot Streaklines
        indices_plot = linspace(1, it, min(it*10, 2000));
        for s = 1:6
             if it > 1
                % Filtro NaN per evitare crash interpolazione
                Xpts = xS{s}(1, 1:it); Ypts = xS{s}(2, 1:it);
                if any(isnan(Xpts)) || any(isnan(Ypts)); continue; end
                
                x_plot = interp1(1:it, Xpts, indices_plot, 'pchip');
                y_plot = interp1(1:it, Ypts, indices_plot, 'pchip');
                col = 'r'; if s > 3; col = 'b'; end
                plot(x_plot, smooth(y_plot), 'LineWidth', 1.5, 'Color', col);
             end
        end
        
        title(['t = ', num2str(time, 2), ' | ', char(method), ' | Solver Time: ', num2str(time_per_step(it), '%.4f'), 's']);
        axis image; axis([0, Lx, 0, Ly]);
        hold off;
        
        subplot(2,1,2); % Plot enstrofia sotto
        plot(t(1:it), enstrophy(1:it), 'r', 'LineWidth', 2.5);
        title('Evoluzione Enstrofia Globale');
        xlabel('t'); ylabel('\Omega'); grid on; xlim([0, T]);
        
        drawnow limitrate;
    end
end

avg_time = mean(time_per_step);
fprintf('Simulazione conclusa. Metodo: %s. Tempo medio solver: %.5f s\n', method, avg_time);