% function isInside = isInsideNACA2418(x, y, c, le, alpha_deg)
%     isInside = false;
%     m = 0.02;
%     p = 0.40;
%     t = 0.18; % thickness ratio relative to chord (adjusted from paper's t=0.18*c context)
% 
%     % Helper functions for camber line
%     yc_fun = @(xc) (m/p^2 * (2*p*xc - xc.^2)) .* (xc >= 0 & xc <= p) + ...
%                    (m/(1-p)^2 * ( (1-2*p) + 2*p*xc - xc.^2 )) .* (xc > p & xc <= 1);
% 
%     dyc_dx_fun = @(xc) (2*m/p^2 * (p - xc)) .* (xc >= 0 & xc <= p) + ...
%                        (2*m/(1-p)^2 * (p - xc)) .* (xc > p & xc <= 1);
% 
%     theta_fun = @(xc) atan(dyc_dx_fun(xc));
% 
%     % Thickness distribution
%     yt_fun = @(xc) 5 * t * c * (0.2969*sqrt(xc) - 0.1260*xc - 0.3516*xc.^2 + ...
%                               0.2843*xc.^3 - 0.1015*xc.^4);
% 
%     % Coordinate transformation (simplified logic for point check)
%     % Inverse rotation to align with airfoil chord
%     dx = x - le(1);
%     dy = y - le(2);
%     xc_rot = dx * cosd(-alpha_deg) - dy * sind(-alpha_deg);
%     yc_rot = dx * sind(-alpha_deg) + dy * cosd(-alpha_deg);
% 
%     xc = xc_rot / c; % Normalized chord position
% 
%     if xc >= 0 && xc <= 1
%         th = theta_fun(xc);
%         yt = yt_fun(xc);
%         yc = yc_fun(xc) * c; % Scale camber by chord
% 
%         yU = yc + yt * cos(th);
%         yL = yc - yt * cos(th);
% 
%         % Check if the point's y-coordinate (in rotated frame) is within Upper and Lower surfaces
%         if yc_rot >= yL && yc_rot <= yU
%             isInside = true;
%         end
%     end
% end



function isInside = isInsideNACA2418(x, y, c, le, alpha_deg)
    % isInsideNACA2418: Verifica se un punto (x,y) è dentro il profilo.
    % x, y: coordinate del punto della griglia
    % c: lunghezza della corda
    % le: vettore [x,y] del bordo d'attacco (Leading Edge)
    % alpha_deg: angolo di attacco in gradi
    
    isInside = false;
    
    % Parametri NACA 2418
    m = 0.02;       % Massima curvatura (2%)
    p = 0.40;       % Posizione mass. curvatura (40%)
    t = 0.18;       % Spessore massimo (18%)
    
    % --- 1. Trasformazione Coordinate (Dal Globale al Locale) ---
    % Traslazione rispetto al Leading Edge
    dx = x - le(1);
    dy = y - le(2);
    
    % Rotazione inversa (riportiamo il punto sull'asse orizzontale)
    % Se l'ala è ruotata di +alpha, ruotiamo il punto di -alpha
    xc_rot = dx * cosd(-alpha_deg) - dy * sind(-alpha_deg);
    yc_rot = dx * sind(-alpha_deg) + dy * cosd(-alpha_deg);
    
    % Normalizzazione sulla corda (xc va da 0 a 1)
    xc = xc_rot / c; 
    
    % --- 2. Check Limiti del Profilo ---
    % Se siamo prima del bordo d'attacco o dopo il bordo d'uscita, è falso
    if xc < 0 || xc > 1
        return;
    end
    
    % --- 3. Calcolo Geometria NACA ---
    % Linea media (Camber Line) adimensionale
    if xc <= p
        yc_adim = (m / p^2) * (2*p*xc - xc^2);
        dyc_dx  = (2*m / p^2) * (p - xc);
    else
        yc_adim = (m / (1-p)^2) * ((1-2*p) + 2*p*xc - xc^2);
        dyc_dx  = (2*m / (1-p)^2) * (p - xc);
    end
    
    theta = atan(dyc_dx);
    
    % Distribuzione spessore (tutto dimensionale qui, moltiplichiamo per c)
    % Nota: usiamo real(sqrt) per sicurezza contro errori di arrotondamento negativi
    yt = 5 * t * c * (0.2969*real(sqrt(xc)) - 0.1260*xc - 0.3516*xc^2 + ...
                      0.2843*xc^3 - 0.1015*xc^4);
    
    % Linea media dimensionale
    yc = yc_adim * c;
    
    % Calcolo coordinate superiore (Upper) e inferiore (Lower) locali
    % Nota: approssimiamo la x locale alla xc per il check verticale, 
    % che per la generazione di maschere (G) è sufficiente.
    yU = yc + yt * cos(theta);
    yL = yc - yt * cos(theta);
    
    % --- 4. Verifica Finale ---
    % Controlliamo se la y del punto (nel sistema locale) è tra dorso e ventre
    if yc_rot >= yL && yc_rot <= yU
        isInside = true;
    end
end