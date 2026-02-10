#include "mex.h"
#include <cmath>
#include <algorithm>

/*
 * PoissonSolver_SOR_CPP.cpp
 * Risolutore SOR in C++ per l'equazione di Poisson nel modello Psi-Zita.
 * * INPUTS (dati da MATLAB):
 * 0: Zita (Nx x Ny double)
 * 1: PSI_old (Nx x Ny double) - Initial Guess
 * 2: G (Nx x Ny double) - Topological Matrix
 * 3: PsiS (Vector double) - BC Sud
 * 4: PsiN (Vector double) - BC Nord
 * 5: PsiW (Vector double) - BC Ovest
 * 6: PsiE (Vector double) - BC Est
 * 7: PsiCyl (Scalar double) - BC Ostacolo
 * 8: hq (Scalar double) - h^2
 * 9: om (Scalar double) - Omega relaxation
 * 10: tol (Scalar double) - Tolerance
 * * OUTPUT:
 * 0: PHI (Nx x Ny double) - Updated stream function
 */

// Macro per accedere agli array lineari come matrici 2D (Column-Major order di MATLAB)
// i = riga (x), j = colonna (y), M = numero di righe (Nx)
#define IDX(i, j, M) ((i) + (j)*(M))

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    // CONTROLLO INPUT 
    if (nrhs != 11) {
        mexErrMsgIdAndTxt("Fluid:SOR:InvalidInput", "Servono 11 argomenti di input.");
    }

    // ESTRAZIONE DATI
    // Puntatori agli array di input
    double *Zita = mxGetPr(prhs[0]);
    double *PSI_In = mxGetPr(prhs[1]);
    double *G = mxGetPr(prhs[2]);
    double *PsiS = mxGetPr(prhs[3]);
    double *PsiN = mxGetPr(prhs[4]);
    double *PsiW = mxGetPr(prhs[5]);
    double *PsiE = mxGetPr(prhs[6]);

    // Scalari
    double PsiCyl = mxGetScalar(prhs[7]);
    double hq = mxGetScalar(prhs[8]);
    double om = mxGetScalar(prhs[9]);
    double tol = mxGetScalar(prhs[10]);

    // Dimensioni (Nx corrisponde alle righe M, Ny alle colonne N)
    const mwSize Nx = mxGetM(prhs[0]);
    const mwSize Ny = mxGetN(prhs[0]);

    // PREPARAZIONE OUTPUT 
    // Creiamo la matrice di uscita duplicando la PSI_old (initial guess)
    plhs[0] = mxDuplicateArray(prhs[1]);
    double *PHI = mxGetPr(plhs[0]);

    // ALGORITMO SOR
    int max_iter = 5000;
    int iter = 0;
    double normres = 1.0;

    // Variabili temporanee per i valori dei vicini
    double valN, valS, valE, valW;

    while (normres > tol && iter < max_iter) {
        iter++;
        double max_res_loc = 0.0;

        // Cicli sui nodi interni (1-based in MATLAB diventa 1..Nx-2 in C++)
        // IMPORTANTE: Loop esterno su J (colonne) e interno su I (righe) 
        // per rispettare la memoria column-major e ottimizzare la cache.
        for (mwSize j = 1; j < Ny - 1; j++) {
            for (mwSize i = 1; i < Nx - 1; i++) {

                // Indice lineare del nodo corrente
                mwSize k = IDX(i, j, Nx);

                // Se siamo nel fluido (G > 0)
                if (G[k] > 0) {

                    // --- EST (i+1) ---
                    mwSize kE = IDX(i + 1, j, Nx);
                    if (G[kE] > 0)       valE = PHI[kE];
                    else if (G[kE] == 0) valE = PsiE[j]; // PsiE(j) in Matlab
                    else                 valE = PsiCyl;

                    // --- OVEST (i-1) ---
                    mwSize kW = IDX(i - 1, j, Nx);
                    if (G[kW] > 0)       valW = PHI[kW];
                    else if (G[kW] == 0) valW = PsiW[j]; // PsiW(j) in Matlab
                    else                 valW = PsiCyl;

                    // --- NORD (j+1) ---
                    mwSize kN = IDX(i, j + 1, Nx);
                    if (G[kN] > 0)       valN = PHI[kN];
                    else if (G[kN] == 0) valN = PsiN[i]; // PsiN(i) in Matlab
                    else                 valN = PsiCyl;

                    // --- SUD (j-1) ---
                    mwSize kS = IDX(i, j - 1, Nx);
                    if (G[kS] > 0)       valS = PHI[kS];
                    else if (G[kS] == 0) valS = PsiS[i]; // PsiS(i) in Matlab
                    else                 valS = PsiCyl;

                    // --- CALCOLO SOR ---
                    double RHS_Term = hq * Zita[k];
                    double Phi_Star = 0.25 * (valE + valW + valN + valS - RHS_Term);

                    // Calcolo residuo (ogni 10 iterazioni per performance)
                    if (iter % 10 == 0) {
                        double res = std::abs(Phi_Star - PHI[k]);
                        if (res > max_res_loc) max_res_loc = res;
                    }

                    // Aggiornamento In-Place
                    PHI[k] = (1.0 - om) * PHI[k] + om * Phi_Star;
                }
            }
        }

        // Aggiorna la norma del residuo ogni 10 iterazioni
        if (iter % 10 == 0) {
            normres = max_res_loc;
        }
    }

    // (Opzionale) Debug: stampare iterazioni se serve
    // mexPrintf("Converged in %d iterations.\n", iter);
}