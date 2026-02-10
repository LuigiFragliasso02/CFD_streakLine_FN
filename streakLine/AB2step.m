function yNew = AB2step(yN, fN, fNm1, Dt)
% Adams-Bashforth 2-step explicit
% y^{n+1} = y^n + Dt*(3/2 f^n - 1/2 f^{n-1})

yNew = yN + Dt*(1.5*fN - 0.5*fNm1);

end