function yNew = AB3step(yN, fN, fNm1, fNm2, Dt)
% Adams-Bashforth 3-step explicit
% y^{n+1} = y^n + Dt*(23/12 f^n - 16/12 f^{n-1} + 5/12 f^{n-2})

yNew = yN + Dt*( (23/12)*fN - (16/12)*fNm1 + (5/12)*fNm2 );

end