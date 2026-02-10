function yNew = RK4(t, yOld, RHSfunction, Dt)
    t1 = t;
    y1 = yOld;
    f1 = RHSfunction(t1, y1);
    
    t2 = t + Dt / 2;
    y2 = yOld + Dt * 0.5 * f1;
    f2 = RHSfunction(t2, y2);
    
    t3 = t + Dt / 2;
    y3 = yOld + Dt * 0.5 * f2;
    f3 = RHSfunction(t3, y3);
    
    t4 = t + Dt;
    y4 = yOld + Dt * f3;
    f4 = RHSfunction(t4, y4);
    
    yNew = yOld + Dt * ((f1 + f4)./6 + (f2 + f3)./3);
end