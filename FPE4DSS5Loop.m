%Copyright 2013 Rice University
%All Rights Reserved
%
%Created by Ali Arya Mokhtarzadeh
%Department of Mechanical Engineering
%

%%Loop
%The Optimization Loop
while(LR > l.tol && count < l.Nmax)
    
    %Stepping to the Next Basis
    count = count + 1;
    
    %Find the Optimization Parameters
    [par fval] = ga(@(par)FPE4DSS6OptSymD2(par,p,d,v),d.numpar,gao);
    %[par fval] = patternsearch(@(par) FPE1DSSnd6OptSymD2(par,p,d,v),rand(d.numpar,1),[],[],[],[],[],[],[],pso);
    %[par fval] = patternsearch(@(par) FPE4DOptSSndSym(par,p,d,v),rand(6,1),[],[],[],[],[],[],[],pso);
    
    %Updating the Approximation
    if (  isnan(fval)  ) %(abs(par(2)) <0.1) || (abs(par(2)) > 5) || abs(-100*par(1)^2) > 100 || 
        count = count - 1;
    else
        FPE4DSS7UpdateSym
    end
    
    %Updating the problem
    %if (fval > 0);
    %    FPE4DUpdateSSnd;
    %else
    %    count = count-1;
    %end
    
end

toc

figure;
%ERROR ANALYSIS
%subplot(2,2,1); plot(counter,rmshist);
%title('rms error (approximation to true solution) vs. basis function')
%subplot(2,2,2); loglog(counter,rmshist);
%title('rms error (approximation to true solution) vs. basis function')
%subplot(2,2,3); plot(counter,LRhist);
%title('norm of the pde residual vector vs. basis function')
%subplot(2,2,4); loglog(counter,LRhist);
%title('norm of the pde residual vector vs. basis function')