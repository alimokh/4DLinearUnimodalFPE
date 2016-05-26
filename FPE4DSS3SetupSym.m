%Copyright 2013 Rice University
%All Rights Reserved
%
%Created by Ali Arya Mokhtarzadeh
%Department of Mechanical Engineering
%

casenumber = 1;
v = FPE4DInput(casenumber);
casenumber = 2;
d = FPE4DInput(casenumber);

%Create the Domain 
%[x1 x2] = ndgrid(d.lowx1:d.disc:d.upx1, d.lowx2:d.disc:d.upx2);
[x1 x2 x3 x4] = ndgrid(d.lowx1:d.disc:d.upx1, d.lowx2:d.disc:d.upx2, d.lowx3:d.disc:d.upx3, d.lowx4:d.disc:d.upx4);
d.x1 = x1;
d.x2 = x2;
d.x3 = x3;
d.x4 = x4;
d.firstlength = length(d.x1);
d.qq = length(d.x1);
d.mid = (d.firstlength+1)/2;

%INITIALIZE true solution, Loop and Optimization Information
casenumber = 3;
pe = FPE4DInput(casenumber,v,d);
casenumber = 4;
l = FPE4DInput(casenumber);


count = 0;
gao = gaoptimset('UseParallel', 'always', 'Vectorized', 'off','MutationFcn',@mutationadaptfeasible,'Display','off','PopInitRange',[0 0 0 0 0 0; 1 1 1 1 1 1],'PopulationSize',40,'stallgenlimit',50);%%'PopulationSize',40)   'UseParallel', 'always', 'Vectorized', 'off',
%gao = gaoptimset('UseParallel', 'always', 'Vectorized', 'off','Display','off','PopInitRange',[0 0; 1 1],'PopulationSize',80,'stallgenlimit',50);%%'PopulationSize',40)   'UseParallel', 'always', 'Vectorized', 'off',
%gao = gaoptimset('UseParallel', 'always', 'Vectorized', 'off','MutationFcn',@mutationadaptfeasible,'Display','off','PopInitRange',[0 0 -10 -10; 1 1 10 10],'PopulationSize',40,'stallgenlimit',50);%Spencer Bergman 1991
pso = psoptimset('Display','off','CompletePoll','on','CompleteSearch','on','UseParallel', 'always', 'Vectorized', 'off');

%Define symbolic variables and function (must satisfy the problem IC)  
casenumber = 5;
x = FPE4DInput(casenumber,v,d);
casenumber = 6;
p = FPE4DInput(casenumber,v,d,x);

%%Take the numerical derivatives
%p.p(1,:) = 0;p.p(:,1) = 0;p.p(length(p.p),:) = 0;p.p(:,length(p.p)) = 0;

[p.px2 p.px1 p.px3 p.px4] = gradient(p.p,d.disc);
[p.px2x2 p.px2x1 p.px2x3 p.px2x4] = gradient(p.px2,d.disc);
[p.px4x2 p.px4x1 p.px4x3 p.px4x4] = gradient(p.px4,d.disc);

%[p.px2 p.px1] = gradient2(p.p,d.disc,d.disc);
%[p.px1x2 p.px1x1] = gradient2(p.px1,d.disc,d.disc);
%[p.px2x2 p.px2x1] = gradient2(p.px2,d.disc,d.disc);

%Combining the derivatives to form the residual of the PDE
casenumber = 7;
p = FPE4DInput(casenumber,v,d,x,p);

%Removing the "ghost points" before the Residual calculation
x1LR = d.x1(2:d.firstlength-1,2:d.firstlength-1);
x2LR = d.x2(2:d.firstlength-1,2:d.firstlength-1);
p.e = p.e(2:d.firstlength-1,2:d.firstlength-1);
p.p = p.p(2:d.firstlength-1,2:d.firstlength-1);
pe = pe(2:d.firstlength-1,2:d.firstlength-1);

% norm of the PDE residual vector
LR = norm(p.e(:));
disp('  LR');disp(LR);

p.initial = p.p;
