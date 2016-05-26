%Copyright 2013 Rice University
%All Rights Reserved
%
%Created by Ali Arya Mokhtarzadeh
%Department of Mechanical Engineering
%

o.w = -100*par(1)^2;%-100*par(1)^2;
%o.w = exp(-1*par(1)^2);%-100*par(1)^2;
%o.w = par(1);
o.c = (1*(par(2))-0.5);
o.x1c = (d.upx1-d.lowx1)*par(3)+d.lowx1;%par(3);
o.x2c = (d.upx1-d.lowx1)*par(4)+d.lowx1;%par(4);
%o.c = par(2);
%o.x1c = par(3);
%o.x2c = par(4);

%newbleh = abs(p.e);
%[maxLRvalue maxLRloc] = max(newbleh(:));
%[x1c_row,x2c_column] = ind2sub(size(newbleh),maxLRloc);
%o.x1c = d.x1(x1c_row);
%o.x2c = d.x1(x2c_column);

x1plotc(count) = o.x1c;%x1c = ucbc.x1(count);
%%%%%%%%%      Numerical Derivatives - No Boundary Enforcers!!! %%%%%
%syms x x2
x = d.x1;

casenumber = 8;
o = FPE4DInput(casenumber,v,d,x,p,o);
%funphi = (w^((x-x1c)^2));

%o.funphi(1,:) = 0;o.funphi(:,1) = 0;o.funphi(length(o.funphi),:) = 0;o.funphi(:,length(o.funphi)) = 0;
%o.funphi = (o.funphi).*cos((pi*d.x1)/16).*cos((pi*d.x2)/16);%.^2; %.*d.x1.*d.x2;                                                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%THIS IS WHAT DID IT

[o.phix2 o.phix1] = gradient2(o.funphi,d.disc,d.disc);
[o.phix1x2 o.phix1x1] = gradient2(o.phix1,d.disc,d.disc);
[o.phix2x2 o.phix2x1] = gradient2(o.phix2,d.disc,d.disc);

%o.funphi(1,:) = 0;o.funphi(:,1) = 0;o.funphi(length(o.funphi),:) = 0;o.funphi(:,length(o.funphi)) = 0;
%o.phix1(1,:) = 0;o.phix1(:,1) = 0;o.phix1(length(o.phix1),:) = 0;o.phix1(:,length(o.phix1)) = 0;
%o.phix2(1,:) = 0;o.phix2(:,1) = 0;o.phix2(length(o.phix2),:) = 0;o.phix2(:,length(o.phix2)) = 0;
%o.phix2x2(1,:) = 0;o.phix2x2(:,1) = 0;o.phix2x2(length(o.phix2x2),:) = 0;o.phix2x2(:,length(o.phix2x2)) = 0;

%[o.phix2 o.phix1] = gradient2(o.funphi,d.disc,d.disc);
%[o.phix1x2 o.phix1x1] = gradient2(o.phix1,d.disc,d.disc);
%[o.phix2x2 o.phix2x1] = gradient2(o.phix2,d.disc,d.disc);
%o.phix1 = 2*(d.x1-o.x1c).*log(o.w).*o.w.^((((d.x1-o.x1c).^2)+(d.x2-o.x2c).^2));
%o.phix2 = 2*(d.x2-o.x2c).*log(o.w).*o.w.^((((d.x1-o.x1c).^2)+(d.x2-o.x2c).^2));
%o.phix2x2 = 2*log(o.w).*(2*((d.x2-o.x2c).^2).*log(o.w)+1).*o.w.^((((d.x1-o.x1c).^2)+(d.x2-o.x2c).^2));

casenumber = 9;
o = FPE4DInput(casenumber,v,d,x,p,o);

%eps1 = dot(o.funpsi(3:d.firstlength-2),p.e(:));eps2 = dot(o.funpsi(3:d.firstlength-2),o.funpsi(3:d.firstlength-2));eps3 = dot(p.e(:),p.e(:));%LR = 1-eps1^2/(eps2*eps3);
%c = -eps1/eps2;

%o.funphi = max(0,o.funphi); %Matt's idea for positivity
newbasis = o.c*o.funphi;%%%%%%%%%%%Turn ON for RBF
%newbasis = o.funphi;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Only for Piecewise Basis Functions because o.c is already used to determine their height

newbasis = newbasis(2:d.firstlength-1,2:d.firstlength-1);
p.p = p.p + newbasis;                                                   %%%%%%%%%%%USE OF abs() TOTALLY HEURISTIC MAYBE NOT PHYSICALLY MEANINGFUL!!!!!
%p.p = real(exp(p.p));

newres = o.c*o.funpsi;%%%%%%%%%%%Turn ON for RBF
%newres = o.funpsi;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Only for Piecewise Basis Functions because o.c is already used to determine their height

newres = newres(2:d.firstlength-1,2:d.firstlength-1);
p.e = p.e + newres;                                                     %%%%%%%%%%%USE OF abs() TOTALLY HEURISTIC MAYBE NOT PHYSICALLY MEANINGFUL!!!!!

%LR = norm(p.e(:));% differential equation RMS


%LR = sum((p.e(:).^2));
LR = norm(p.e(:));
%p.testerr = p.e(12:d.firstlength-11,12:d.firstlength-11);
%LR = norm(p.testerr(:));
%LR = max(abs(p.e(:)));%LR = trapz(d.x1,p.e.^2);

%LR = sqrt( mse(p.pe - p.p) );

%pe(d.mid-4) = 0.6307;
%p.p(d.mid-4) = 0.6307;

%NORMALIZE THE INTEGRAL TO 1
%peintpdf = trapz(x1LR,pe);
%penormconst = 1/peintpdf;
%penew = pe*penormconst;
%intpdf = trapz(x1LR,p.p);
%normconst = 1/intpdf;
%perr = p.p*normconst;

%penew = real(pe);
%perr = real(p.p);

%record results
%rms_error = sqrt( mse(perr - penew) );
%rms_error = sqrt( mse(p.p - pe) );
%rmshist(count) = rms_error;
LRhist(count) = LR;
counter(count) = count;
chist(count) = o.c;
whist(count) = o.w;
x1hist(count) = o.x1c;
x2hist(count) = o.x2c;
%covariance
%truecovar(count) = cov(penew);
%appcovar(count) = cov(perr);
%display results
disp('count');disp(count);
disp('  w');disp(o.w);
disp('  c');disp(o.c);
disp('  LR');disp(LR);
%disp('  RMS');disp(rms_error);


%plot results
%FPE1DSSnd9Plot

%subplot(3,3,1);mesh(x1LR,x2LR,p.initial);title('Initial Condition');
subplot(3,3,2);mesh(x1LR,x2LR,p.p);title('Approximation');
subplot(3,3,3);mesh(x1LR,x2LR,p.e.^2);title('Equation Residual Squared');
subplot(3,3,4);mesh(x1LR,x2LR,newbasis);title('Current Basis Function');
subplot(3,3,5);plot(counter,LRhist,'-+');ylabel('LR');
subplot(3,3,6);loglog(counter,LRhist,'-+');ylabel('LR');
%subplot(3,3,7);loglog(counter,rmshist,'-*');ylabel('rms');
subplot(3,3,8);plot(counter,chist,'-^');ylabel('c');
subplot(3,3,9);plot(counter,whist,'-o');ylabel('w');

%M(count) = getframe(gcf);
%movie2avi(M,'run1.avi','compression','none','fps',30);

pause(0.1);

