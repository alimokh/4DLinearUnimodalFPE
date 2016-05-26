%Copyright 2013 Rice University
%All Rights Reserved
%
%Created by Ali Arya Mokhtarzadeh
%Department of Mechanical Engineering
%

function [LR, o] = FPE4DSS6OptSymD2(par,p,d,v)

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

%x1c = ucbc.x1(count);
%%%%%%%%%      Numerical Derivatives - No Boundary Enforcers!!! %%%%%
%syms x
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

%eps1 = dot(o.funpsi(2:d.firstlength-1,2:d.firstlength-1),p.e(:));eps2 = dot(o.funpsi(2:d.firstlength-1,2:d.firstlength-1),o.funpsi(2:d.firstlength-1,2:d.firstlength-1));eps3 = dot(p.e(:),p.e(:));
%LR = 1-eps1^2/(eps2*eps3);
%c = -eps1/eps2;

newoptbasis = o.c*o.funpsi;%%%%%%%%%%%Turn ON for RBF
%newoptbasis = o.funpsi;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Only for Piecewise Basis Functions because o.c is already used to determine their height

newoptbasis = newoptbasis(2:d.firstlength-1,2:d.firstlength-1);
p.e = p.e + newoptbasis;


%LR = sum((p.e(:).^2));

%p.testerr = p.e(12:d.firstlength-11,12:d.firstlength-11);
%LR = norm(p.testerr(:));
%LR = max(abs(p.e(:)));%LR = trapz(d.x1,p.e.^2);

%p.pe
%p.e

%LR = sqrt( mse(p.pe - p.p) );

%o.funphi = max(0,o.funphi); %Matt's idea for positivity
%Following Lines: AAM 8/19/13 Possible Idea for Positivity Enforcement
%add basis func to approx
%find locations where value of approx is zero or negative
%set basis function value to zero at those locations
%update the approx with the new basis function


newbasis = o.c*o.funphi;%%%%%%%%%%%Turn ON for RBF
%newbasis = o.funphi;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Only for Piecewise Basis Functions because o.c is already used to determine their height

newbasis = newbasis(2:d.firstlength-1,2:d.firstlength-1);
penalty = p.p + newbasis; 
if ((any(penalty(:) < 0)) ||  o.x1c > d.upx1+1 || o.x1c < d.lowx1-1 || o.x2c > d.upx2+1 || o.x2c < d.lowx2-1 ) % || o.w <= 0 || o.w >= 1 || o.c < 0.001
    %the one about o.w is questionable AAM 6/9/13     || o.w >= 1  
    %but the o.w term seems to prevent the error of the basis function having complex
    %terms AAM 6/12/13      ??? Error using ==> mesh at 80     X, Y, Z, and C cannot be complex.
    LR = 1e50;
    
    % || penalty(d.mid-3,d.mid-3) ~= 1 %try use this to fix one point in
    % domain AAM 6/28/13
else
    LR = norm(p.e(:));
    
end
 
%newbasis = c*funphi;
%newbasis = newbasis(3:d.firstlength-2);
%p.p = p.p + newbasis;
%newres = c*funpsi;
%newres = newres(3:d.firstlength-2);
%p.e = p.e + newres;
%LR = norm(p.e);% differential equation RMS

if (isnan(LR) || isinf(LR) || ~isreal(LR) )
    LR = 1e50;
end

%dimension = length(ind);
%if order == 1
%    for thiscount = 1:dimension
%        temp = (diff(o.funphi,ind(thiscount)));
%        fd(thiscount) = temp;
%    end
%end