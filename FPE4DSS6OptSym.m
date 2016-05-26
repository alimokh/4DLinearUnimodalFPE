%Copyright 2013 Rice University
%All Rights Reserved
%
%Created by Ali Arya Mokhtarzadeh
%Department of Mechanical Engineering
%

function [LR, o] = FPE1DSSnd6OptSym(par,p,d,v)

o.w = -100*par(1)^2;
o.x1c = (d.upx1-d.lowx1)*par(2)+d.lowx1;

%x1c = ucbc.x1(count);
%%%%%%%%%      Numerical Derivatives - No Boundary Enforcers!!! %%%%%
syms x

casenumber = 8;
o = TestCopyD1Input(casenumber,v,d,x,p,o);
%funphi = (w^((x-x1c)^2));


o.phix1 = (diff(o.funphi,x));
o.phix1x1 = diff(o.phix1,x);

x = d.x1;

o.funphi = subs(o.funphi);
o.phix1 = subs(o.phix1);
o.phix1x1 = subs(o.phix1x1);
  

casenumber = 9;
o = TestCopyD1Input(casenumber,v,d,x,p,o);

eps1 = dot(o.funpsi(3:d.firstlength-2),p.e(:));eps2 = dot(o.funpsi(3:d.firstlength-2),o.funpsi(3:d.firstlength-2));eps3 = dot(p.e(:),p.e(:));
LR = 1-eps1^2/(eps2*eps3);
c = -eps1/eps2;
 
%newbasis = c*funphi;
%newbasis = newbasis(3:d.firstlength-2);
%p.p = p.p + newbasis;
%newres = c*funpsi;
%newres = newres(3:d.firstlength-2);
%p.e = p.e + newres;
%LR = norm(p.e);% differential equation RMS

if (isnan(LR) || isinf(LR) || ~isreal(LR))
    LR = 2;
end

%dimension = length(ind);
%if order == 1
%    for thiscount = 1:dimension
%        temp = (diff(o.funphi,ind(thiscount)));
%        fd(thiscount) = temp;
%    end
%end