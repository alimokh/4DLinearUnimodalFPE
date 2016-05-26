%Copyright 2013 Rice University
%All Rights Reserved
%
%Created by Ali Arya Mokhtarzadeh
%Department of Mechanical Engineering
%

newdisc = 0.05;
[x1 x2] = ndgrid(d.lowx1:newdisc:d.upx1, d.lowx2:newdisc:d.upx2);
newd.x1 = x1;
newd.x2 = x2;

initapprox = (1/(sqrt(1/9)*sqrt(2*pi)))*exp(-(9/2)*(((newd.x1-1).^2)+(newd.x2-1).^2));

newapproximation = initapprox;

newcount = 0;

while(count<l.Nmax)
    newcount = newcount+1;
    radial = rbfunc(chist, whist, x1hist, x2hist, newd, newcount);
    newapproximation = newapproximation + radial;
end

mesh(newd.x1,newd.x2,newapproximation);title('New Approximation');