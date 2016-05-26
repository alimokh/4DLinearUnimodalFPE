%Copyright 2013 Rice University
%All Rights Reserved
%
%Created by Ali Arya Mokhtarzadeh
%Department of Mechanical Engineering
%

function [radial] = rbfunc(chist, whist, x1hist, x2hist, newd, newcount)

newc = chist(newcount);
neww = whist(newcount);
newx1c = x1hist(newcount);
newx2c = x2hist(newcount);

radial = newc.*exp(neww*(((newd.x1-newx1c).^2)+(newd.x2-newx2c).^2));


end
