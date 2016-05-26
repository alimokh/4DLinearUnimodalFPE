%Copyright 2013 Rice University
%All Rights Reserved
%
%Created by Ali Arya Mokhtarzadeh
%Department of Mechanical Engineering
%

function [initapprox] = initfunc(newd)

initapprox = (1/(sqrt(1/9)*sqrt(2*pi)))*exp(-(9/2)*(((newd.x1-1).^2)+(newd.x2-1).^2));


end