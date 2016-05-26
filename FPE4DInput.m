%Copyright 2013 Rice University
%All Rights Reserved
%
%Created by Ali Arya Mokhtarzadeh
%Department of Mechanical Engineering
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%          SOMA INPUT FILE               %%% %TestCopyD1Input.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output] = FPE4DInput(casenumber,v,d,x,p,o)

if casenumber == 1
  %Enter the equation variables and specify the values                       %v.field
  output.k1 = 1; output.k2 = 1; output.k3 = 1; output.c1 = 0.4; output.c2 = 0.4; output.D1 = 0.2; output.D2 = 0.2;
end
if casenumber == 2
  %Enter the size of your domain                                             %d.field
  output.lowx1  = -4;output.upx1   =  4;output.upx2   = output.upx1;output.lowx2   = output.lowx1;
  output.upx3   = output.upx1;output.lowx3   = output.lowx1;output.upx4   = output.upx1;output.lowx4   = output.lowx1;
  output.disc   =  0.1; output.numpar =  6;    %Number of optimization parameters (1 + D + (1 if nonlinear))
end
if casenumber == 3
  %Enter the True Solution
  Dtrusol = 4;
  A = [0 1 0 0; -(v.k1 + v.k2) -v.c2 v.k2 0; 0 0 0 1; v.k2 0 -(v.k2 + v.k3) -v.c2];
  B = zeros(Dtrusol,Dtrusol);
  G = [0 0; 1 0; 0 0; 0 1];
  Q = [2*v.D1 0; 0 2*v.D2];
  [Pstat, Lstat, Gstat] = care(A, B, G*Q*G');        % Pstat is the stationary covariance matrix
  Mstat = zeros(Dtrusol,1);        % Mstat is the stationary mean. Clearly this is zero because we have a stable linear system
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Sigma = Pstat(1:2,1:2);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  mu = [0 0];
  %Sigma = [.25 .3; .3 1];
  x1a = -4:d.disc:4; x2a = -4:d.disc:4;
  [X1A,X2A] = meshgrid(x1a,x2a);
  F = mvnpdf([X1A(:) X2A(:)],mu,Sigma);
  output = reshape(F,length(x2a),length(x1a));
end
if casenumber == 4                                                                                                                                              %l.field
  output.Nmax   =  30;   %Maximum number of basis functions
  output.tol    =  1e-9; %Maximum allowable error in the Galerkin Integral
end
if casenumber == 5
  %Enter the independent variables (t, x, y, x1, x2, etc)
  %syms x
  x1 = d.x1; x2 = d.x2;x3 = d.x3; x4 = d.x4;
  output.x1 = x1; output.x2 = x2;output.x3 = x3; output.x4 = x4;
end
if casenumber == 6
  %Enter the Initial Conditions
  output.p = (1/sqrt(2*pi))*exp(-((d.x1.^2)+(d.x2.^2)+(d.x3.^2)+(d.x4.^2)));
end
if casenumber == 7
  %Enter the Equation to be solved (in such a form that the left side is 0)
  output.e = v.D1*p.px2x2 + v.D2*p.px4x4 - d.x2.*p.px1 - d.x4.*p.px3 + (v.k1+v.k2)*d.x1.*p.px2 + v.c1*d.x2.*p.px2 - v.k2*d.x3.*p.px2 - v.k2*d.x1.*p.px4 + (v.k2+v.k3)*d.x3.*p.px4 + v.c2*d.x2.*p.px4 + (v.c1+v.c2)*p.p;
  output.p = p.p;
end
if casenumber == 8
  output.funphi = exp(o.w*(((d.x1-o.x1c).^2)+((d.x2-o.x2c).^2)+((d.x3-o.x3c).^2)+((d.x4-o.x4c).^2)));
  %output.funphi = o.w.^((((d.x1-o.x1c).^2)+(d.x2-o.x2c).^2));
  output.w = o.w;output.x1c = o.x1c;output.x2c = o.x2c;output.x3c = o.x3c;output.x4c = o.x4c;output.c = o.c;
end
if casenumber == 9
  %Enter the Equation again, whith the "p's" replaced by "phi"
  output.funpsi = v.D1*o.phix2x2 + v.D2*o.phix4x4 - d.x2.*o.phix1 - d.x4.*o.phix3 + (v.k1+v.k2)*d.x1.*o.phix2 + v.c1*d.x2.*o.phix2 - v.k2*d.x3.*o.phix2 - v.k2*d.x1.*o.phix4 + (v.k2+v.k3)*d.x3.*o.phix4 + v.c2*d.x2.*o.phix4 + (v.c1+v.c2)*o.funphi;
  output.funphi = o.funphi;
  output.w = o.w;output.x1c = o.x1c;output.x2c = o.x2c;output.x3c = o.x3c;output.x4c = o.x4c;output.c = o.c;
end

