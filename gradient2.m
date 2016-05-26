function [dfx,dfy]=gradient2(f,x,y)
%GRADIENT2 Gradient Based on Quadratic Fitting.
% GRADIENT2(F,X) returns the numerical gradient dF/dX for vectors F and X.
% X need NOT be EQUALLY spaced. If X is not given, X=1:length(F) is
% assumed.
%
% GRADIENT2(F,X) when F is a matrix returns the numerical gradient dF/dX
% over each column of F. X need NOT be equally spaced. The number of
% elements in X must equal the number of rows in F, i.e., size(F,1). If X
% is not given, X=1:size(F,1) is assumed.
%
% [dFx,dFy]=GRADIENT2(F,X,Y) when F is a matrix returns the 2D numerical
% gradients dFx = dF/dx and dFy = dF/dy, where vector X contains the x-axis
% data points and vector Y contains the y-axis data points. X and Y need
% NOT be EQUALLY spaced. The number of elements in X must equal the number
% of columns in F. The number of elements in Y must equal the number of
% rows in F. If X and Y are not given, X=1:size(F,2) and Y=1:size(F,1)
% are assumed.
% [dFx,dFy]=GRADIENT2(F,X,Y) when X and Y are matrices the same size as F
% assumes that X and Y are 2D plaid as produced by [X,Y]=MESHGRID(x,y).
%
% [dFx,dFy]=GRADIENT2(F,Dx,Dy) for scalar Dx and Dy, assumes that Dx and Dy
% are the X-axis and Y-axis data spacings respectively.
%
% Algorithm: To support non equally spaced X and MATLAB vectorization, a
% quadratic polynomial fitting approach is taken. A second order polynomial
% is fit to each sequence of three consecutive data points, e.g., X(i-1),
% X(i), and X(i+1) for i=2:length(X)-1. The slope of this polynomial is
% then computed at X(i-1), X(i), and X(i+1). Using this approach, one slope
% is computed at X(1) and X(end), giving a second order gradient on the
% edges. Two polynomial slopes are computed at X(2) and X(end). Three
% polynomial slopes are computed at all other interior points. These slopes
% are optimally weighted for equally spaced data to give the gradient at
% interior points. For unequally spaced data, the weighting is close to
% optimal. This doesn't fit any textbook algorithm, but is easily
% implemented in MATLAB. Accuracy is typically an order of magnitude
% greater than the function GRADIENT.
%
% See also GRADIENT, DIFF, DEL2.

% D.C. Hanselman, University of Maine, Orono, ME 04469
% MasteringMatlab@yahoo.com
% Mastering MATLAB 7
% 2005-11-14

if nargin==0
   error('At Least One Input is Required.')
else
   fisrow=false;
   [rf,cf]=size(f);
   if max(rf,cf)<4
      error('F Must Contain at Least 4 Data Points.')
   elseif rf==1 % row vector vector F
      fisrow=true;
      f=f.';
      rf=cf;
      cf=1;
   end
end
if nargin==1                       % GRADIENT2(F) or [dFx,dFy]=GRADIENT2(F)
   x=(1:cf)';
   y=(1:rf)';
elseif nargin==2               % GRADIENT2(F,X) or [dFx,dFy]=GRADIENT2(F,X)
   x=x(:);
   if length(x)~=rf
      error('Length of X Must Match Length of F or Rows of Matrix F.')
   end
   y=(1:cf)';
   if nargout==2
      error('GRADIENT2(F,X) Requires a Single Output Variable')
   end
elseif min(rf,cf)>1               % Matrix F and [dFx,dFy]=GRADIENT2(F,X,Y)
   if all(size(x)==size(f)) % plaid X
      x=x(1,:)'; % get first row
   elseif isscalar(x) % GRADIENT(F,Dx,Dy)
      x=(0:cf-1)'*x;
   else % vector x
      x=x(:);
   end
   if all(size(y)==size(f)) % plaid Y
      y=y(:,1); % get first column
   elseif isscalar(y) % GRADIENT(F,Dx,Dy)
      y=(0:rf-1)'*y;
   else
      y=y(:);
   end
   if length(x)~=cf
      error('Length of X Must Match Number of Columns in Matrix F.')
   elseif length(y)~=rf
      error('Length of Y Must Match Number of Rows in Matrix F.')
   end
end
% now have column x, column y, and column or matrix f

if nargout<2            % 1D gradient (potentially across multiple columns)  
   dfx=local_getgrad(f,x);
   if fisrow
      dfx=dfx.';
   end
else                                                          % 2D gradient
   dfx=local_getgrad(f.',x).';
   dfy=local_getgrad(f,y);
end
%------------------------------------------------------------
function df=local_getgrad(f,z)
% core of gradient algorithm
% f = a*z^2 + b*z + c
% df = 2*a*z + b

dz=repmat(diff(z),1,size(f,2)); % [z(2)-z(1) z(3)-z(2) ... z(N)-z(N-1)]
df=diff(f);                     % [f(2)-f(1) f(3)-f(2) ... f(N)-f(N-1)]

dz1=dz(1:end-1,:);% [z(2)-z(1) z(3)-z(2) ... z(N-1)-z(N-2)], N-2 pts
dz2=dz(2:end,:);  % [z(3)-z(2) z(4)-z(3) ... z(N)-z(N-1)],   N-2 pts
dzs=dz1+dz2;      % [z(3)-z(1) z(4)-z(2) ... z(N)-z(N-2)],   N-2 pts
df1=df(1:end-1,:)./(dz1.*dzs); % normalized df1, N-2 pts
df2=df(2:end,:)./(dz2.*dzs);   % normalized df2, N-2 pts

a=2*(df2 - df1);               % a coef., N-2 pts
b=df2.*dz1 + df1.*dz2;         % b coef., N-2 pts

df1=-a.*dz1 + b; % slope at left polynomial point
df2=a.*dz2 + b;  % slope at right polynomial point
df0=2*b;         % twice slope at middle polynomial point

% approximate gradient using optimum weightings of polynomial slopes

df=[df1(1,:)                                           % left slope at z(1)
    (df0(1,:) + df1(2,:))/3             % 2/3 weight on center, 1/3 on left
    (df0(2:end-1,:) + (df1(3:end,:)+df2(1:end-2,:))/2)/3      % middle data
    (df0(end,:) + df2(end-1,:))/3      % 2/3 weight on center, 1/3 on right
    df2(end,:)];                                    % right slope at z(end)

 
% use the following to get standard central differences in the middle
% and quadratic slope on the edges
% this is typically an order of magnitude less accurate than the above

% df=[df1(1,:)                        % left slope at z(1)
%     df0(1:end,:)/2                  % middles on middles
%     df2(end,:)];                    % right slope at z(end)
