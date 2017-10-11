function [t,F,p,df,R2,cR2,B,r,aR2,iR2] = glm(y,X,c,pflag);

% [t,F,p,df,R2,cR2,B] = glm(y,X,c,pflag);
%
% Generic function for General Linear Model (GLM)
% Note: assumes error is spherical (residuals white)
% Rik Henson, 2004
%
% (Note requires Matlab stats toolbox, or SPM functions to be on path)
%
% Input:
%    y = n column vector of observations
%    X = n x p (design) matrix 
%    c = p x e contrast matrix for T- (e=1) or F- (e>1) contrast
%    pflag = whether to print fit
%
% Output
%    t = T-value
%    F = F-value
%    p = p-value
%    df = degrees of freedom 
%    R2  = model fit (% variance explained)
%    cR2 = contrast fit (% variance explained)
%    B = p betas (parameter estimates)

if nargin<4
	pflag=0;
end

B = pinv(X)*y;

Y = X*B;

r = y - Y;

if size(c,2) > size(c,1)
     warning('Transposing c!')
     c = c';
end

l = size(c,1);
if l < size(X,2)
     c = [c; zeros(size(X,2)-l,size(c,2))];
     warning('Padding c with zeros!')
end

df = length(y) - rank(X);

% T and F's could be combined, but done separately for pedagogical reasons

if size(c,2)==1
  s = r'*r / df;
  t = c'*B / sqrt(s*c'*pinv(X'*X)*c);
  p = t2p(t,df);
  F = t.^2; 
  cR2 = F2R(F,1,df);
  R2 = 1 - (r'*r) / (y'*y); 
  aR2 = 1 - (r'*r/df(end)) / (y'*y/(length(y)-1));
  fprintf('T(%d)=%f, p(one-tailed)=%f  (R2=%3.2f; overall R2=%3.2f (adjusted R2=%3.2f))\n',df,t,p,cR2,R2,aR2);
else
  c_0 = eye(size(X,2)) - c*pinv(c);
  X_0 = X*c_0;
  R   = eye(size(X,1)) - X*pinv(X);
  R_0 = eye(size(X,1)) - X_0*pinv(X_0);
  M = R_0 - R;
  df = [rank(X)-rank(X_0) size(X,1)-rank(X)];
  F  = ((B'*X'*M*X*B)/df(1)) / ((y'*R*y)/df(2)); 
  p  = F2p(F,df(1),df(2));
  t = [];
  cR2 = F2R(F,df(1),df(2));
  R2 = 1 - (r'*r) / (y'*y); 
  aR2 = 1 - (r'*r/df(end)) / (y'*y/(length(y)-1));
  fprintf('F(%d,%d)=%f, p(two-tailed)=%f  (R2=%3.2f; overall R2=%3.2f (adjusted R2=%3.2f))\n',df(1),df(2),F,p,cR2,R2,aR2);
end

%Below doesn't work if constant not included
%R2 = 1 - (r'*r) / (y'*y);   %R2 = (Y'*Y)/(y'*y) equivalent
%disp(sprintf('Overall R2 = %f',R2))
%aR2 = 1 - (r'*r/df(end)) / (y'*y/(length(y)-1));
%disp(sprintf('Overall Adjusted R2 (assumes a constant in model) = %f',R2))

% Hacky way to calculate unique contribution of contrast 
c_0 = eye(size(X,2)) - c*pinv(c);
X_0 = X*c_0;
y_0 = X_0*(pinv(X_0)*y);
r_0 = y - y_0;
iR2 = R2 - (1 - (r_0'*r_0) / (y'*y));


if pflag
    figure(pflag), clf
    Ne = size(c,2);
    for e = 1:Ne
        subplot(1,Ne,e), hold on
        Xc = X*c(:,e);
        Yc = Xc*(c(:,e)'*B);
        yc = Yc+r;
        plot(Xc,yc,'r.')
        plot(Xc,Yc,'b-')
    end
end

return

%%%%%%%%%%%%%%%%%%
function p=t2p(t,df);

% one-tailed p-value

try
    p=tcdf(t,df);
catch
    try 
        p=spm_Tcdf(t,df);
    catch
        error('Need Matlab Stats toolbox (tcdf) or SPM on path')
    end
end

f=find(p>0.5);
p(f)=1-p(f);

return

%%%%%%%%%%%%%%%%%%%%

function p=F2p(F,df1,df2);
% one-tailed p-value

try
    p=1-fcdf(F,df1,df2);
catch
    try 
        p=1-spm_Fcdf(F,df1,df2);
    catch
        error('Need Matlab Stats toolbox (fcdf) or SPM on path')
    end
end

return

%%%%%%%%%%%%%%%%%%%%
function R2=F2R(F,df1,df2);

R2=df1*F/(df2+df1*F);

return
