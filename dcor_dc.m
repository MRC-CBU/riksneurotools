function [R] = dcor_dc(X,Y)
% [R] = dcor_dc(X,Y)
% Computes the double centered distance correlation between X and Y. 
% Rows represent the examples, and columns the variables.

% Based on: http://www.mathworks.com/matlabcentral/fileexchange/49968-dcorr--x--y--
% and, the R package energy and papers by Szekely and Rizzo (2007; 2013 and 2014). 

% Author: Linda Geerligs (lindageerligs@gmail.com), Date: 22-10-2015

a = pdist2(X, X);
b = pdist2(Y, Y);
n=size(X,1);

A = Dcenter(a);
B = Dcenter(b);

dcovXY = sum(sum(A.*B)) ./ (n.^2);
dvarX = sum(sum(A.*A)) ./ (n.^2);
dvarY = sum(sum(B.*B)) ./ (n.^2);

R=sqrt(dcovXY / sqrt(dvarX * dvarY));

    function A=Dcenter(a) 
        A = a - bsxfun(@plus,mean(a),mean(a,2))+mean(a(:));
    end

end