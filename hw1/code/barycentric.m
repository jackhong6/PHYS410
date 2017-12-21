function values = barycentric(VecX,VecY,x)
% Barycentric Lagrange interpolation.
% VecX, VecY: interpolation data. x: observation points.
% PHYS 410: Homework 1
% Modified from given by Jack Hong, 30935134

%input checking
n=length(VecX);
if length(VecY) ~=n; disp('Interpolation data inconsistent'); return;end

values = zeros(1,length(x));
iter = 1;
% Calculate the interpolant value for each point in x
for x0 = x
    
    % calculating the weights
    V=bsxfun(@minus,VecX,VecX.');
    V(1:n+1:end) = 1;
    VV = exp(sum(log(abs(V))));
    w = 1./(prod(sign(V)).*VV);
    
    % calculating the interpolant
    num=sum(w.*VecY./(x0-VecX));
    den=sum(w./(x0-VecX));
    
    values(iter) = num/den;
    
    iter = iter + 1;
end

end