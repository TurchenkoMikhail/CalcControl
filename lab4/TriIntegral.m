function I = TriIntegral(f, Tx, Ty)
% I = TriIntegral(f, Tx, Ty)
% 2D integration of f on a triangle
% INPUTS:
%   - f is the vectorized function handle that when calling f(x,y) returns
%       function value at (x,y), x and y are column vectors
%   - Tx,Ty are is two vectors of length 3, coordinates of the triangle
% OUTPUT
%   I: integral of f in T

T = [Tx(:), Ty(:)];
I = integral2(@(s,t) fw(f,s,t,T),0,1,0,1);
A = det(T(2:3,:)-T(1,:));
I = I*abs(A);
end % TriIntegral
