function y = fw(f, s, t, T)
sz = size(s);
w1 = (1-s); % Bug fix
w2 = s.*t;
w3 = 1-w1-w2;
P = [w1(:),w2(:),w3(:)] * T;
y = feval(f,P(:,1),P(:,2));
y = s(:).*y(:);
y = reshape(y,sz);
end