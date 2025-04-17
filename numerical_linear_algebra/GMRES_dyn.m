function [x,k] = GMRES_dyn(A,b,epsilon)

[n,~] = size(A);
v = b/norm(b);
V = v;
H = [];
u = [norm(b);0];
k = 1;
r = 10;

while r/norm(b) > epsilon
    w = A*v;
    h = V'*w;
    % (1) Modified Gram-Schmidt
    for p=1:k
        w = (eye(n) - V(:,p)*V(:,p)') * w;
    end
    g = norm(w);
    v = w/g;
    V = [V v];
    H = [H,h; zeros(1,k-1),g];
    
    % (2,3,4) Dynamic updating
    G = planerot([H(k,k); H(k+1,k)]);
    H(k:k+1,k) = G*H(k:k+1,k);
    u(k:k+1) = G*u(k:k+1);
    
    x = V(:,1:k)*(H(1:k,1:k)\u(1:k));
    r = abs(u(k+1));
    
    k = k+1;
    u = [u;0];
end

