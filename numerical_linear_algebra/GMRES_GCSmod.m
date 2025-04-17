function [x,k] = GMRES_GCSmod(A,b,tol)
% %GMRES old
% v = b/norm(b); V = v; H = [];
% for j=1:k
%     w = A*v;
%     h = V'*w;
%     w = w-V*h;
%     g = norm(w);
%     H = [H,h; zeros(1,j-1),g];
%     v = w/g;
%     V = [V v];
% end

%% GMRES with MGS
v = b/norm(b); V = v;
H = [];
k = 0;
[m,n] = size(v);
%r_k = b - A*v
r=1

while r/norm(b) > tol
    k = k+1
    w = A*v;
    h = V'*w;
    % CGS version
    % w = w - V*h;
    % MGS version
    for i=1:k
        w = (eye(m) - V(:,i)*V(:,i)') * w;
    end
    g = norm(w);
    H = [H, h; zeros(1, k-1), g];
    v = w/g;
    V = [V v];
    %r_k = b - A*v;

    u = [norm(b);zeros(k,1)];

    x=[H(i,i);H(i+1,i)];
    [G,y] = planerot(x);
    H(i:i+1,:) = G*H(i:i+1,:);
    u(i:i+1,:) = G*u(i:i+1,:);
    x = V(:,1:k)*(H(1:k,1:k)\u(1:k));
    %r_k = b - A*x;
    r = abs(u(k+1))
end

V
H
u
%x=5
x = V(:,1:k)*(H(1:k,1:k)\u(1:k));
% ;
