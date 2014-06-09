%This function finds the tform matrix and params vector for a
%set of 4 input and 4 output points only

function [params tform] = GetProjectiveTransformMatrix(x, y, u, v)
    
    M = zeros(8,8);
    targetVector = zeros(8,1);
    
    for n = 1:4
        M(2*n - 1,:) = [x(n) y(n) 1 0 0 0 -x(n)*u(n) -y(n)*u(n)];
        M(2*n    ,:) = [0 0 0 x(n) y(n) 1 -x(n)*v(n) -y(n)*v(n)]';
        
        targetVector(2*n - 1) = u(n);
        targetVector(2*n    ) = v(n);
    end
    
    params = M\targetVector;
    
    [a b c d e f g h] = deal(params(1), ...
    params(2), ...
    params(3), ...
    params(4), ...
    params(5), ...
    params(6), ...
    params(7), ...
    params(8));

    tform = [a, d, g; ...
             b, e, h; ...
             c, f, 1];
end