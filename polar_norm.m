function [ norm ] = polar_norm(v,r)
%snake_norm Compute the norm of column vector
%   with the snake metric tensor given alpha and beta
    dims = size(v);
    assert( dims(2) == 1);
    switch dims(1)
        case 2 %induced metric on spacelike 2D hyperplane
             norm = v(1)^2 + r^2*v(2)^2 ;
        case 3 %induced cylindrical metric on spacelike 3D hyperplane
             norm = v(1)^2 + r^2*v(2)^2 + v(3)^2 ;
        otherwise
            error('Vector has incorrect dimensions!');
    end
end