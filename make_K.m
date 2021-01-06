function [K] = make_K(k1, k2, dim)
%Gives the STIFNESS Matrix
if dim ==3
    
K = diag(k2*ones(1, dim)) + diag([k1 k2 k1]);
K = K - diag(k2*ones(1,dim-1),1) - diag(k2*ones(1,dim-1),-1);
end
if dim ==2
K = diag([k1 k2+k1]);
K(dim-1, dim) = -k1;
K(dim, dim-1) = -k1;
    
end
end

