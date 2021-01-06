function [C] = make_C(c1, c2, dim)
%Gives the DAMPING Matrix
if dim ==3
    
C = diag(c2*ones(1, dim)) + diag([c1 c2 c1]);
C = C - diag(c2*ones(1,dim-1),1) - diag(c2*ones(1,dim-1),-1);
end
if dim ==2
C = diag([c1 c2+c1]);
C(dim-1, dim) = -c1;
C(dim, dim-1) = -c1;
end
end