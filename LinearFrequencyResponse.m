function [ Hha, Hhta, ph, V, D, fdom ] = LinearFrequencyResponse( M, K, C, f, r, k_t)
%linear frequency response for a MDOF of masses in series
%fdom=Frequency domain.
dofs=length(diag(M));

%% EIGEN PROBLEM
[V,D] = eig(K,M);%V:eigen vectors, D:eigen values.
d=diag(D);
W = sqrt(d);

fdom=0.001:0.001:1.2*(max(W));%Frequency domain

%Transfer function Multi forcing
a=f*V;


for i=1:dofs    
    NUM(i)=V(r,i)*a(i);%Normalized eigenvectors V
end

%Damping model
%C=bK+gM


c = C(1, 1)/2;
g = 0;


k = (K(1, 1)-k_t)/2;
b = c/k;
nr= b+g./d;%damping loss factor d IS NOT the vector of damped natural frequencies.

for j=1:dofs
    for i=1:length(fdom)
        denhd=complex((d(j)-fdom(i)^2), nr(j)*d(j));%denominator.
        Hh(i,j)=NUM(j)/denhd;%Mode contribution SPD:Simple Prop. Damp.
        
    end
end
Hha=abs(Hh);%Absolute part

for i=1:length(fdom)
    Hhtc(i)=sum((Hh(i,:)));%Complex total transfer function.
end

for i=1:length(fdom)
    Hhta(i)=abs(Hhtc(i));%Absolute part of the total transfer function.
end

for i=1:length(fdom)
    ph(i)=angle(Hhtc(i));%Phase part of the total complex transfer function.
end

%ph=unwrap(ph);
ph=ph*180/pi;

end
%% plot
% 
% figure(1)
% for i=1:dofs
%     
%     txt = ['Mode_' , num2str(i)];
%     semilogy(fdom,abs(Hha(:,i)), 'DisplayName', txt)
%     %[maxv(i), maxp(i)]=max(x(:,i));%Value and position of the maximium of every mode contribution
%     hold on
% end
% grid on
% 
% 
% semilogy(fdom,abs(Hhta))
% legend show
%legend('Total Transfer function')
