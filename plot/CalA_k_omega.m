% ====== The defination of correlation mode ====== %
%
%   size_t CorrelationMode; //0 : zz <S_j^z(t) S_i^z(0)> with i fixed in
%   the middle of the chain
%                           //1 : pm <S_j^-(t) S_i^+(0)>
%                           //2 : mp <S_j^+(t) S_i^-(0)>


% delta_kx is a small interval space between continue kx values
% 
% output: spectral function at ky=0 and ky=pi
function [A0,Api] = CalA_k_omega(k_set, omega_set, time, x_set, G_t_x, gaussian_factor)
s = 1;
time_sign = 1;

% a/b for 2 legs
G_t_k_a = zeros(numel(time), numel(k_set));
G_t_k_b = zeros(numel(time), numel(k_set));

for i = 1: numel(time)
    G_t_k_a(i,:) = SymmetrizedFourierTrans(s * k_set, x_set, G_t_x(i,1:2:end)) * exp(-gaussian_factor * time(i)^2);
    G_t_k_b(i,:) = SymmetrizedFourierTrans(s * k_set, x_set, G_t_x(i,2:2:end)) * exp(-gaussian_factor * time(i)^2);
end

A_k_omega_a = zeros(numel(omega_set), numel(k_set));
A_k_omega_b = zeros(numel(omega_set), numel(k_set));
for i = 1: numel(k_set)
    G_t_on_the_k = G_t_k_a(:,i).';
    G_omega_on_the_k = MyFourierTrans(omega_set, time_sign*time, G_t_on_the_k);
    A_k_omega_a(:, i) = real(G_omega_on_the_k.');
end
for i = 1: numel(k_set)
    G_t_on_the_k = G_t_k_b(:,i).';
    G_omega_on_the_k = MyFourierTrans(omega_set, time_sign*time, G_t_on_the_k);
    A_k_omega_b(:, i) = real(G_omega_on_the_k.');
end
A0 = (A_k_omega_a+A_k_omega_b)/pi * time(2);
Api = (A_k_omega_a-A_k_omega_b)/pi * time(2);
end