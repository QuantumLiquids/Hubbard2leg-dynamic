function n_k=MyFourierTrans(k, x, n_x)
% n_k= sum_x n_x exp(1i*kx)
%x, n_x are row vectors
%k is a number or a arbirary length row vector(do not need to equal the
%length of x)
if size(x)~=size(n_x)
    error('the size of x and n_x is different');
end
if size(x, 1) ~= 1
    error('x is not a row vector');
end
n_k = sum(n_x.' .* exp(1i * k .* x.')); %every column is a k, sum is follow the column direction
end