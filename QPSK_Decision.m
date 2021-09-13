function bk_hat = QPSK_Decision(yk,h)
yk = yk ./ h;
no_rand_bits = length(yk);
bk_hat = ones(1,no_rand_bits);
for v = 1:no_rand_bits
    if real(yk(v)) <= 0 && imag(yk(v)) <= 0
        bk_hat(v) = -1 - 1i;
    elseif  real(yk(v)) <= 0 && imag(yk(v)) >= 0
        bk_hat(v) = -1 + 1i;
    elseif  real(yk(v)) >= 0 && imag(yk(v)) <= 0
        bk_hat(v) = 1 - 1i;
    elseif  real(yk(v)) >= 0 && imag(yk(v)) >= 0
        bk_hat(v) = 1 + 1i;
    end
end
end