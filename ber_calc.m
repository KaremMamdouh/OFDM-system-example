

function bit_error = ber_calc(xk,bk_hat,rep)
bit_error = 0;
count = 0;
no_rand_bits = length(xk);
for v = 1:rep:no_rand_bits
    for s = 0 : rep-1
        if xk(v+s) ~= bk_hat(v+s)
            count = count + 1;
        end
    end
    if count >= rep / 2
        bit_error = bit_error + 1;
    end
    count = 0;
end
bit_error = bit_error / no_rand_bits;
end

