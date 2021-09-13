function [ output] = QAM_dec( data )
for k = 1: length(data)
    if real(data(k)) <= 0 && real(data(k)) >= -2 &&  imag(data(k)) <= 0 && imag(data(k)) >= -2
        output(k,:) = [0   1   0    1];
    elseif  real(data(k)) <= 0 && real(data(k)) <= -2 &&  imag(data(k)) <= 0 && imag(data(k)) >= -2
        output(k,:) = [0    0    0   1];
    elseif  real(data(k)) >= 0 && real(data(k)) >= 2 &&  imag(data(k)) <= 0 && imag(data(k)) >= -2
        output(k,:) = [1    0    0   1];
    elseif  real(data(k)) >= 0 && real(data(k)) <= 2 &&  imag(data(k)) <= 0 && imag(data(k)) >= -2
        output(k,:) = [1    1    0   1];
    elseif  real(data(k)) <= 0 && real(data(k)) >= -2 &&  imag(data(k)) >= 0 && imag(data(k)) <= 2
        output(k,:) =[0    1    1   1];
    elseif  real(data(k)) <= 0 && real(data(k)) <= -2 &&  imag(data(k)) >= 0 && imag(data(k)) <= 2
        output(k,:) = [0    0    1   1];
    elseif  real(data(k)) >= 0 && real(data(k)) >= 2 &&  imag(data(k)) >= 0 && imag(data(k)) <= 2
        output(k,:) = [1    0    1   1];
    elseif  real(data(k)) >= 0 && real(data(k)) <= 2 &&  imag(data(k)) >= 0 && imag(data(k)) <= 2
        output(k,:) = [1    1    1   1];
    elseif  real(data(k)) <= 0 && real(data(k)) >= -2 &&  imag(data(k)) <= 0 && imag(data(k)) <= -2
        output(k,:) = [0    1    0   0];
    elseif  real(data(k)) <= 0 && real(data(k)) <= -2 &&  imag(data(k)) <= 0 && imag(data(k)) <= -2
        output(k,:) = [0    0    0   0];
    elseif  real(data(k)) >= 0 && real(data(k)) >= 2 &&  imag(data(k)) <= 0 && imag(data(k)) <= -2
        output(k,:) = [1    0    0   0];
    elseif  real(data(k)) >= 0 && real(data(k)) <= 2 &&  imag(data(k)) <= 0 && imag(data(k)) <= -2
        output(k,:) = [1    1    0   0];
    elseif  real(data(k)) <= 0 && real(data(k)) >= -2 &&  imag(data(k)) >= 0 && imag(data(k)) >= 2
        output(k,:) = -[0    1    1   0];
    elseif  real(data(k)) <= 0 && real(data(k)) <= -2 &&  imag(data(k)) >= 0 && imag(data(k)) >= 2
        output(k,:) = [0    0    1   0];
    elseif  real(data(k)) >= 0 && real(data(k)) >= 2 &&  imag(data(k)) >= 0 && imag(data(k)) >= 2
        output(k,:) = [1    0    1  0];
    elseif  real(data(k)) >= 0 && real(data(k)) <= 2 &&  imag(data(k)) >= 0 && imag(data(k)) >= 2
        output(k,:) = [1    1    1   0];
    end
end


end

