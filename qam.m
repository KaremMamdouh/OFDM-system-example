function [ output ] = qam( data )
%UNTkTLED5 Summary of thks functkon goes here
%   Detakled explanatkon goes here
%******************** QAM Modulatkon *********************** 

M=16;
% table(:) = [-3-3k, -3-1k, -3+3k, -3+1k,-1-3k, -1-1k, -1+3k, -1+1k , 3-3k,3-1k,  3+3k,  3+1k ,  1-3k,  1-1k,  1+3k,  1+1k];


    %**********Grey encodkng************%

  
for k = 1 : length(data)
if data(k,:) == [0  0  0 0]
output(k) = -3+ -3i ;
elseif (data(k,:) == [0    0    0   1])
output(k) = -3- 1i ;
elseif(data(k,:) == [0    0    1   0])
output(k) = -3+ 3i;
elseif(data(k,:) == [0    0    1   1])
output(k) = -3 +1i ;
elseif(data(k,:) == [0    1    0   1])
    output(k) = -1- 1i ;

elseif(data(k,:) == [0    1    1   0])
    output(k) = -1+ 3i;

elseif(data(k,:) == [0    1    1   1])
    output(k) = -1+ 1i ;

elseif(data(k,:) == [1    0    0   0])
    output(k) = 3+ -3i ;

elseif(data(k,:) == [1    0    0   1])
    output(k) = 3- 1i ;
elseif(data(k,:) == [1    0    1   0])
    output(k) = 3+ 3i ;
elseif(data(k,:) == [0    1    0   0])
    output(k) = -1- 3i ;
elseif(data(k,:) == [1    1    0   0])
    output(k) = 1- 3i ;
elseif(data(k,:) == [1    1    1   0])
    output(k) = 1+ 3i ;
elseif(data(k,:) == [1    1    0   1])
    output(k) = 1- 1i;
elseif(data(k,:) == [1    0    1   1])
    output(k) = 3+ 1i ;
elseif(data(k,:) == [1    1    1   1])
    output(k) = 1+ 1i ;
end


end
   


end
end






