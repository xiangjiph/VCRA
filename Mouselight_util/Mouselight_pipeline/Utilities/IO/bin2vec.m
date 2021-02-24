function vec = bin2vec(num, dig)
num = num-1;
vec = zeros(1, dig);
for i = 1:dig
    fac = 2.^(dig-i);
    vec(i) = floor(num./fac);
    num = num-vec(i)*fac;
end


