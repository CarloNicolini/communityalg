function sum = sumFactorial(n)
sum = 0;
if(n > 1)
    for i = 2:n
        sum = sum + log(i);
    end
end
