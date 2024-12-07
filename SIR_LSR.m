function LSRoutput = SIR_LSR(trueInfect,N,S0,gamma,dayNum)
    
    %Assumes trueInfect is a 2xM matrix with 1st row being days, 2nd row being infected population
    if dayNum > size(trueInfect(1,2:end),2)
        error('Number of days to use is greater than number of possible days in matrix');
    end

    spacing = floor(size(trueInfect(1,2:end),2)/dayNum);
    remainder = mod(size(trueInfect(1,2:end),2),dayNum);
    
    %Code to check if the number of days to use is a factor of the total number of days
    %Spacing can only be an integer, hence why this is needed.
    %
    %if it is not the case, last day used is the last day in the set minus the remainder.
    %
    %otherwise, uses last day with the spacing expecteds
    if (remainder ~= 0)
        %excludes day 0 and any day not divisible by spacing. Last day in set unused
        days = trueInfect(1, 1+spacing:spacing:dayNum+1-remainder); 
        infectList = trueInfect(2,1+spacing:spacing:dayNum+1-remainder);
    else
        %excludes day 0 and any day not divisible by spacing. Last day used
        days = trueInfect(1, 1+spacing:spacing:dayNum+1); 
        infectList = log(trueInfect(2,1+spacing:spacing:dayNum+1));
    end

    
    

    daySum = sum(days);
    infectSum = sum(infectList);
    dayInfectSum = sum(days .* infectList);
    daySqrSum = sum(days.^2);

    %ln I(t) = ln(I0) + kt -> y = a0 + a1 * t
    
    m = size(days,2);
    a1 = (m * dayInfectSum - daySum * infectSum) / (m * daySqrSum - (daySum)^2);
    a0 = infectSum / m - a1 * daySum / m;
    
    I0 = exp(a0);
    %k = (beta * S0 /N - gamma) -> beta = (k + gamma) * N/S0; k = a1
    beta = (a1 + gamma) * N/S0;

    LSRoutput = [I0,beta];

end