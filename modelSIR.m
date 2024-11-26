function SIRmatrix = modelSIR(S0, I0, R0, T0, T, h, beta, gamma)
    
    N = S0 + I0 + R0;
    days = T0:h:T;
    
    %Preallocates matrix size based on size of "days".
    susepList = zeros(1,size(days,2)); 
    infectList = zeros(1,size(days,2));
    recoverList = zeros(1,size(days,2));
    
    %Sets values for initial date.
    susepList(1) = S0;
    infectList(1) = I0;
    recoverList(1) = R0;
    
    %Models with respect to time
    dSdt = @(S,I) -beta/N * S * I;
    dIdt = @(S,I) beta/N * S * I - gamma * I;
    dRdt = @(I) gamma*I;
    
    K = zeros(3,4); %Rows correspond to SIR, Columns correspond to Runge Kutta constants K1-K4

    for i = 1:size(days,2)-1
        K(1,1) = dSdt(susepList(i),infectList(i));
        K(2,1) = dIdt(susepList(i),infectList(i));
        K(3,1) = dRdt(infectList(i));

        K(1,2) = dSdt(susepList(i) +K(1,1) * h/2, infectList(i) + K(2,1) * h/2);
        K(2,2) = dIdt(susepList(i) +K(1,1) * h/2, infectList(i) + K(2,1) * h/2);
        K(3,2) = dRdt(infectList(i) + K(2,1) * h/2);

        K(1,3) = dSdt(susepList(i) +K(1,2) * h/2, infectList(i) + K(2,2) * h/2);
        K(2,3) = dIdt(susepList(i) +K(1,2) * h/2, infectList(i) + K(2,2) * h/2);
        K(3,3) = dRdt(infectList(i) + K(2,2) * h/2);

        K(1,4) = dSdt(susepList(i) +K(1,3) * h, infectList(i) + K(2,3) * h);
        K(2,4) = dIdt(susepList(i) +K(1,3) * h, infectList(i) + K(2,3) * h);
        K(3,4) = dRdt(infectList(i) + K(2,3) * h);

        susepList(i+1) = susepList(i) + 1/6 * (K(1,1) + 2*K(1,2) + 2*K(1,3) + K(1,4))*h;
        infectList(i+1) = infectList(i) + 1/6 * (K(2,1) + 2*K(2,2) + 2*K(2,3) + K(2,4))*h;
        recoverList(i+1) = recoverList(i) + 1/6 * (K(3,1) + 2*K(3,2) + 2*K(3,3) + K(3,4))*h;
    end

    %Returns a 4 row matrix.
    %Row 1: days
    %Row 2: Susceptible population
    %Row 3: Infected population
    %Row 4: Recovered population

    SIRmatrix = [days;susepList;infectList;recoverList];
end