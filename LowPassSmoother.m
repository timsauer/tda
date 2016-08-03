function x = LowPassSmoother(Xi,L,dt,smootherwidth)
    
    Fs = 1/dt;
    f = Fs*(0:ceil(L/2))/L;
    x = zeros(size(f));
    
    for j = 1:length(Xi)
        x = x + exp(-2*pi*sqrt(-1)*Xi(j)*f);
    end
    
     x(f>1/smootherwidth)=0;
    % x(f>1/smootherwidth) = x(f>1/smootherwidth).*(1-exp(-1./(f(f>1/smootherwidth)-smootherwidth).^2));
    % x(f>1/smootherwidth) = x(f>1/smootherwidth)./(1+(f(f>1/smootherwidth)-smootherwidth).^2);
    x(1)=0;
    x(end+1:end+length(x)-1) = fliplr(conj(x(2:end)));
   
    x = ifft(x);

end

