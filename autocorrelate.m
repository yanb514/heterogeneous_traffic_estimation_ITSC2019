function R = autocorrelate(tau,den)
R = exp(-abs(tau)/den);
end

