function R = covariance(pf,den)
% autocovariance matrix R
tau = (pf.meas_pt(end)-pf.meas_pt(1))/(length(pf.meas_pt)-1);
tau_m = tau*[0 inf inf inf inf inf
    inf 0 inf inf inf inf
    inf inf 0 inf inf inf
    inf inf inf 0 inf inf
    inf inf inf inf 0 inf
    inf inf inf inf inf 0];

R = pf.meas_stdev^2 .* autocorrelate(tau_m,den);
end

