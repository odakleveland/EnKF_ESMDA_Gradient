function w = ricker_wavelet_zero_phased(dt,t,centr_frq)

nt = length(t);

w = (1 - 2 .* pi^2 .* centr_frq.^2 .* t.^2) .* exp(-pi^2 .* centr_frq.^2 .* t.^2);

w = w(1:end);

end