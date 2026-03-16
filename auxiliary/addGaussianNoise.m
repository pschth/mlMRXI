function [b_noisy, noise] = addGaussianNoise(b, SNR, noise)
% adds gaussian white noise to the signal b with amplitude according to the
% SNR in dB.
% 
% INPUT
% b - vector; signal to be corrupted with noise
% SNR - scalar; the demanded signal-to-noise ratio in dB
% noise - vector (optional); noise signal
% 
% OUTPUT
% b_noisy - vector; signal corrupted with noise
% noise - vector; applied noise signal

if nargin < 3
    noise = randn(size(b));
elseif isempty(noise)
    noise = randn(size(b));
end

noise = noise./norm(noise, 'fro') * norm(b, 'fro') / 10^(SNR/20);
b_noisy = b+noise;
end