function out = ifft3(in)
%FFT3 Summary of this function goes here
%   Detailed explanation goes here
    out = ifft(ifft(ifft(in, [], 1), [], 2), [], 3);
end

