function out = fftshift3(in)
%FFTSHIFT3 Summary of this function goes here
%   Detailed explanation goes here
    out = fftshift(fftshift(fftshift(in,1),2),3);
end

