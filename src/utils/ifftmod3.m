function out = ifftmod3(in)
    out = ifftshift3(ifft3(ifftshift3(in)));
end
