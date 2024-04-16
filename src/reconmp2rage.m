function reconmp2rage(ks_path,out_path)

[filepath,name,ext] = fileparts(ks_path);
twix_obj = mapVBVD(ks_path);

image_obj = twix_obj{end}.image;
sizeData = image_obj.sqzSize; %kx, nch, ky, slices, partitions
dimsData = image_obj.sqzDims;
dimsallData = image_obj.dataDims;

Npar = twix_obj{end}.image.NPar;
Ncha = twix_obj{end}.image.NCha;
Npe = twix_obj{end}.image.NLin;
Nset = twix_obj{end}.image.NSet;
Ncol = twix_obj{end}.image.NCol; %sample 

data = twix_obj{end}.image(:,:,:,:,:);

data  = permute(data,[1,3,4,2,5]);

ksp_filter = false;

combine_ks2im(data,1,name,out_path);
combine_ks2im(data,2,name,out_path);

save_name = char(fullfile(out_path,splitbefore(name) + "_MP2RAGE.json"));
savejson('',twix_obj{end}.hdr.Dicom,save_name);

end

function img_combo = coilcomb(img)
    
    num_svd = 8;                   % no of SVD channels for compression (num_svd = 16 works well for 32 chan array)

    N = size(img);
    N = N(1:3);
    temp = reshape(img, [prod(N), size(img,4)]);
    
    [V,D] = eig(temp'*temp);
    V = flipdim(V,2);
    
    % coil compressed image, where 1st chan is the virtual body coil to be used as phase reference:
    img_svd = reshape(temp * V(:,1:num_svd), [N, num_svd]);
    
    bart_path = '/opt/local/bin/bart';

    num_acs = 16;                   % size of calibration region for sensitivity estimation (doesn't change the result too much)
    c = 0.4;                        % mask size for coil sensitivities, smaller "c" provides larger mask        

writecfl('img_svd', single(img_svd))

% run ESPIRiT:
system([bart_path, ' fft 7 ', 'img_svd ', 'kspace_svd']);
system([bart_path,' ecalib -r ', num2str(num_acs), ' -c ', num2str(c), ' kspace_svd ', 'calib_svd']);
system([bart_path,' slice 4 0 ', 'calib_svd ', 'sens_svd']);


system('rm calib_svd.hdr');
system('rm calib_svd.cfl');

system('rm kspace_svd.hdr');
system('rm kspace_svd.cfl');

system('rm img_svd.hdr');
system('rm img_svd.cfl');


sens_svd = single(readcfl('sens_svd'));

% use ESPIRiT sensitivities for Roemer/SENSE coil combination:
img_combo = sum(img_svd .* conj(sens_svd), 4) ./ (eps + sum(abs(sens_svd).^2, 4));

% clean up space:
system('rm sens_svd.hdr');
system('rm sens_svd.cfl');

end


function combine_ks2im(data,ti,fname,filepath)

    data = squeeze(data(:,:,:,:,ti));
    im = ifftmod3(data);
    img_combo = coilcomb(im);

    nii = make_nii(double(real(img_combo).*1e6), [2,2,2], [0,0,0],64);
    save_name = char(fullfile(filepath,splitbefore(fname) + "_inv-" + num2str(ti) + "_part-real_MP2RAGE.nii.gz"));
    save_nii(nii,save_name);

    nii = make_nii(double(imag(img_combo).*1e6), [2,2,2], [0,0,0],64);
    save_name = char(fullfile(filepath,splitbefore(fname) + "_inv-" + num2str(ti) + "_part-imag_MP2RAGE.nii.gz"));
    save_nii(nii,save_name);

    nii = make_nii(double(abs(img_combo).*1e6), [2,2,2], [0,0,0],64);
    save_name = char(fullfile(filepath,splitbefore(fname) + "_inv-" + num2str(ti) + "_part-mag_MP2RAGE.nii.gz"));
    save_nii(nii,save_name);

    nii = make_nii(double(angle(img_combo).*1e6), [2,2,2], [0,0,0],64);
    save_name = char(fullfile(filepath,splitbefore(fname) + "_inv-" + num2str(ti) + "_part-phase_MP2RAGE.nii.gz"));
    save_nii(nii,save_name);

end

function result = splitbefore(name)

name = char(name);
idx = max(strfind(name,'_'))-1;
result = string(name(1:idx));

end