% IMPORTANT NOTE
% This script is written to handle cases specific to this 
% particular dataset (ds-pulseqmp2rage).
% Metadata, file formats etc. will have to be  properly handled to 
% generalize this solution to any pulseq mp2rage acquisition.

% num_acs:     Size of the calibration region for ESPRIT
%              Affects T1 sanity. Lower values may revert 
%              the distribution. 
% c            ESPRIT mask size. Smaller c --> larger mask
%              This is the reason why reconstructed images 
%              has a noise distribution within a circular FOV.
% num_svd:     Number of SVD channels for coil compression.
%              Larger SVD will take much longer to process.
% Resources:   https://martinos.org/~berkin/2849.html
%              https://www.dropbox.com/s/d295o03z1n2ul9s/SVD_ESPIRiT_Coil_Combine_Toolbox_v2.zip?e=1&dl=0



function reconmp2rage(ks_path,out_path)

[filepath,name,ext] = fileparts(ks_path);

recon_params = struct();

if ext == ".dat"
% TWIX from siemens
    twix_obj = mapVBVD(ks_path);
    image_obj = twix_obj{end}.image;
    sizeData = image_obj.sqzSize; %kx, nch, ky, slices, partitions
    dimsData = image_obj.sqzDims;
    dimsallData = image_obj.dataDims;
    Npar = twix_obj{end}.image.NPar;
    Ncha = twix_obj{end}.image.NCha;
    Npe = twix_obj{end}.image.NLin;
    Nset = twix_obj{end}.image.NSet;
    Ncol = twix_obj{end}.image.NCol;
    data = twix_obj{end}.image(:,:,:,:,:);
    % Permute such that Nsamp, Npe, Npar, Ncha, Nset
    % (i.e., 128  120  96  16  2)
    data  = permute(data,[1,3,4,2,5]);
    % This option works OK for 16 channels, which is the 
    % case with Siemens (both skyra and prisma).
    recon_params.num_svd = 12;
    recon_params.num_acs = 20;
    recon_params.c = 0.2;
elseif ext == ".mat"
% Parsed GE data saved as mat file 
    disp(["Loading " + name]);
    Npe = 120;
    Nset = 2;
    Npar = 96;
    Ncha = 32;
    data = load(ks_path).ks;

    oversampled = false;
    % Oversample along the readout hence 256
    if oversampled
    % Crop along readout
    % !!!!!!! 4X 
    data = reshape(data,[512,Ncha,Npe,Npar,Nset]);
    data = data(1:4:512,:,:,:,:);
    % !!!!!!! 2X 
    data = reshape(data,[256,Ncha,Npe,Npar,Nset]);
    data = data(1:2:256,:,:,:,:);
    else
    Ncha = 32;
    % data = reshape(data,[256,Ncha,Npe,Npar,Nset]);
    % data = data(1:2:256,:,:,:,:);
    end
    % Permute such that Nsamp, Npe, Npar, Ncha, Nset
    data = permute(data,[1,3,4,2,5]);
    recon_params.num_svd = 16;
    %recon_params.num_svd = 24; % There are 48 channels.
    recon_params.num_acs = 16;
    recon_params.c = 0.4;
end

% Parallelize TI and T2 recon.
parfor ii=1:2
    combine_recon(data,ii,recon_params,sprintf('%s%d',name,ii),out_path,false);
end

save_name = char(fullfile(out_path,splitbefore(name) + "_MP2RAGE.json"));
if ext == ".dat"
    savehdr = twix_obj{end}.hdr.Dicom;
else
    savehdr = struct();
end

savehdr.SequenceName = 'Pulseq-mp2rage';
savehdr.RepetitionTimePreparation = 5;
savehdr.RepetitionTimeExcitation = 0.008;
savehdr.InversionTime1 = 0.8;
savehdr.InversionTime2 = 2.1;
savehdr.FlipAngle1 = 7;
savehdr.FlipAngle2 = 5;
savehdr.FieldStrength = 3;
savehdr.NumberShots = [48,48];
savehdr.ReconParams = recon_params;
savehdr.ReconSoftware = 'BART v0.8';
savejson('',savehdr,save_name);

end

function img_combo = coilcomb(img,recon_params,ti,lbl,use_existing)

    num_svd = recon_params.num_svd;
    num_acs = recon_params.num_acs;
    c = recon_params.c;

    N = size(img);
    N = N(1:3);
    temp = reshape(img, [prod(N), size(img,4)]);
    
    [V,D] = eig(temp'*temp);
    V = flipdim(V,2);
    
    % coil compressed image, where 1st chan is the virtual body coil 
    % to be used as phase reference:
    img_svd = reshape(temp * V(:,1:num_svd), [N, num_svd]);
    
    % Default path
    bart_path = 'bart';

    % Alternative path 
    if system(bart_path)~=0
        disp('Using alternative bart path.');
        bart_path = '/opt/local/bin/bart';
    end

writecfl(sprintf('img_svd_%s',lbl), single(img_svd))

if ~use_existing
    % run ESPIRiT:
    disp('Calculating coil sens...');
    disp('fft');
    system([bart_path, ' fft 7 ', sprintf('img_svd_%s ',lbl), sprintf('kspace_svd_%s',lbl)]);
    disp('ecalib');
    system([bart_path, ' ecalib -r ', num2str(num_acs), ' -c ', num2str(c), sprintf(' kspace_svd_%s ',lbl), sprintf('calib_svd_%s',lbl)]);
    disp('slice');
    system([bart_path, ' slice 4 0 ', sprintf('calib_svd_%s ',lbl), sprintf('sens_svd_%s',lbl)]);
end

disp('reading svd');
sens_svd = single(readcfl(sprintf('sens_svd_%s',lbl)));

disp('combination');
% use ESPIRiT sensitivities for Roemer/SENSE coil combination:
img_combo = sum(img_svd .* conj(sens_svd), 4) ./ (eps + sum(abs(sens_svd).^2, 4));

% HARDCODED 
% Set to false if one coil sens map will be used
% to recon more than 1 image.
override_delete = true;

if (use_existing || override_delete)
% Clean up space after 
system([sprintf('rm sens_svd_%s.hdr',lbl)])
system([sprintf('rm sens_svd_%s.cfl',lbl)])
system([sprintf('rm calib_svd_%s.hdr',lbl)])
system([sprintf('rm calib_svd_%s.cfl',lbl)]);
system([sprintf('rm kspace_svd_%s.hdr',lbl)])
system([sprintf('rm kspace_svd_%s.cfl',lbl)]);
system([sprintf('rm img_svd_%s.hdr',lbl)]);
system([sprintf('rm img_svd_%s.cfl',lbl)]);
end

end


function combine_recon(data,ti,recon_params,fname,filepath,use_existing_calib)

    disp(["Reconstructing for " + fname + " TI " + num2str(ti) ]);

    data = squeeze(data(:,:,:,:,ti));
    im = ifftmod3(data);
    
    img_combo = coilcomb(im,recon_params,ti,fname,use_existing_calib);

    nii = make_nii(double(real(img_combo)), [2,2,2], [0,0,0],64);
    save_name = char(fullfile(filepath,splitbefore(fname) + "_inv-" + num2str(ti) + "_part-real_MP2RAGE.nii.gz"));
    save_nii(nii,save_name);

    nii = make_nii(double(imag(img_combo)), [2,2,2], [0,0,0],64);
    save_name = char(fullfile(filepath,splitbefore(fname) + "_inv-" + num2str(ti) + "_part-imag_MP2RAGE.nii.gz"));
    save_nii(nii,save_name);

    nii = make_nii(double(abs(img_combo)), [2,2,2], [0,0,0],64);
    save_name = char(fullfile(filepath,splitbefore(fname) + "_inv-" + num2str(ti) + "_part-mag_MP2RAGE.nii.gz"));
    save_nii(nii,save_name);

    nii = make_nii(double(angle(img_combo)), [2,2,2], [0,0,0],64);
    save_name = char(fullfile(filepath,splitbefore(fname) + "_inv-" + num2str(ti) + "_part-phase_MP2RAGE.nii.gz"));
    save_nii(nii,save_name);

end

function result = splitbefore(name)

    name = char(name);
    idx = max(strfind(name,'_'))-1;
    result = string(name(1:idx));

end


% CANON

% rawData = rawData(1:2:end,:) + 1i*rawData(2:2:end,:);


% GEPARSER

% blockSize = 96;
% numBlocks = 240;

% blocks = cell(1, numBlocks);

% for i = 1:numBlocks
%     startIndex = (i - 1) * blockSize + 1;
%     endIndex = i * blockSize;
%     blocks{i} = rr(:,:,startIndex:endIndex);
% end

% gedat = zeros(128,64,120,96,2);

% idx = 1;
% for ii=1:2:240
% gedat(:,:,idx,:,1) = blocks{ii};
% idx = idx + 1; 
% end

% idx = 1;
% for ii=2:2:240
% gedat(:,:,idx,:,2) = blocks{ii};
% idx = idx + 1; 
% end