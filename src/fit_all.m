
function fit_all(bidsDir)

derivsDir = fullfile(bidsDir,'derivatives','qMRLab');

if ~exist(derivsDir, 'dir')
    mkdir(derivsDir);
 end

% Load MP2RAGE qmrlab protocol for fitting
mp2rage_prot = load('protocol.mat');
mp2rage_prot = mp2rage_prot.Protocol;

subs = struct();
subs.phantom.sessions = ["neutPRI","neutSKY","natSKY"];
subs.invivo.sessions = ["neutSKY","natSKY"];

% Just for the sake of it.
flds = fields(subs);
for i=1:length(flds)
    cur_sessions = subs.(flds{i}).sessions;
    subs.(flds{i}).ks = [];
    subs.(flds{i}).out = [];
    for j = 1:length(cur_sessions)
        cursub = "sub-" + flds{i};
        curses  = "ses-" + cur_sessions(j);
        common = cursub + "_" + curses + "_inv-";
        indir = fullfile(bidsDir,cursub,curses,'anat');
        outdir = fullfile(derivsDir,cursub,curses,'anat');

        if ~exist(outdir, 'dir')
            mkdir(outdir);
         end

        data = struct();
        % Please don't hate me for this.
        data.INV1mag = double(load_nii(char(fullfile(indir,common + num2str(1) + "_part-mag_MP2RAGE.nii.gz"))).img);
        data.INV1phase = double(load_nii(char(fullfile(indir,common + num2str(1) + "_part-phase_MP2RAGE.nii.gz"))).img);
        data.INV2mag = double(load_nii(char(fullfile(indir,common + num2str(2) + "_part-mag_MP2RAGE.nii.gz"))).img);
        data.INV2phase = double(load_nii(char(fullfile(indir,common + num2str(2) + "_part-phase_MP2RAGE.nii.gz"))).img);

        % Set preset prots 
        Model = mp2rage;
        Model.Prot = mp2rage_prot;

        FitResults = FitData(data,Model,0);

        FitResultsSave_nii(FitResults, char(fullfile(indir,common + num2str(1) + "_part-mag_MP2RAGE.nii.gz")),char(outdir));

        addField = struct();
        addField.EstimationAlgorithm =  'qMRLab MP2RAGE';
        addField.Protocol =  mp2rage_prot;
        addField.BasedOn = [{char(fullfile(indir,common + num2str(1) + "_part-mag_MP2RAGE.nii.gz"))};
                            {char(fullfile(indir,common + num2str(1) + "_part-phase_MP2RAGE.nii.gz"))};
                            {char(fullfile(indir,common + num2str(2) + "_part-mag_MP2RAGE.nii.gz"))};
                            {char(fullfile(indir,common + num2str(2) + "_part-phase_MP2RAGE.nii.gz"))}
                            ];
        provenance = Model.getProvenance('extra',addField);

        save_name = char(fullfile(outdir, "qmrlab_provenance.json"));
        savejson('',provenance,save_name);
    end
end

disp('DONE');
end