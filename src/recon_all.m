
function recon_all(bidsDir)

subs = struct();
%subs.phantom.sessions = ["neutPRI","neutSKY"];
subs.phantom.sessions = ["neutPRI","neutSKY","neut750"];
subs.invivo.sessions = ["neutSKY"];

flds = fields(subs);
for i=1:length(flds)
    cur_sessions = subs.(flds{i}).sessions;
    subs.(flds{i}).ks = [];
    subs.(flds{i}).out = [];
    for j = 1:length(cur_sessions)
        cursub = "sub-" + flds{i};
        curses  = "ses-" + cur_sessions(j);
        if ~isfile(fullfile(bidsDir,cursub,curses,'ks', cursub + '_' + curses + '_MP2RAGE.dat'))
            disp('.dat file not found, assuming .mat');
            subs.(flds{i}).ks = [subs.(flds{i}).ks, fullfile(bidsDir,cursub,curses,'ks', cursub + '_' + curses + '_MP2RAGE.mat')];
        else
            subs.(flds{i}).ks = [subs.(flds{i}).ks, fullfile(bidsDir,cursub,curses,'ks', cursub + '_' + curses + '_MP2RAGE.dat')];
        end
        subs.(flds{i}).out = [subs.(flds{i}).out, fullfile(bidsDir,cursub,curses,'anat')];
    end
end

for i=1:length(flds)
    for j = 1:length(subs.(flds{i}).ks)
        reconmp2rage(char(subs.(flds{i}).ks(j)),char(subs.(flds{i}).out(j)));
    end
end

disp('DONE');
end