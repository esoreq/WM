clear; close all;clc;
root = '/Users/eyalsoreq/Dropbox/PhD/CURRENT/papers/WM/wm-decoding-project/results/report';
outDir = [root,'/Parcelation'];
mkdir(outDir);
fn = c3nl.select('pth','/Users/eyalsoreq/Dropbox/PhD/CURRENT/papers/WM/wm-decoding-project/results/Data/conjunction','name','*.nii');
[Y,Go] = load.vol(fn);
YS1 = Y{1};
YT1 = Y{2};
LS = atlas.watershed(YS1,1,50);
LT = atlas.watershed(YT1,1,50);
TT = atlas.toTable(LT,Go{1});
TS = atlas.toTable(LS,Go{1});
TT = TT(:,[1:7,13,16]);
TS = TS(:,[1:7,13,16]);

caption = ['Data-drive parcelation set based on Domain general conjunction.'];
plot.mdl2Table(TT,[outDir,'/DGset.tex'],'latex',caption,'DGset');
caption = ['Data-drive parcelation set based on Stage general conjunction.'];
plot.mdl2Table(TS,[outDir,'/SGset.tex'],'latex',caption,'SGset');
save.vol(LT,Go{1},[outDir,'/DGset.nii'],'nii','uint16');
save.vol(LS,Go{1},[outDir,'/SGset.nii'],'nii','uint16');

save([outDir,'/DGset.mat'],'TT')
save([outDir,'/SGset.mat'],'TS')
