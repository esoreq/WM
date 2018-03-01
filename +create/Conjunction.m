clear; close all;clc;
root = '/Users/eyalsoreq/Dropbox/PhD/CURRENT/papers/WM/wm-decoding-project/results/report';
outDir = [root,'/Conjunction'];
mkdir(outDir);
fn = c3nl.select('pth','/Users/eyalsoreq/Dropbox/PhD/CURRENT/papers/WM/wm-decoding-project/results/Data/conjunction','name','*.nii');
[Y,Go] = load.vol(fn);
C.bb=[47,139,221]./256;
C.o1=[250,64,64]./256;
YS1 = Y{1};
YT1 = Y{2};
yg = min(cat(4,YS1,YT1),[],4);
m1 = bwareaopen(YT1>0,50,6);
m2 = bwareaopen(YS1>0,50,6);
ds1 = cat(4,m1&m2,m1&~m2).*cat(4,yg,yg);
plot.layersMultiEffect(ds1,Go{1},min(YT1(:)),1,[C.bb;C.o1],{'Domain \cap Stage','Domain \cap Stage^c'},[outDir filesep 'conj_domain_general_task1.ai'],6,6);

YS2 = Y{4};
YT2 = Y{3};
yg = min(cat(4,YT2,YS2),[],4);
m1 = bwareaopen(YT2>0,50,6);
m2 = bwareaopen(YS2>0,50,6);
ds2 = cat(4,m1&m2,m1&~m2).*cat(4,yg,yg);
plot.layersMultiEffect(ds2,Go{1},min(YT2(:)),1,[C.bb;C.o1],{'Domain \cap Stage','Domain \cap Stage^c'},[outDir filesep 'conj_domain_general_task1.ai'],6,6);


QS = @(a,b) 2*(nnz(a&b))/(nnz(a)+nnz(b));
sprintf('Correspondence between stage conjunction across studies is %f',   QS(YS2>0,YS1>0))
sprintf('Correspondence between domain conjunction across studies is %f',   QS(YT2>0,YT1>0))

