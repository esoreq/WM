clear; close all;clc;
root = '/Users/eyalsoreq/Dropbox/PhD/CURRENT/papers/WM/wm-decoding-project/results';
inputDir = [root,'/Data/MinedData'];
outDir = [root,'/report/HCL'];
mkdir(outDir)
load([inputDir '/DG_gpp.mat']);

X1 = grpstats(data.task1.BA,{'subj','domain','state'},{'mean'},'DataVars',data.task1.BA.Properties.VariableNames(8:end));% aggregating over dfact
X2 = grpstats(data.task2.BA,{'subj','domain','state'},{'mean'},'DataVars',data.task2.BA.Properties.VariableNames(9:end));% aggregating over dfact
ix1 = ~c3nl.strDetect(X1.state,'Cue');
ix2 = ~c3nl.strDetect(X2.state,'Cue');
Y2 = categorical([X1.domain(ix1);X2.domain(ix2)]);
Y3 = categorical([X1.state(ix1);X2.state(ix2)]);
load([root '/report/Parcelation/DGset.mat']);

ixh = 1:height(TT);
xx = [X1{ix1,5:end};X2{ix2,5:end}];
cmap = [188,209,244;73,139,234;242,97,43;224,135,184;92,186,71;234,61,106]./255;
clf;[measuresOrder,featuresOrder,cid,measuresClusters,croi]=plot.HCLvert(xx,...
		struct('clim',[-3,3],'Y1',[ones(nnz(ix1),1);zeros(nnz(ix2),1)],...
		'Y2',double(Y2),'Y3',double(Y3),...
		'xlabel',{arrayfun(@(a,b) sprintf('%02d. %s',a,b{1}),ixh', TT.name_AAL2(ixh),'un',0)},...
		'Yc1',struct('type','Study','cat',categorical({'Study. 1','Study. 2'})),...
		'Yc2',struct('type','Domain','cat',unique(Y2)),...
		'Yc3',struct('type','Stage','cat',unique(Y3))),'euclidean','ward',cmap([4,5,6,1,2],:));
    
save.pdf( [outDir filesep 'domainGeneral_hcl.ai'],30,21);close all
[L,G] = load.vol([root,'/report/Parcelation/DGset.nii']);
plot.roi(L,G,TT,cid,cmap([4,5,6,1,2],:),1,9)
save.pdf( [outDir filesep 'domainGeneral_roi_hcl.ai'],34,12);close all
T = table();
T.cid = cid;
T.name = matlab.lang.makeUniqueStrings(TT.name_AAL2);
T.number = TT.ROIid;
plot.labels(T,500,800,2,cmap([4,5,6,1,2],:).^2)
save.pdf( [outDir filesep 'domainGeneral_roi_labels.ai'],6,12);close all
save([outDir filesep 'DGT.mat'],'T');

M = crosstab(grp2idx(Y3),measuresClusters);
ix = c3nl.assignment(M,'cols','max'); % match label to maximum count
a = dummyvar(measuresClusters);
a = sum(a(:,ix).*(1:3),2);
yfit = categorical(a,1:3,categories(Y3));
plot.confusiongrad(Y3,yfit,categories(Y3),'',600,600,0.05,1)
save.pdf( [outDir filesep 'domainGeneral_purity.ai'],5,5);close all

