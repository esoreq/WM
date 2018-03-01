%% load beahvioural data, perform repeated measures anova on all tests 
clear; close all;clc;
load('/Users/eyalsoreq/Dropbox/PhD/CURRENT/papers/WM/wm-decoding-project/results/Data/behaviour/results.mat');
root = '/Users/eyalsoreq/Dropbox/PhD/CURRENT/papers/WM/wm-decoding-project/results/report';
outDir = [root,'/behavioural'];
mkdir(outDir);


%% Accuracies repaeted measures anova with post hoc on load for study 1

wideT = Study1.Full_acc;
factors = table(Study1.acc.Load,Study1.acc.Stimulus);
factors.Properties.VariableNames = {'Load','Domain'};
factors.Load = categorical(factors.Load,unique(factors.Load),{'Low','Med','High'});
mes = wideT.Properties.VariableNames(2:end);
rm = fitrm(wideT,sprintf('%s-%s~1',mes{1},mes{end}),'WithinDesign',factors);
mauchly(rm)
epsilon(rm)
% using Greenhouse Geisser correction 
tbl_acc = ranova(rm,'WithinModel',strjoin({'Load','Domain'},'*'));

tbl_acc_posthoc = c3nl.addt(multcompare(rm,'Load'));
mean(Study1.acc.mean_Response)
std(Study1.acc.mean_Response)

caption = ['Accuracy was examined using a 3 x 3 (domain x load) repeated-measures ANOVA',...
    'There was a significant main effect of load ($F_{2,36}=11.85, p_{GG}<0.0001$)'...
    'There was no significant main effect of domain and no significant interaction between load and domain',...
    '($F_{2,36}=2.667, p=0.083;F_{4,72}=1.7515, p=0.106$)'];

plot.mdl2Table(c3nl.roundTable(tbl_acc(:,[1:4,6]),4),[outDir,'/rmANOVA_acc_study1.tex'],'latex',caption,'behStudy1Acc')

caption = ['Post-hoc multicomparison tests were perfomred to assess the direction of the load effects',...
           'significantly lower accuracy during high compared to both low ($p<0.001$) and medium load ($p<0.018$) was found.'];

plot.mdl2Table(c3nl.roundTable(tbl_acc_posthoc,4),[outDir,'/rmANOVA_acc_ph_study1.tex'],'latex',caption,'behPHStudy1Acc')




wideT = Study2.Full_acc;
factors = Study2.acc(:,1:3);
factors.Load = categorical(factors.Load,unique(factors.Load),{'Low','High'});
factors.Stimulus = categorical(factors.Stimulus,{'Pos','Num'},{'Pos','Num'});
factors.Manipulation = categorical(factors.Manipulation,{'Off','On'},{'Off','On'});
mes = wideT.Properties.VariableNames(2:end);
rm = fitrm(wideT,sprintf('%s-%s~1',mes{1},mes{end}),'WithinDesign',factors);
mauchly(rm)
epsilon(rm)
% using Greenhouse Geisser correction 
tbl_acc = ranova(rm,'WithinModel',strjoin({'Load','Stimulus','Manipulation'},'*'));
tbl_acc_posthoc = c3nl.addt(multcompare(rm,'Load'));
mean(Study2.acc.mean_Response)
std(Study2.acc.mean_Response)

caption = ['Accuracy was examined using a 2 x 2 X 2 (domain x load x manipulation) repeated-measures ANOVA',...
    'There was a significant main effect of load ($F_{1,15}=31.532, p<0.00001$)'...
    'There was no significant main effect of domain or manipulation and no significant interactions'];

plot.mdl2Table(c3nl.roundTable(tbl_acc(:,[1:4,6]),4),[outDir,'/rmANOVA_acc_study2.tex'],'latex',caption,'behStudy2Acc')

caption = ['Post-hoc multicomparison tests were perfomred to assess the direction of the load effects',...
           'significantly lower accuracy during high compared was found ($p<0.001$).'];

plot.mdl2Table(c3nl.roundTable(tbl_acc_posthoc,4),[outDir,'/rmANOVA_acc_ph_study2.tex'],'latex',caption,'behPHStudy2Acc')

tbl_acc_posthoc_inter = multcompare(rm,'Load','By','Stimulus');

caption = ['Post-hoc multicomparison tests were perfomred to assess the direction of the interaction between load and domain effects',...
           'significantly lower accuracy for number condition during high load compared was found ($p<0.0001$).'];

plot.mdl2Table(c3nl.roundTable(tbl_acc_posthoc_inter,4),[outDir,'/rmANOVA_acc_ph_inter_study2.tex'],'latex',caption,'rmANOVA_acc_ph_inter_study2')


%% RT with repaeted measures anova for study 1
factors = Study1.acc(:,1:2);
Study1.RawRT.logRT= log(Study1.RawRT.RT);
Study2.RawRT.logRT= log(Study2.RawRT.RT);
xTm1 = grpstats(Study1.RawRT(Study1.RawRT.fail~='fail',:),{'subj','Load','Type'},{'mean'},'DataVars',{'logRT'});% aggregating over runs
%xTm1.mean_logRT = xTm1.mean_RT/1000;
wideT1 = unstack(table(xTm1.subj,xTm1.Load.*xTm1.Type,xTm1.mean_logRT,'Variablenames',{'subj','group','MeanRT'}),'MeanRT','group');
mes = wideT1.Properties.VariableNames(2:end);
rm = fitrm(wideT1,sprintf('%s-%s~1',mes{1},mes{end}),'WithinDesign',factors);
mauchly(rm)
epsilon(rm)
tbl_rt_study1 = ranova(rm,'WithinModel',strjoin({'Load','Stimulus'},'*'));
tbl_rt_study1_load_posthoc = c3nl.addt(multcompare(rm,'Load'));
tbl_rt_study1_domain_posthoc = c3nl.addt(multcompare(rm,'Stimulus'));
tbl_rt_study1_loadbydomain_posthoc = c3nl.addt(multcompare(rm,'Stimulus','By','Load'));

caption = ['Reaction time was examined using a 3 x 3 (domain x load) repeated-measures ANOVA',...
    'There was a significant main effect of load and domain ($F_{2,36}=168.59, p<0.00001;F_{2,36}=24.507, p<0.00001$)'...
    'There was a significant interaction effect between load and domain',...
    '($F_{4,72}=4.508, p<0.005$)'];

plot.mdl2Table(c3nl.roundTable(tbl_rt_study1(:,1:5),4),[outDir,'/rmANOVA_RT_study1.tex'],'latex',caption,'behStudy1RT')

caption = ['Post-hoc multicomparison tests were perfomred to assess the direction of the load effects',...
           'significantly slower RTs during high compared to both low ($\Delta=1.854s,p<0.0001$) and medium load ($\Delta=0.9s,p<0.0001$) were found.'];

plot.mdl2Table(c3nl.roundTable(tbl_rt_study1_load_posthoc,4),[outDir,'/rmANOVA_rt_ph_load_study1.tex'],'latex',caption,'behStudy1RTphLoad')


caption = ['Post-hoc multicomparison tests were perfomred to assess the direction of the domain effects',...
           'significantly faster RTs during Position compared to both Number ($\Delta=0.949s,p<0.0001$) and Object ($\Delta=0.99s,p<0.0001$) were found.'];

plot.mdl2Table(c3nl.roundTable(tbl_rt_study1_domain_posthoc,4),[outDir,'/rmANOVA_rt_ph_domain_study1.tex'],'latex',caption,'behStudy1RTphdomain')


caption = ['Post-hoc multicomparison tests were perfomred to assess the direction of the domain by load effects',...
           'Multiple significant effects were found, it is evident that as load increases so does the gap between position and both object and number.'];

plot.mdl2Table(c3nl.roundTable(tbl_rt_study1_loadbydomain_posthoc,4),[outDir,'/rmANOVA_rt_ph_loadbydomain_study1.tex'],'latex',caption,'behStudy1RTphloadbydomain')

%% RT with repaeted measures anova for study 2
Study2.RawRT.Manipulation = categorical(c3nl.strDetect(cellstr(Study2.RawRT.Type),'Add|Rot'),[0,1],{'Off','On'});
Study2.RawRT.Type = categorical(c3nl.strDetect(cellstr(Study2.RawRT.Type),'Num|Add'),[0,1],{'Pos','Num'});
xTm2 = grpstats(Study2.RawRT(Study2.RawRT.fail~='fail',:),{'subj','Load','Type','Manipulation'},{'mean'},'DataVars',{'logRT'});% aggregating over runs
%xTm2.mean_RT = xTm2.mean_RT/1000;
wideT2 = unstack(table(xTm2.subj,xTm2.Load.*xTm2.Type.*xTm2.Manipulation,xTm2.mean_logRT,'Variablenames',{'subj','group','mean_logRT'}),'mean_logRT','group');
factors = Study2.acc(:,1:3);
mes = wideT2.Properties.VariableNames(2:end);
rm = fitrm(wideT2,sprintf('%s-%s~1',mes{1},mes{end}),'WithinDesign',factors);
mauchly(rm)
epsilon(rm)
tbl_rt_study2 = ranova(rm,'WithinModel',strjoin({'Load','Stimulus','Manipulation'},'*'));
tbl_rt_study2_load_posthoc = c3nl.addt(multcompare(rm,'Load'));
tbl_rt_study2_loadbydomain_posthoc = c3nl.addt(multcompare(rm,'Load','By','Stimulus'));

%tbl,fn,output,caption,label
caption = ['Reaction time was examined using a 2 x 2 x 2 ( domain x load x manipulation) repeated-measures ANOVA',...
    'There were significant main effect of load and interaction between load and domain ($F_{1,15}=74.262, p<0.00001;F_{1,15}=17.049, p<0.001$)'...
    'There were no other significant effects'];

plot.mdl2Table(c3nl.roundTable(tbl_rt_study2(:,[1:4,6]),4),[outDir,'/rmANOVA_RT_study2.tex'],'latex',caption,'behStudy2RT')


caption = ['Post-hoc multicomparison tests were perfomred to assess the direction of the load effects.',...
           ' Significantly slower RTs during high compared to low ($\Delta=1.375s,p<0.0001$) were found.'];

plot.mdl2Table(c3nl.roundTable(tbl_rt_study2_load_posthoc,4),[outDir,'/rmANOVA_rt_ph_load_study2.tex'],'latex',caption,'behStudy2RTphLoad')

caption = ['Post-hoc multicomparison tests were perfomred to assess the direction of the domain by load effects.',...
           ' Significantly slower RTs during high compared to low regardless of domain were found. However, the differences was twice as long for number compared to location.'];

plot.mdl2Table(c3nl.roundTable(tbl_rt_study2_loadbydomain_posthoc,4),[outDir,'/rmANOVA_rt_ph_loadbydomain_study2.tex'],'latex',caption,'behStudy2RTphloadbydomain')

c3nl.copy(c3nl.select('name','*.tex'),'/Users/eyalsoreq/Dropbox/PhD/CURRENT/papers/WM/wm-decoding-project/supplamentryMaterial')
