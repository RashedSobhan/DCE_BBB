BBB_no_drift=xlsread('BBB_results_ADA_no_drift_on_T1_std_enh_timed100_Michael_exact_t1vif_kappavif.xlsx')

% plot for WM_vp (top left, 3 by 3) 

subplot(3,3,1)
WM_vp_SMI_con=100.*BBB_no_drift([1:8],7)
WM_vp_SMI_act=100.*BBB_no_drift([9:17],7)
WM_vp_MCI_con=100.*BBB_no_drift([18:22],7)



group = [repmat({'SMI-con'}, 8, 1); repmat({'SMI-act'}, 9, 1); repmat({'MCI-con'}, 5, 1)];
WM_vp_boxplot=boxplot([WM_vp_SMI_con;WM_vp_SMI_act;WM_vp_MCI_con], group)
set(WM_vp_boxplot, 'linewidth' ,2)

x= get(gca, 'XTick')
hold on 

scatter (repmat(x(1),8,1), WM_vp_SMI_con, 'mo','linewidth', 1)
scatter (repmat(x(2),9,1), WM_vp_SMI_act, 'ro','linewidth', 1)
scatter (repmat(x(3),5,1), WM_vp_MCI_con, 'k*','linewidth', 1)

title ('WM v_p')
ylabel('v_p(\times 10^-^2^ )')
ylim([-1 4])

% plot for WM_ps (top middle, 3 by 3) 
subplot(3,3,2)

WM_ps_SMI_con=1000*BBB_no_drift([1:8],8)
WM_ps_SMI_act=1000*BBB_no_drift([9:17],8)
WM_ps_MCI_con=1000*BBB_no_drift([18:22],8)



group = [repmat({'SMI-con'}, 8, 1); repmat({'SMI-act'}, 9, 1); repmat({'MCI-con'}, 5, 1)];
WM_ps_boxplot=boxplot([WM_ps_SMI_con;WM_ps_SMI_act;WM_ps_MCI_con], group)
set(WM_ps_boxplot, 'linewidth' ,2)

x= get(gca, 'XTick')
hold on 

scatter (repmat(x(1),8,1), WM_ps_SMI_con, 'mo','linewidth', 1)
scatter (repmat(x(2),9,1), WM_ps_SMI_act, 'ro','linewidth', 1)
scatter (repmat(x(3),5,1), WM_ps_MCI_con, 'k*','linewidth', 1)

title ('WM PS')
ylabel('PS(\times 10^-^3^ min^-^1^ )')
ylim([-1 2])

% plot for WM_T1 (top right, 3 by 3) 

subplot(3,3,3)
WM_T1_SMI_con=BBB_no_drift([1:8],9)
WM_T1_SMI_act=BBB_no_drift([9:17],9)
WM_T1_MCI_con=BBB_no_drift([18:22],9)



group = [repmat({'SMI-con'}, 8, 1); repmat({'SMI-act'}, 9, 1); repmat({'MCI-con'}, 5, 1)];
WM_T1_boxplot=boxplot([WM_T1_SMI_con;WM_T1_SMI_act;WM_T1_MCI_con], group)
set(WM_T1_boxplot, 'linewidth' ,2)

x= get(gca, 'XTick')
hold on 

scatter (repmat(x(1),8,1), WM_T1_SMI_con, 'mo','linewidth', 1)
scatter (repmat(x(2),9,1), WM_T1_SMI_act, 'ro','linewidth', 1)
scatter (repmat(x(3),5,1), WM_T1_MCI_con, 'k*','linewidth', 1)

title ('WM T1')
ylabel('T1(ms)')

% plot for GM_vp (middle left, 3 by 3) 

subplot(3,3,4)
GM_vp_SMI_con=100*BBB_no_drift([1:8],12)
GM_vp_SMI_act=100*BBB_no_drift([9:17],12)
GM_vp_MCI_con=100*BBB_no_drift([18:22],12)



group = [repmat({'SMI-con'}, 8, 1); repmat({'SMI-act'}, 9, 1); repmat({'MCI-con'}, 5, 1)];
GM_vp_boxplot=boxplot([GM_vp_SMI_con;GM_vp_SMI_act;GM_vp_MCI_con], group)
set(GM_vp_boxplot, 'linewidth' ,2)

x= get(gca, 'XTick')
hold on 

scatter (repmat(x(1),8,1), GM_vp_SMI_con, 'mo','linewidth', 1)
scatter (repmat(x(2),9,1), GM_vp_SMI_act, 'ro','linewidth', 1)
scatter (repmat(x(3),5,1), GM_vp_MCI_con, 'k*','linewidth', 1)

title ('GM v_p')
ylabel('v_p(\times 10^-^2^ )')
ylim([0 7])

% plot for GM_ps (middle of the middle, 3 by 3) 
subplot(3,3,5)

GM_ps_SMI_con=1000*BBB_no_drift([1:8],13)
GM_ps_SMI_act=1000*BBB_no_drift([9:17],13)
GM_ps_MCI_con=1000*BBB_no_drift([18:22],13)



group = [repmat({'SMI-con'}, 8, 1); repmat({'SMI-act'}, 9, 1); repmat({'MCI-con'}, 5, 1)];
GM_ps_boxplot=boxplot([GM_ps_SMI_con;GM_ps_SMI_act;GM_ps_MCI_con], group)
set(GM_ps_boxplot, 'linewidth' ,2)

x= get(gca, 'XTick')
hold on 

scatter (repmat(x(1),8,1), GM_ps_SMI_con, 'mo','linewidth', 1)
scatter (repmat(x(2),9,1), GM_ps_SMI_act, 'ro','linewidth', 1)
scatter (repmat(x(3),5,1), GM_ps_MCI_con, 'k*','linewidth', 1)

title ('GM PS')
ylabel('PS(\times 10^-^3^ min^-^1)')
 ylim([-2 6])

% plot for GM_T1 (middle right, 3 by 3) 

subplot(3,3,6)
GM_T1_SMI_con=BBB_no_drift([1:8],14)
GM_T1_SMI_act=BBB_no_drift([9:17],14)
GM_T1_MCI_con=BBB_no_drift([18:22],14)



group = [repmat({'SMI-con'}, 8, 1); repmat({'SMI-act'}, 9, 1); repmat({'MCI-con'}, 5, 1)];
GM_T1_boxplot=boxplot([GM_T1_SMI_con;GM_T1_SMI_act;GM_T1_MCI_con], group)
set(GM_T1_boxplot, 'linewidth' ,2)

x= get(gca, 'XTick')
hold on 

scatter (repmat(x(1),8,1), GM_T1_SMI_con, 'mo','linewidth', 1)
scatter (repmat(x(2),9,1), GM_T1_SMI_act, 'ro','linewidth', 1)
scatter (repmat(x(3),5,1), GM_T1_MCI_con, 'k*','linewidth', 1)

title ('GM T1')
ylabel('T1(ms)')

% CSF_vp plot (bottom left, 3 by 3)
subplot(3,3,7)
CSF_vp_SMI_con=100*BBB_no_drift([1:8],17)
CSF_vp_SMI_act=100*BBB_no_drift([9:17],17)
CSF_vp_MCI_con=100*BBB_no_drift([18:22],17)



group = [repmat({'SMI-con'}, 8, 1); repmat({'SMI-act'}, 9, 1); repmat({'MCI-con'}, 5, 1)];
CSF_vp_boxplot=boxplot([CSF_vp_SMI_con;CSF_vp_SMI_act;CSF_vp_MCI_con], group)
set(CSF_vp_boxplot, 'linewidth' ,2)

x= get(gca, 'XTick')
hold on 

scatter (repmat(x(1),8,1), CSF_vp_SMI_con, 'mo','linewidth', 1)
scatter (repmat(x(2),9,1), CSF_vp_SMI_act, 'ro','linewidth', 1)
scatter (repmat(x(3),5,1), CSF_vp_MCI_con, 'k*','linewidth', 1)

title ('CSF v_p')
ylabel('v_p(\times 10^-^2^ )')
 ylim([1 20])

% plot for CSF_ps (bottom middle, 3 by 2) 
subplot(3,3,8)

CSF_ps_SMI_con=1000*BBB_no_drift([1:8],18)
CSF_ps_SMI_act=1000*BBB_no_drift([9:17],18)
CSF_ps_MCI_con=1000*BBB_no_drift([18:22],18)



group = [repmat({'SMI-con'}, 8, 1); repmat({'SMI-act'}, 9, 1); repmat({'MCI-con'}, 5, 1)];
CSF_ps_boxplot=boxplot([CSF_ps_SMI_con;CSF_ps_SMI_act;CSF_ps_MCI_con], group)
set(CSF_ps_boxplot, 'linewidth' ,2)

x= get(gca, 'XTick')
hold on 

scatter (repmat(x(1),8,1), CSF_ps_SMI_con, 'mo','linewidth', 1)
scatter (repmat(x(2),9,1), CSF_ps_SMI_act, 'ro','linewidth', 1)
scatter (repmat(x(3),5,1), CSF_ps_MCI_con, 'k*','linewidth', 1)

title ('CSF PS')
ylabel('PS(\times 10^-^3^ min^-^1)')
ylim([-2 4])

% plot for CSF_T1 ( bottom right, 3 by 3) 

subplot(3,3,9)
CSF_T1_SMI_con=BBB_no_drift([1:8],19)
CSF_T1_SMI_act=BBB_no_drift([9:17],19)
CSF_T1_MCI_con=BBB_no_drift([18:22],19)



group = [repmat({'SMI-con'}, 8, 1); repmat({'SMI-act'}, 9, 1); repmat({'MCI-con'}, 5, 1)];
CSF_T1_boxplot=boxplot([CSF_T1_SMI_con;CSF_T1_SMI_act;CSF_T1_MCI_con], group)
set(WM_T1_boxplot, 'linewidth' ,2)

x= get(gca, 'XTick')
hold on 

scatter (repmat(x(1),8,1), CSF_T1_SMI_con, 'mo','linewidth', 1)
scatter (repmat(x(2),9,1), CSF_T1_SMI_act, 'ro','linewidth', 1)
scatter (repmat(x(3),5,1), CSF_T1_MCI_con, 'k*','linewidth', 1)

title ('CSF T1')
ylabel('T1(ms)')

%saveas(gcf,'BBB_WM+GM_withdrift_Michael.png', '-dpng','-r300')
saveas (gcf, 'BBB_MATLAB_result/BBB_WM+GM+CSF_no_drift_Michael_exact_t1vif_kappa_vif.fig')

exportgraphics(gcf,'BBB_MATLAB_result/BBB_WM+GM+CSF_no_drift_Michael_exact_t1vif_kappavif.png','Resolution',300)

figure 

% plot for Hippo_vp (top left, 3 by 3) 

subplot(3,3,1)
Hippo_vp_SMI_con=100*BBB_no_drift([1:8],22)
Hippo_vp_SMI_act=100*BBB_no_drift([9:17],22)
Hippo_vp_MCI_con=100*BBB_no_drift([18:22],22)



group = [repmat({'SMI-con'}, 8, 1); repmat({'SMI-act'}, 9, 1); repmat({'MCI-con'}, 5, 1)];
Hippo_vp_boxplot=boxplot([Hippo_vp_SMI_con;Hippo_vp_SMI_act;Hippo_vp_MCI_con], group)

set(Hippo_vp_boxplot, 'linewidth' ,2)

x= get(gca, 'XTick')
hold on 

scatter (repmat(x(1),8,1), Hippo_vp_SMI_con, 'mo','linewidth', 1)
scatter (repmat(x(2),9,1), Hippo_vp_SMI_act, 'ro','linewidth', 1)
scatter (repmat(x(3),5,1), Hippo_vp_MCI_con, 'k*','linewidth', 1)

title ('Hippocampus v_p')
ylabel('v_p(\times 10^-^2^ )')
ylim([-1 6])

% plot for Hippo_ps (top right, 3 by 2) 

subplot(3,3,2)
Hippo_ps_SMI_con=1000*BBB_no_drift([1:8],23)
Hippo_ps_SMI_act=1000*BBB_no_drift([9:17],23)
Hippo_ps_MCI_con=1000*BBB_no_drift([18:22],23)



group = [repmat({'SMI-con'}, 8, 1); repmat({'SMI-act'}, 9, 1); repmat({'MCI-con'}, 5, 1)];


Hippo_ps_boxplot=boxplot([Hippo_ps_SMI_con;Hippo_ps_SMI_act;Hippo_ps_MCI_con], group)
set(Hippo_ps_boxplot, 'linewidth' ,2)

x= get(gca, 'XTick')
hold on 

scatter (repmat(x(1),8,1), Hippo_ps_SMI_con, 'mo','linewidth', 1)
scatter (repmat(x(2),9,1), Hippo_ps_SMI_act, 'ro','linewidth', 1)
scatter (repmat(x(3),5,1), Hippo_ps_MCI_con, 'k*','linewidth', 1)

title ('Hippocampus PS')
ylabel('PS(\times 10^-^3^ min^-^1)')
ylim([-2 4])


% plot for Hippo_T1 (top right, 3 by 3) 

subplot(3,3,3)
Hippo_T1_SMI_con=BBB_no_drift([1:8],24)
Hippo_T1_SMI_act=BBB_no_drift([9:17],24)
Hippo_T1_MCI_con=BBB_no_drift([18:22],24)



group = [repmat({'SMI-con'}, 8, 1); repmat({'SMI-act'}, 9, 1); repmat({'MCI-con'}, 5, 1)];
Hippo_T1_boxplot=boxplot([Hippo_T1_SMI_con;Hippo_T1_SMI_act;Hippo_T1_MCI_con], group)
set(WM_T1_boxplot, 'linewidth' ,2)

x= get(gca, 'XTick')
hold on 

scatter (repmat(x(1),8,1), Hippo_T1_SMI_con, 'mo','linewidth', 1)
scatter (repmat(x(2),9,1), Hippo_T1_SMI_act, 'ro','linewidth', 1)
scatter (repmat(x(3),5,1), Hippo_T1_MCI_con, 'k*','linewidth', 1)

title ('Hippocampus T1')
ylabel('T1(ms)')

% plot for Thalamus_vp (middle left, 3 by 3) 

subplot(3,3,4)
Thal_vp_SMI_con=100*BBB_no_drift([1:8],27)
Thal_vp_SMI_act=100*BBB_no_drift([9:17],27)
Thal_vp_MCI_con=100*BBB_no_drift([18:22],27)



group = [repmat({'SMI-con'}, 8, 1); repmat({'SMI-act'}, 9, 1); repmat({'MCI-con'}, 5, 1)];
Thal_vp_boxplot=boxplot([Thal_vp_SMI_con;Thal_vp_SMI_act;Thal_vp_MCI_con], group)
set(Thal_vp_boxplot, 'linewidth' ,2)

x= get(gca, 'XTick')
hold on 

scatter (repmat(x(1),8,1), Thal_vp_SMI_con, 'mo','linewidth', 1)
scatter (repmat(x(2),9,1), Thal_vp_SMI_act, 'ro','linewidth', 1)
scatter (repmat(x(3),5,1), Thal_vp_MCI_con, 'k*','linewidth', 1)

title ('Thalamus v_p')
ylabel('v_p(\times 10^-^2^ )')
ylim([-1 6])

% plot for Thalamus_ps (middle middle, 3 by 2) 

subplot(3,3,5)
Thal_ps_SMI_con=1000*BBB_no_drift([1:8],28)
Thal_ps_SMI_act=1000*BBB_no_drift([9:17],28)
Thal_ps_MCI_con=1000*BBB_no_drift([18:22],28)



group = [repmat({'SMI-con'}, 8, 1); repmat({'SMI-act'}, 9, 1); repmat({'MCI-con'}, 5, 1)];
Thal_ps_boxplot=boxplot([Thal_ps_SMI_con;Thal_ps_SMI_act;Thal_ps_MCI_con], group)
set(Thal_ps_boxplot, 'linewidth' ,2)

x= get(gca, 'XTick')
hold on 

scatter (repmat(x(1),8,1), Thal_ps_SMI_con, 'mo','linewidth', 1)
scatter (repmat(x(2),9,1), Thal_ps_SMI_act, 'ro','linewidth', 1)
scatter (repmat(x(3),5,1), Thal_ps_MCI_con, 'k*','linewidth', 1)

title ('Thalamus PS')
ylabel('PS(\times 10^-^3^ min^-^1)')
ylim([-2 4])


% plot for Thalamus_T1 (middle right, 3 by 3) 

subplot(3,3,6)
Thal_T1_SMI_con=BBB_no_drift([1:8],29)
Thal_T1_SMI_act=BBB_no_drift([9:17],29)
Thal_T1_MCI_con=BBB_no_drift([18:22],29)



group = [repmat({'SMI-con'}, 8, 1); repmat({'SMI-act'}, 9, 1); repmat({'MCI-con'}, 5, 1)];
Thal_T1_boxplot=boxplot([Thal_T1_SMI_con;Thal_T1_SMI_act;Thal_T1_MCI_con], group)
set(WM_T1_boxplot, 'linewidth' ,2)

x= get(gca, 'XTick')
hold on 

scatter (repmat(x(1),8,1), Thal_T1_SMI_con, 'mo','linewidth', 1)
scatter (repmat(x(2),9,1), Thal_T1_SMI_act, 'ro','linewidth', 1)
scatter (repmat(x(3),5,1), Thal_T1_MCI_con, 'k*','linewidth', 1)

title ('Thalamus T1')
ylabel('T1(ms)')


% Amygdala_vp plot (bottom left, 3 by 2)
subplot(3,3,7)
Amyg_vp_SMI_con=100*BBB_no_drift([1:8],32)
Amyg_vp_SMI_act=100*BBB_no_drift([9:17],32)
Amyg_vp_MCI_con=100*BBB_no_drift([18:22],32)



group = [repmat({'SMI-con'}, 8, 1); repmat({'SMI-act'}, 9, 1); repmat({'MCI-con'}, 5, 1)];
Amyg_vp_boxplot=boxplot([Amyg_vp_SMI_con;Amyg_vp_SMI_act;Amyg_vp_MCI_con], group)
set(Amyg_vp_boxplot, 'linewidth' ,2)

x= get(gca, 'XTick')
hold on 

scatter (repmat(x(1),8,1), Amyg_vp_SMI_con, 'mo','linewidth', 1)
scatter (repmat(x(2),9,1), Amyg_vp_SMI_act, 'ro','linewidth', 1)
scatter (repmat(x(3),5,1), Amyg_vp_MCI_con, 'k*','linewidth', 1)

title ('Amygdala v_p')
ylabel('v_p(\times 10^-^2^ )')
ylim([-1 6])


% Amygdala_ps plot (bottom right, 3 by 2)
subplot(3,3,8)
Amyg_ps_SMI_con=1000*BBB_no_drift([1:8],33)
Amyg_ps_SMI_act=1000*BBB_no_drift([9:17],33)
Amyg_ps_MCI_con=1000*BBB_no_drift([18:22],33)



group = [repmat({'SMI-con'}, 8, 1); repmat({'SMI-act'}, 9, 1); repmat({'MCI-con'}, 5, 1)];
Amyg_ps_boxplot=boxplot([Amyg_ps_SMI_con;Amyg_ps_SMI_act;Amyg_ps_MCI_con], group)
set(Amyg_ps_boxplot, 'linewidth' ,2)

x= get(gca, 'XTick')
hold on 

scatter (repmat(x(1),8,1), Amyg_ps_SMI_con, 'mo','linewidth', 1)
scatter (repmat(x(2),9,1), Amyg_ps_SMI_act, 'ro','linewidth', 1)
scatter (repmat(x(3),5,1), Amyg_ps_MCI_con, 'k*','linewidth', 1)

title ('Amygdala PS')
ylabel('PS(\times 10^-^3^ min^-^1)')
ylim([-2 4])

% plot for Amygdala_T1 (middle right, 3 by 3) 

subplot(3,3,9)
Amyg_T1_SMI_con=BBB_no_drift([1:8],34)
Amyg_T1_SMI_act=BBB_no_drift([9:17],34)
Amyg_T1_MCI_con=BBB_no_drift([18:22],34)



group = [repmat({'SMI-con'}, 8, 1); repmat({'SMI-act'}, 9, 1); repmat({'MCI-con'}, 5, 1)];
Amyg_T1_boxplot=boxplot([Amyg_T1_SMI_con;Amyg_T1_SMI_act;Amyg_T1_MCI_con], group)
set(WM_T1_boxplot, 'linewidth' ,2)

x= get(gca, 'XTick')
hold on 

scatter (repmat(x(1),8,1), Amyg_T1_SMI_con, 'mo','linewidth', 1)
scatter (repmat(x(2),9,1), Amyg_T1_SMI_act, 'ro','linewidth', 1)
scatter (repmat(x(3),5,1), Amyg_T1_MCI_con, 'k*','linewidth', 1)

title ('Amygdala T1')
ylabel('T1(ms)')
saveas (gcf, 'BBB_MATLAB_result/BBB_Hippo+Thal+Amyg_no_drift_Michael_exact_t1vif_kappa_vif.fig')

exportgraphics(gcf,'BBB_MATLAB_result/BBB_Hippo+Thal+Amyg_no_drift_Michael_exact_t1vif_kappa_vif.png','Resolution',300)