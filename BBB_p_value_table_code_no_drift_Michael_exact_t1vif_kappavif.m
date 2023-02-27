

% code for creating p_table with no drift data 
% before running this code, run "BBB_no_drift_MATLAB_presentation_Michael_exact_t1vif_kappavif"

p_table_no_drift=zeros(18,3);

p_table_row_header= [{'SMI-con vs SMI-act'}, {'SMI-con vs MCI-con'}, {'SMI-act vs MCI-con'}]
p_table_column_header= [{'WM_vp'},{'WM_ps'},{'GM_vp'},{'GM_ps'},{'CSF_vp'},{'CSF_ps'},{'Hippo_vp'},{'Hippo_ps'}, {'Thal_vp'},{'Thal_ps'}, {'Amyg_vp'}, {'Amyg_ps'}]; 

[p,h,stats] = ranksum(WM_vp_SMI_con,WM_vp_SMI_act,'alpha',0.01)
p_table_no_drift(1,1)=p;
[p,h,stats] = ranksum(WM_ps_SMI_con,WM_ps_SMI_act,'alpha',0.01)
p_table_no_drift(2,1)=p;
[p,h,stats] = ranksum(WM_T1_SMI_con,WM_T1_SMI_act,'alpha',0.01)
p_table_no_drift(3,1)=p;

[p,h,stats] = ranksum(GM_vp_SMI_con,GM_vp_SMI_act,'alpha',0.01)
p_table_no_drift(4,1)=p;
[p,h,stats] = ranksum(GM_ps_SMI_con,GM_ps_SMI_act,'alpha',0.01)
p_table_no_drift(5,1)=p;
[p,h,stats] = ranksum(GM_T1_SMI_con,GM_T1_SMI_act,'alpha',0.01)
p_table_no_drift(6,1)=p;

[p,h,stats] = ranksum(CSF_vp_SMI_con,CSF_vp_SMI_act,'alpha',0.01)
p_table_no_drift(7,1)=p;
[p,h,stats] = ranksum(CSF_ps_SMI_con,CSF_ps_SMI_act,'alpha',0.01)
p_table_no_drift(8,1)=p;
[p,h,stats] = ranksum(CSF_T1_SMI_con,CSF_T1_SMI_act,'alpha',0.01)
p_table_no_drift(9,1)=p;

[p,h,stats] = ranksum(Hippo_vp_SMI_con,Hippo_vp_SMI_act,'alpha',0.01)
p_table_no_drift(10,1)=p;
[p,h,stats] = ranksum(Hippo_ps_SMI_con,Hippo_ps_SMI_act,'alpha',0.01)
p_table_no_drift(11,1)=p;
[p,h,stats] = ranksum(Hippo_T1_SMI_con,Hippo_T1_SMI_act,'alpha',0.01)
p_table_no_drift(12,1)=p;

[p,h,stats] = ranksum(Thal_vp_SMI_con,Thal_vp_SMI_act,'alpha',0.01)
p_table_no_drift(13,1)=p;
[p,h,stats] = ranksum(Thal_ps_SMI_con,Thal_ps_SMI_act,'alpha',0.01)
p_table_no_drift(14,1)=p;
[p,h,stats] = ranksum(Thal_T1_SMI_con,Thal_T1_SMI_act,'alpha',0.01)
p_table_no_drift(15,1)=p;

[p,h,stats] = ranksum(Amyg_vp_SMI_con,Amyg_vp_SMI_act,'alpha',0.01)
p_table_no_drift(16,1)=p;
[p,h,stats] = ranksum(Amyg_ps_SMI_con,Amyg_ps_SMI_act,'alpha',0.01)
p_table_no_drift(17,1)=p;
[p,h,stats] = ranksum(Amyg_T1_SMI_con,Amyg_T1_SMI_act,'alpha',0.01)
p_table_no_drift(18,1)=p;



[p,h,stats] = ranksum(WM_vp_SMI_con,WM_vp_MCI_con,'alpha',0.01)
p_table_no_drift(1,2)=p;
[p,h,stats] = ranksum(WM_ps_SMI_con,WM_ps_MCI_con,'alpha',0.01)
p_table_no_drift(2,2)=p;
[p,h,stats] = ranksum(WM_T1_SMI_con,WM_T1_MCI_con,'alpha',0.01)
p_table_no_drift(3,2)=p;

[p,h,stats] = ranksum(GM_vp_SMI_con,GM_vp_MCI_con,'alpha',0.01)
p_table_no_drift(4,2)=p;
[p,h,stats] = ranksum(GM_ps_SMI_con,GM_ps_MCI_con,'alpha',0.01)
p_table_no_drift(5,2)=p;
[p,h,stats] = ranksum(GM_T1_SMI_con,GM_T1_MCI_con,'alpha',0.01)
p_table_no_drift(6,2)=p;

[p,h,stats] = ranksum(CSF_vp_SMI_con,CSF_vp_MCI_con,'alpha',0.01)
p_table_no_drift(7,2)=p;
[p,h,stats] = ranksum(CSF_ps_SMI_con,CSF_ps_MCI_con,'alpha',0.01)
p_table_no_drift(8,2)=p;
[p,h,stats] = ranksum(CSF_T1_SMI_con,CSF_T1_MCI_con,'alpha',0.01)
p_table_no_drift(9,2)=p;

[p,h,stats] = ranksum(Hippo_vp_SMI_con,Hippo_vp_MCI_con,'alpha',0.01)
p_table_no_drift(10,2)=p;
[p,h,stats] = ranksum(Hippo_ps_SMI_con,Hippo_ps_MCI_con,'alpha',0.01)
p_table_no_drift(11,2)=p;
[p,h,stats] = ranksum(Hippo_T1_SMI_con,Hippo_T1_MCI_con,'alpha',0.01)
p_table_no_drift(12,2)=p;


[p,h,stats] = ranksum(Thal_vp_SMI_con,Thal_vp_MCI_con,'alpha',0.01)
p_table_no_drift(13,2)=p;
[p,h,stats] = ranksum(Thal_ps_SMI_con,Thal_ps_MCI_con,'alpha',0.01)
p_table_no_drift(14,2)=p;
[p,h,stats] = ranksum(Thal_T1_SMI_con,Thal_T1_MCI_con,'alpha',0.01)
p_table_no_drift(15,2)=p;



[p,h,stats] = ranksum(Amyg_vp_SMI_con,Amyg_vp_MCI_con,'alpha',0.01)
p_table_no_drift(16,2)=p;
[p,h,stats] = ranksum(Amyg_ps_SMI_con,Amyg_ps_MCI_con,'alpha',0.01)
p_table_no_drift(17,2)=p;
[p,h,stats] = ranksum(Amyg_T1_SMI_con,Amyg_T1_MCI_con,'alpha',0.01)
p_table_no_drift(18,2)=p;


[p,h,stats] = ranksum(WM_vp_SMI_act,WM_vp_MCI_con,'alpha',0.01)
p_table_no_drift(1,3)=p;
[p,h,stats] = ranksum(WM_ps_SMI_act,WM_ps_MCI_con,'alpha',0.01)
p_table_no_drift(2,3)=p;
[p,h,stats] = ranksum(WM_T1_SMI_act,WM_T1_MCI_con,'alpha',0.01)
p_table_no_drift(3,3)=p;


[p,h,stats] = ranksum(GM_vp_SMI_act,GM_vp_MCI_con,'alpha',0.01)
p_table_no_drift(4,3)=p;
[p,h,stats] = ranksum(GM_ps_SMI_act,GM_ps_MCI_con,'alpha',0.01)
p_table_no_drift(5,3)=p;
[p,h,stats] = ranksum(GM_T1_SMI_act,GM_T1_MCI_con,'alpha',0.01)
p_table_no_drift(6,3)=p;

[p,h,stats] = ranksum(CSF_vp_SMI_act,CSF_vp_MCI_con,'alpha',0.01)
p_table_no_drift(7,3)=p;
[p,h,stats] = ranksum(CSF_ps_SMI_act,CSF_ps_MCI_con,'alpha',0.01)
p_table_no_drift(8,3)=p;
[p,h,stats] = ranksum(CSF_T1_SMI_act,CSF_T1_MCI_con,'alpha',0.01)
p_table_no_drift(9,3)=p;

[p,h,stats] = ranksum(Hippo_vp_SMI_act,Hippo_vp_MCI_con,'alpha',0.01)
p_table_no_drift(10,3)=p;
[p,h,stats] = ranksum(Hippo_ps_SMI_act,Hippo_ps_MCI_con,'alpha',0.01)
p_table_no_drift(11,3)=p;
[p,h,stats] = ranksum(Hippo_T1_SMI_act,Hippo_T1_MCI_con,'alpha',0.01)
p_table_no_drift(12,3)=p;

[p,h,stats] = ranksum(Thal_vp_SMI_act,Thal_vp_MCI_con,'alpha',0.01)
p_table_no_drift(13,3)=p;
[p,h,stats] = ranksum(Thal_ps_SMI_act,Thal_ps_MCI_con,'alpha',0.01)
p_table_no_drift(14,3)=p;
[p,h,stats] = ranksum(Thal_T1_SMI_act,Thal_T1_MCI_con,'alpha',0.01)
p_table_no_drift(15,3)=p;

[p,h,stats] = ranksum(Amyg_vp_SMI_act,Amyg_vp_MCI_con,'alpha',0.01)
p_table_no_drift(16,3)=p;
[p,h,stats] = ranksum(Amyg_ps_SMI_act,Amyg_ps_MCI_con,'alpha',0.01)
p_table_no_drift(17,3)=p;
[p,h,stats] = ranksum(Amyg_T1_SMI_act,Amyg_T1_MCI_con,'alpha',0.01)
p_table_no_drift(18,3)=p;

row_header={'WM_vp','WM_ps','WM_T1','GM_vp','GM_ps','GM_T1','CSF_vp','CSF_ps','CSF_T1','Hippo_vp','Hippo_ps','Hippo_T1', 'Thal_vp','Thal_ps', 'Thal_T1','Amyg_vp', 'Amyg_ps', 'Amyg_T1'}
col_header= {'SMI_con_vs_SMI_act', 'SMI_con_vs_MCI_con', 'SMI_act_vs_MCI_con'}

xlswrite('BBB_MATLAB_result/p_table_no_drift_Michael_exact_t1vif_kappa_vif.xlsx',p_table_no_drift,'Sheet1','B2');     %Write data
xlswrite('BBB_MATLAB_result/p_table_no_drift_Michael_exact_t1vif_kappa_vif.xlsx',col_header,'Sheet1','B1:D1');     %Write column header
xlswrite('BBB_MATLAB_result/p_table_no_drift_Michael_exact_t1vif_kappa_vif.xlsx',row_header,'Sheet1','A2:A19');      %Write row header



close all