%plot experiment data and answer returned by SES
type='Ub2b7.6';
pcs_raw = ['data_PCS/' type '.txt'];
t = readtable(pcs_raw, 'HeaderLines', 2);
res_pcs = [t.Var1,t.Var4];
res_pcs(any(isnan(res_pcs),2)|any(res_pcs==0,2),:) = [];
pcs_exp = res_pcs(:,2);
pcs_exp = pcs_exp(1:length(res_pcs)-3);
A = dlmread(['/Users/yuexichen/Desktop/nmr/experiment/SES/data/pcs/Amat/Amat_' type '.txt'],'\t');
columns = [182	149	12];
weight = [0.67153993	0.24303594	0.19373275];
pcs_ses = sum(A(:, columns) .* weight, 2);
scatter(pcs_exp, pcs_ses);
xlabel("PCS experiment")
ylabel("PCS returned by SES")
num_conf = num2str(length(columns));
title(['reconstruct ' type ' PCS data: ' num_conf ' conformation'])