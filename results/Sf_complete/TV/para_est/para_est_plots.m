y= RMSE_list;
x= lambda_list;
p = plot(1:length(x), y);
set(p, 'marker', 'o', 'Color', 'b');
set(p, 'LineWidth', 1, 'markersize', 4, 'MarkerEdgeColor', 'b');
hold on;
p= plot(length(x)-1, y(end-1));
set(p, 'marker', 'o', 'Color', 'r', 'markersize',4,'MarkerFaceColor', 'r');

xlim([1 length(x)+3])
xticks([1,5 10,15 20,25,30,35,40,45]);
xticklabels({'(0.9^{1}','0.9^{5}', '0.9^{10}','0.9^{15}', ...
    '0.9^{20}', '0.9^{25}','0.9^{30}','0.9^{35}','0.9^{40}','0.9^{45})'});
ylabel('RMSE', 'FontSize', 15);
xlabel('\lambda', 'FontSize', 15);
print -depsc2 rmse_vs_lambda.eps
close;

y = SNR_list;
x = lambda_list;
p= plot(1:length(x), y);
set(p, 'marker', 'o', 'Color', 'b');
set(p, 'LineWidth', 1, 'markersize', 4, 'MarkerEdgeColor', 'b');
hold on;
p= plot(length(x)-1, y(end-1));
set(p, 'marker', 'o', 'Color', 'r', 'markersize',4,'MarkerFaceColor', 'r');

xlim([1 length(x)+3])
xticks([1,5 10,15 20,25,30,35,40,45]);
xticklabels({'(0.9^{1}','0.9^{5}', '0.9^{10}','0.9^{15}', ...
    '0.9^{20}', '0.9^{25}','0.9^{30}','0.9^{35}','0.9^{40}','0.9^{45})'});
ylabel('SNR [dB]', 'FontSize', 15);
xlabel('\lambda', 'FontSize', 15);
print -depsc2 snr_vs_lambda.eps
close;

rmse = zeros(2000, 4);
snr  = zeros(2000,4);

for i=1:4
    fname= sprintf(['santaFe_inputs_n_results_rho_1e',num2str(i+1),'.mat']);
    load(fname);
    Nite=length(result.all_rmse);
    rmse(1:Nite,i)=result.all_rmse;
    snr(1:Nite,i)=result.all_snr;
    ind = find(rmse(:,i)==0); rmse(ind,i)=NaN; snr(ind, i)=NaN;
end

max_itr = 2000;
rmse(400:2000,4)=NaN;
figure(1); clf;
h=plot(1:max_itr, rmse(:, 2), 'k');
set(h,'LineWidth',1);
hold on; 
h=plot(1:max_itr, rmse(:, 3),':k');
set(h,'LineWidth', 1);
h=plot(1:max_itr, rmse(:, 4),'--k');
set(h,'LineWidth',1);
max_x = 400;
max_y = max(max(rmse))+.005;
min_y = min(min(rmse))-.005;
axis([ 0 max_x min_y max_y]); 
xlabel('Iterations','Fontsize',15);
ylabel('RMSE','Fontsize',15);
lh=legend({'\rho=1e3', '\rho=1e4', '\rho=1e5'},'Location', 'southeast', ...
    'Fontsize', 12);
lh.Position(2) = 0.5 - lh.Position(4)/2;
print -depsc2 rho_rmse.eps

snr(400:2000,4)=NaN;
max_itr = 2000;
figure(1); clf; 
h=plot(1:max_itr, snr(:, 2),'k');
set(h,'LineWidth',1);
hold on;
h=plot(1:max_itr, snr(:, 3),':k');
set(h,'LineWidth', 1);
hold on
h=plot(1:max_itr, snr(:, 4),'--k');
set(h,'LineWidth',1);
max_x = 400;
max_y = max(max(snr))+2;
min_y = min(min(snr))-2;
xlabel('Iterations','Fontsize',15);
ylabel('SNR [dB]','Fontsize',15);
lh=legend({'\rho=1e3', '\rho=1e4', '\rho=1e5'},'Location', 'southeast', ...
    'Fontsize', 12);
lh.Position(2) = 0.5 - lh.Position(4)/2;
print -depsc2 rho_snr.eps

subplot(2,2,1)
h=plot(1:17, snr(1:17, 1),'k');
set(h,'LineWidth',1);
max_y = max(max(snr))+2;
min_y = min(min(snr))-2;
xlabel('Iterations','Fontsize',12);
ylabel('SNR [dB]','Fontsize',12);
title(rho_string{i}, 'Fontsize', 12, 'FontWeight', 'normal');

subplot(2,2,2)
h=plot(1:50, snr(1:50, 2),'k');
set(h,'LineWidth',1);
max_y = max(max(snr))+2;
min_y = min(min(snr))-2;
xlabel('Iterations','Fontsize',12);
ylabel('SNR [dB]','Fontsize',12);
title(rho_string{2}, 'Fontsize', 12, 'FontWeight', 'normal');

subplot(2,2,3)
h=plot(1:196, snr(1:196, 2),'k');
set(h,'LineWidth',1);
max_y = max(max(snr))+2;
min_y = min(min(snr))-2;
xlabel('Iterations','Fontsize',12);
ylabel('SNR [dB]','Fontsize',12);
title(rho_string{3}, 'Fontsize', 12, 'FontWeight', 'normal');

rho_string={'\rho=1e2','\rho=1e3', '\rho=1e4', '\rho=1e5'};
txt_xg= [10, 12, 40, 1700];
for i=1:4
    fname= sprintf(['santaFe_inputs_n_results_rho_1e',num2str(i+1),'.mat']);
    load(fname);
    Nite=length(result.all_rmse);
    snr=result.all_snr;
    
    subplot(2,2,i)
    h=plot(1:Nite, snr,'k');
    set(h,'LineWidth',1);
    max_y = max(max(snr))+0.5;
    min_y = min(min(snr))-1;
    axis([ 0 Nite min_y max_y]); 
    xlabel('Iterations','Fontsize',10);
    ylabel('SNR [dB]','Fontsize',10);
    title(rho_string{i}, 'Fontsize', 10, 'FontWeight', 'normal');
    %txt=sprintf(['max SNR=', num2str(max(snr))]);
    %text(Nite, max(snr)-1, txt);  
end
print -depsc2 rho_vs_snr_sep.eps

for i=1:4
    fname= sprintf(['santaFe_inputs_n_results_rho_1e',num2str(i+1),'.mat']);
    load(fname);
    Nite=length(result.all_rmse);
    rmse=result.all_rmse;
    
    subplot(2,2,i)
    h=plot(1:Nite, rmse,'k');
    set(h,'LineWidth',1);
    max_y = max(max(rmse))+0.05;
    min_y = min(min(rmse))-0.05;
    axis([ 0 Nite min_y max_y]); 
    xlabel('Iterations','Fontsize',10);
    ylabel('RMSE','Fontsize',10);
    title(rho_string{i}, 'Fontsize', 10, 'FontWeight', 'normal');
    %txt=sprintf(['max RMSE=', num2str(max(rmse))]);
    %t=text(Nite, max(rmse)-0.09, txt);  
    %t.Position(2) = 0.5 - t.Position(4)/2;

end
print -depsc2 rho_vs_rmse_sep.eps
%lambda_min=0.000042
%admm_opts.lambda        = 0.9^50*0.0000001*norm(dataset.A'*dataset.b, inf); %10.747299; %45.02 %32.940500s %8.078978; %1e-8*norm(dataset.A'*dataset.b, inf); 8.078978; %79.872s
%lambda_max=a'b
