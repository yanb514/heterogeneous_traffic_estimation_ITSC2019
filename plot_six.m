function plot_six(U_true_c1,U_sim_c1,U_est_c1,U_true_c2,U_sim_c2,U_est_c2,model_true)
res = 350;


figure;
plot_sc(U_true_c1,model_true);
% title('True density','interpreter', 'latex');
set(gca,'fontsize',20)
set(gcf,'paperpositionmode','auto')
set(gcf,'position',[100  100 res*1.2 res])

figure;
plot_sc(U_true_c2,model_true);
% title('True density','interpreter', 'latex');
set(gca,'fontsize',20)
set(gcf,'paperpositionmode','auto')
set(gcf,'position',[100  100 res*1.2 res])


% 
% 
% figure;
% subplot(122);plot_sc(U_sim_c1,model_true);
% title('Small vehicle density (simulated)','interpreter', 'latex');
% set(gca,'fontsize',20)
% subplot(121);plot_sc(U_sim_c2,model_true);
% title('Large vehicle density (simulated)','interpreter', 'latex');
% set(gca,'fontsize',20)
% res = 1400;
% set(gcf,'paperpositionmode','auto')
% set(gcf,'position',[50  50 res res/3])
%

figure;
plot_sc(U_est_c1,model_true);
% title('Estimated density','interpreter', 'latex');
set(gca,'fontsize',20)
set(gcf,'paperpositionmode','auto')
set(gcf,'position',[100  100 res*1.2 res])

figure;
plot_sc(U_est_c2,model_true);
% title('Estimated density','interpreter', 'latex');
set(gca,'fontsize',20)
set(gcf,'paperpositionmode','auto')
set(gcf,'position',[100  100 res*1.2 res])
end

