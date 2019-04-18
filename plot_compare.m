function fig = plot_compare(n,U_true,U_est,model)

fig = subplot(2,1,1);
plot(model.x(2:1:end),U_true{1}(1,2:1:end),'--','color',[0.8,0,0],'linewidth',2)
hold on;
plot(model.x(2:end),U_true{n}(1,2:end),'-','color',[.6,.6,.6],'linewidth',3)
plot(model.x(2:end),U_est{n}(1,2:end),'-','color',[0.8,0.61,0],'linewidth',3)
h = legend('Initial state','True Model','Approximate Model');

% title(sprintf('Small vehicle density at t=%d',n),'interpreter', 'latex')
title(sprintf('timestep=%d',n),'interpreter', 'latex');
axis([0 model.len 0 1.5])
% set(h,'Location','NorthWest','interpreter', 'latex',...
%     'position',[0.11 0.80 0.5 0.14])
    set(gca,'linewidth',1.2,'FontName','Times')
    set(gca,'fontsize',30)    
    set(gca,'xtick',[])   
set(gca,'position',[.14 .53 .84 .38])
    ylabel('$\rho_1$','interpreter', 'latex') 

subplot(2,1,2),
plot(model.x(2:1:end),U_true{1}(2,2:1:end),'--','color',[0.8,0,0],'linewidth',2)
hold on
plot(model.x(2:end),U_true{n}(2,2:end),'-','color',[0.6,0.6,.6],'linewidth',3)

plot(model.x(2:end),U_est{n}(2,2:end),'-','color',[0.8,0.61,0],'linewidth',3)
% h = legend('Initial state','True Model','Estimated Model');
% title(sprintf('Large vehicle density at t=%d',n),'interpreter', 'latex')
% set(h,'Location','NorthWest','interpreter', 'latex',...
%     'position',[0.11 0.33 0.5 0.14])
set(gca,'linewidth',1.2)
    set(gca,'fontsize',30)
    axis([0 model.len 0 1.4])
ylabel('$\rho_2$','interpreter', 'latex')
xlabel('space','interpreter', 'latex')
  xh = get(gca,'xlabel'); % handle to the label object
    p = get(xh,'position'); % get the current position property
    p(2) = p(2)*0.7;       % negative values put the label below the axis
    set(xh,'position',p)   % set the new position
%     set(gca,'position',[.09 .07 .87 .41])
  set(gca,'position',[.14 .15 .84 .30],'FontName','Times')
  res = 600;
    set(gcf,'paperpositionmode','auto')
    set(gcf,'position',[10  100 res res*0.70])
  
hold off


%-----------------------------------
% fig = subplot(2,1,1);
% plot(model.x(2:1:end),U_true{1}(1,2:1:end),'--','color',[0.8,0,0],'linewidth',2)
% ylim([0,1]);
% ylabel('$\rho_1$','interpreter', 'latex','Fontsize', 22)
% xlabel('x','interpreter', 'latex','Fontsize', 22)
% title('Small vehicle density at k=0','interpreter', 'latex','Fontsize', 22);
% set(gca,'fontsize',22)
% 
% subplot(2,1,2),
% plot(model.x(2:1:end),U_true{1}(2,2:1:end),'--','color',[0.8,0,0],'linewidth',2)
% ylim([0,1]);
% ylabel('$\rho_2$','interpreter', 'latex','Fontsize', 22)
% xlabel('x','interpreter', 'latex','Fontsize', 22)
% title('Large vehicle density at k=0','interpreter', 'latex','Fontsize', 22);
% set(gca,'fontsize',22)
end

