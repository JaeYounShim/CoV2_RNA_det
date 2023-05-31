%%
clear Tms_mtx_nonzeros
for i=1:size(Tms_mtx,2)
    Tms_mtx_nonzeros(1,i)=mean(nonzeros(Tms_mtx(:,i)));
end

plot(Tms_x,Tms_mtx,'r*')
hold on 
plot((1:1:size(Tms_x,2)),Tms_mtx_nonzeros,'ko','MarkerSize',12)
ylabel('Temp(C)')
xlabel('Length(bps)')
xlim([1.5 size(Tms_x,2)+0.5])
ylim([0 90])

%% 
errorbar(T(:,1),T(:,2),T(:,3),'-*k')
hold on 
errorbar(T(:,1),T(:,4),T(:,5),'-*b')

ylabel('T (C)')
xlabel('Length of probes (bps)')
xlim([5.5 20.5])
% ylim([0 10])

legend('JYS','PNAbio')

%%
% bar([1 2],[mean(T(:,1)),mean(T(:,2))],'FaceColor',[0.3 0.3 0.3])
% hold on
scatter(T_x(:,1),T(:,1),'*k','SizeData',100,'LineWidth',.5)
hold on 
scatter(T_x(:,2),T(:,2),'*b','SizeData',100,'LineWidth',.5)
hold on 
errorbar(1,mean(T(:,1)),std(T(:,1)),'or','LineWidth',2)
hold on 
errorbar(1.5,mean(T(:,2)),std(T(:,2)),'or','LineWidth',2)

ylabel('T (C)')
xlim([0.7 1.8])
xticks([1 1.5])
xticklabels({'T (JYS)','T (PNAbio)'})
box off

%%
for i=1:size(T,1)
plot([T_x(i,1),T_x(i,2)],[T(i,1),T(i,2)],'-o',LineWidth=1,Color=[0.5 0.5 0.5])
hold on
end
hold on
scatter(T_x(:,1),T(:,1),'*k','SizeData',100,'LineWidth',.5)
hold on 
scatter(T_x(:,2),T(:,2),'*b','SizeData',100,'LineWidth',.5)
hold on 
errorbar(1,mean(T(:,1)),std(T(:,1)),'or','LineWidth',2)
hold on 
errorbar(1.5,mean(T(:,2)),std(T(:,2)),'or','LineWidth',2)

ylabel('T (C)')
xlim([0.7 1.8])
% ylim([24 44])

xticks([1 1.5])
xticklabels({'T (JYS)','T (PNAbio)'})
box off
