function analysis1(x1,x2,x3,x4,t,equ,i)
    % response information
    x1_infos = lsiminfo(x1(:,1),t,0); x1_infophi = lsiminfo(x1(:,3),t,0);
    x2_infos = lsiminfo(x2(:,1),t,0); x2_infophi = lsiminfo(x2(:,3),t,0);
    x3_infos = lsiminfo(x3(:,1),t,0); x3_infophi = lsiminfo(x3(:,3),t,0);
    x4_infos = lsiminfo(x4(:,1),t,0); x4_infophi = lsiminfo(x4(:,3),t,0);

    figure(i)
    subplot(2,2,1); yyaxis left; plot(t,x1(:,1)); hold on; plot(t,equ,'k--','LineWidth',0.05);
    title('Reponse of y(t) to Initial Condition: [-0.5 0 0 0]'); ylabel('y_1: s(t) (m)');
    p=find(x1(:,1)==max(x1(:,1))); 
    text(t(p),x1(p,1),['(',num2str(t(p)),',',num2str(x1(p,1)),')'],'color','b');
    p=find(x1(:,1)==min(x1(:,1))); 
    text(t(p),x1(p,1),['(',num2str(t(p)),',',num2str(x1(p,1)),')'],'color','b');
    text(x1_infos.SettlingTime,0,(num2str(x1_infos.SettlingTime)),'color','b');
    axis([0 max(t) -inf inf],'auto y');
    
    yyaxis right; plot(t,x1(:,3)); xlabel('Time (seconds)'); ylabel('y_2: phi(t) (rad)');
    hold on; plot(t,equ,'k--','LineWidth',0.05); % axis([-inf inf -pi/4 pi/4])
    p=find(x1(:,3)==max(x1(:,3))); 
    text(t(p),x1(p,3),['(',num2str(t(p)),',',num2str(x1(p,3)),')'],'color','r');
    p=find(x1(:,3)==min(x1(:,3))); 
    text(t(p),x1(p,3),['(',num2str(t(p)),',',num2str(x1(p,3)),')'],'color','r');
    text(x1_infophi.SettlingTime,0,(num2str(x1_infophi.SettlingTime)),'color','r');
    axis([0 max(t) -inf inf],'auto y');
    
    subplot(2,2,2); yyaxis left; plot(t,x2(:,1)); hold on; plot(t,equ,'k--','LineWidth',0.05);
    title('Reponse of y(t) to Initial Condition: [0 -0.5 0 0]'); ylabel('y_1: s(t) (m)');
    p=find(x2(:,1)==max(x2(:,1))); 
    text(t(p),x2(p,1),['(',num2str(t(p)),',',num2str(x2(p,1)),')'],'color','b');
    p=find(x2(:,1)==min(x2(:,1))); 
    text(t(p),x2(p,1),['(',num2str(t(p)),',',num2str(x2(p,1)),')'],'color','b');
    text(x2_infos.SettlingTime,0,(num2str(x2_infos.SettlingTime)),'color','b');
    axis([0 max(t) -inf inf],'auto y');
    
    yyaxis right; plot(t,x2(:,3)); xlabel('Time (seconds)'); ylabel('y_2: phi(t) (rad)');
    hold on; plot(t,equ,'k--','LineWidth',0.05);
    p=find(x2(:,3)==max(x2(:,3))); 
    text(t(p),x2(p,3),['(',num2str(t(p)),',',num2str(x2(p,3)),')'],'color','r');
    p=find(x2(:,3)==min(x2(:,3))); 
    text(t(p),x2(p,3),['(',num2str(t(p)),',',num2str(x2(p,3)),')'],'color','r');
    text(x2_infophi.SettlingTime,0,(num2str(x2_infophi.SettlingTime)),'color','r');
    axis([0 max(t) -inf inf],'auto y');
    
    subplot(2,2,3); yyaxis left; plot(t,x3(:,1)); hold on; plot(t,equ,'k--','LineWidth',0.05);
    title('Reponse of y(t) to Initial Condition: [0 0 -0.7 0]'); ylabel('y_1: s(t) (m)');
    p=find(x3(:,1)==max(x3(:,1))); 
    text(t(p),x3(p,1),['(',num2str(t(p)),',',num2str(x3(p,1)),')'],'color','b');
    p=find(x3(:,1)==min(x3(:,1))); 
    text(t(p),x3(p,1),['(',num2str(t(p)),',',num2str(x3(p,1)),')'],'color','b');
    text(x3_infos.SettlingTime,0,(num2str(x3_infos.SettlingTime)),'color','b');
    axis([0 max(t) -inf inf],'auto y');
    
    yyaxis right; plot(t,x3(:,3)); xlabel('Time (seconds)'); ylabel('y_2: phi(t) (rad)');
    hold on; plot(t,equ,'k--','LineWidth',0.05);
    p=find(x3(:,3)==max(x3(:,3))); 
    text(t(p),x3(p,3),['(',num2str(t(p)),',',num2str(x3(p,3)),')'],'color','r');
    p=find(x3(:,3)==min(x3(:,3))); 
    text(t(p),x3(p,3),['(',num2str(t(p)),',',num2str(x3(p,3)),')'],'color','r');
    text(x3_infophi.SettlingTime,0,(num2str(x3_infophi.SettlingTime)),'color','r');
    axis([0 max(t) -inf inf],'auto y');
    
    subplot(2,2,4); yyaxis left; plot(t,x4(:,1)); hold on; plot(t,equ,'k--','LineWidth',0.05);
    title('Reponse of y(t) to Initial Condition: [0 0 0 -0.5]'); ylabel('y_1: s(t) (m)');
    p=find(x4(:,1)==max(x4(:,1))); 
    text(t(p),x4(p,1),['(',num2str(t(p)),',',num2str(x4(p,1)),')'],'color','b');
    p=find(x4(:,1)==min(x4(:,1))); 
    text(t(p),x4(p,1),['(',num2str(t(p)),',',num2str(x4(p,1)),')'],'color','b');
    text(x4_infos.SettlingTime,0,(num2str(x4_infos.SettlingTime)),'color','b');
    axis([0 max(t) -inf inf],'auto y');
    
    yyaxis right; plot(t,x4(:,3)); xlabel('Time (seconds)'); ylabel('y_2: phi(t) (rad)');
    hold on; plot(t,equ,'k--','LineWidth',0.05);
    p=find(x4(:,3)==max(x4(:,3))); 
    text(t(p),x4(p,3),['(',num2str(t(p)),',',num2str(x4(p,3)),')'],'color','r');
    p=find(x4(:,3)==min(x4(:,3))); 
    text(t(p),x4(p,3),['(',num2str(t(p)),',',num2str(x4(p,3)),')'],'color','r');
    text(x4_infophi.SettlingTime,0,(num2str(x4_infophi.SettlingTime)),'color','r');
    axis([0 max(t) -inf inf],'auto y');
end