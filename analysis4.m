function analysis4(x1,x2,x3,t1,t2,t3,equ1,equ2,equ3,i)
    % response information
    x1_infos = lsiminfo(x1(:,1),t1,0); x1_infophi = lsiminfo(x1(:,3),t1,0);
    x2_infos = lsiminfo(x2(:,1),t2,0); x2_infophi = lsiminfo(x2(:,3),t2,0);
    x3_infos = lsiminfo(x3(:,1),t3,0); x3_infophi = lsiminfo(x3(:,3),t3,0);

    figure(i)
    subplot(3,1,1); yyaxis left; stairs(t1,x1(:,1)); hold on; stairs(t1,equ1,'k--','LineWidth',0.05);
    title('K'); ylabel('y_1: s(t) (m)');
    p=find(x1(:,1)==max(x1(:,1))); 
    text(t1(p),x1(p,1),['(',num2str(t1(p)),',',num2str(x1(p,1)),')'],'color','b');
    p=find(x1(:,1)==min(x1(:,1))); 
    text(t1(p),x1(p,1),['(',num2str(t1(p)),',',num2str(x1(p,1)),')'],'color','b');
    text(x1_infos.SettlingTime,0,(num2str(x1_infos.SettlingTime)),'color','b');
    axis([0 max(t1) -inf inf],'auto y');
    
    yyaxis right; stairs(t1,x1(:,3)); xlabel('Time (seconds)'); ylabel('y_2: phi(t) (rad)');
    hold on; stairs(t1,equ1,'k--','LineWidth',0.05); % axis([-inf inf -pi/4 pi/4])
    p=find(x1(:,3)==max(x1(:,3))); 
    text(t1(p),x1(p,3),['(',num2str(t1(p)),',',num2str(x1(p,3)),')'],'color','r');
    p=find(x1(:,3)==min(x1(:,3))); 
    text(t1(p),x1(p,3),['(',num2str(t1(p)),',',num2str(x1(p,3)),')'],'color','r');
    text(x1_infophi.SettlingTime,0,(num2str(x1_infophi.SettlingTime)),'color','r');
    axis([0 max(t1) -inf inf],'auto y');
    
    subplot(3,1,2); yyaxis left; stairs(t2,x2(:,1)); hold on; stairs(t2,equ2,'k--','LineWidth',0.05);
    title('K_d'); ylabel('y_1: s(t) (m)');
    p=find(x2(:,1)==max(x2(:,1))); 
    text(t2(p),x2(p,1),['(',num2str(t2(p)),',',num2str(x2(p,1)),')'],'color','b');
    p=find(x2(:,1)==min(x2(:,1))); 
    text(t2(p),x2(p,1),['(',num2str(t2(p)),',',num2str(x2(p,1)),')'],'color','b');
    text(x2_infos.SettlingTime,0,(num2str(x2_infos.SettlingTime)),'color','b');
    axis([0 max(t2) -inf inf],'auto y');
    
    yyaxis right; stairs(t2,x2(:,3)); xlabel('Time (seconds)'); ylabel('y_2: phi(t) (rad)');
    hold on; stairs(t2,equ2,'k--','LineWidth',0.05);
    p=find(x2(:,3)==max(x2(:,3))); 
    text(t2(p),x2(p,3),['(',num2str(t2(p)),',',num2str(x2(p,3)),')'],'color','r');
    p=find(x2(:,3)==min(x2(:,3))); 
    text(t2(p),x2(p,3),['(',num2str(t2(p)),',',num2str(x2(p,3)),')'],'color','r');
    text(x2_infophi.SettlingTime,0,(num2str(x2_infophi.SettlingTime)),'color','r');
    axis([0 max(t2) -inf inf],'auto y');

    subplot(3,1,3); yyaxis left; stairs(t3,x3(:,1)); hold on; stairs(t3,equ3,'k--','LineWidth',0.05);
    title('K_d^*'); ylabel('y_1: s(t) (m)');
    p=find(x3(:,1)==max(x3(:,1))); 
%     text(t3(p),x3(p,1),['(',num2str(t3(p)),',',num2str(x3(p,1)),')'],'color','b');
    p=find(x3(:,1)==min(x3(:,1))); 
    text(t3(p),x3(p,1),['(',num2str(t3(p)),',',num2str(x3(p,1)),')'],'color','b');
    text(x3_infos.SettlingTime,0,(num2str(x3_infos.SettlingTime)),'color','b');
    axis([0 max(t3) -inf inf],'auto y');
    
    yyaxis right; stairs(t3,x3(:,3)); xlabel('Time (seconds)'); ylabel('y_2: phi(t) (rad)');
    hold on; % stairs(t3,equ3,'k--','LineWidth',0.05);
    p=find(x3(:,3)==max(x3(:,3))); 
    text(t3(p),x3(p,3),['(',num2str(t3(p)),',',num2str(x3(p,3)),')'],'color','r');
    p=find(x3(:,3)==min(x3(:,3))); 
    text(t3(p),x3(p,3),['(',num2str(t3(p)),',',num2str(x3(p,3)),')'],'color','r');
    text(x3_infophi.SettlingTime,0,(num2str(x3_infophi.SettlingTime)),'color','r');
    axis([0 max(t3) -inf inf],'auto y');
end