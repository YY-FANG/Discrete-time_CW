function analysis2(s1r,s2r,s3r,s4r,t,equ,i)
    %% response information
    s1r_infos = lsiminfo(s1r(:,1),t,0); s1r_infophi = lsiminfo(s1r(:,3),t,0);
    s2r_infos = lsiminfo(s2r(:,1),t,0); s2r_infophi = lsiminfo(s2r(:,3),t,0);
    s3r_infos = lsiminfo(s3r(:,1),t,0); s3r_infophi = lsiminfo(s3r(:,3),t,0);
    s4r_infos = lsiminfo(s4r(:,1),t,0); s4r_infophi = lsiminfo(s4r(:,3),t,0);
    
    figure(i)
    subplot(2,2,1); yyaxis left; stairs(t,s1r(:,1)); hold on; stairs(t,equ,'k--','LineWidth',0.05);
    title('Reponse of y(t) to Initial Condition: [-0.5 0 0 0]'); ylabel('y_1: s(t) (m)');
    p=find(s1r(:,1)==max(s1r(:,1))); 
    text(t(p),s1r(p,1),['(',num2str(t(p)),',',num2str(s1r(p,1)),')'],'color','b');
    p=find(s1r(:,1)==min(s1r(:,1))); 
    text(t(p),s1r(p,1),['(',num2str(t(p)),',',num2str(s1r(p,1)),')'],'color','b');
    text(s1r_infos.SettlingTime,-0.05,(num2str(s1r_infos.SettlingTime)),'color','b');
    axis([0 max(t) -inf inf],'auto y');
    
    yyaxis right; stairs(t,s1r(:,3)); xlabel('Time (seconds)'); ylabel('y_2: phi(t) (rad)');
    hold on; stairs(t,equ,'k--','LineWidth',0.05); 
    p=find(s1r(:,3)==max(s1r(:,3))); 
    text(t(p),s1r(p,3),['(',num2str(t(p)),',',num2str(s1r(p,3)),')'],'color','r');
    p=find(s1r(:,3)==min(s1r(:,3))); 
    text(t(p),s1r(p,3),['(',num2str(t(p)),',',num2str(s1r(p,3)),')'],'color','r');
    text(s1r_infophi.SettlingTime,0,(num2str(s1r_infophi.SettlingTime)),'color','r');
    axis([0 max(t) -inf inf],'auto y');
    
    subplot(2,2,2); yyaxis left; stairs(t,s2r(:,1)); hold on; stairs(t,equ,'k--','LineWidth',0.05);
    title('Reponse of y(t) to Initial Condition: [0 -0.5 0 0]'); ylabel('y_1: s(t) (m)');
    p=find(s2r(:,1)==max(s2r(:,1))); axis([-inf inf -1.8 0.5])
    text(t(p),s2r(p,1),['(',num2str(t(p)),',',num2str(s2r(p,1)),')'],'color','b');
    p=find(s2r(:,1)==min(s2r(:,1))); 
    text(t(p),s2r(p,1),['(',num2str(t(p)),',',num2str(s2r(p,1)),')'],'color','b');
    text(s2r_infos.SettlingTime,-0.05,(num2str(s2r_infos.SettlingTime)),'color','b');
    axis([0 max(t) -inf inf],'auto y');
    
    yyaxis right; stairs(t,s2r(:,3)); xlabel('Time (seconds)'); ylabel('y_2: phi(t) (rad)');
    hold on; stairs(t,equ,'k--','LineWidth',0.05);
    p=find(s2r(:,3)==max(s2r(:,3))); 
    text(t(p),s2r(p,3),['(',num2str(t(p)),',',num2str(s2r(p,3)),')'],'color','r');
    p=find(s2r(:,3)==min(s2r(:,3))); 
    text(t(p),s2r(p,3),['(',num2str(t(p)),',',num2str(s2r(p,3)),')'],'color','r');
    text(s2r_infophi.SettlingTime,0,(num2str(s2r_infophi.SettlingTime)),'color','r');
    axis([0 max(t) -inf inf],'auto y');
    
    subplot(2,2,3); yyaxis left; stairs(t,s3r(:,1)); hold on; stairs(t,equ,'k--','LineWidth',0.05);
    title('Reponse of y(t) to Initial Condition: [0 0 -0.7 0]'); ylabel('y_1: s(t) (m)');
    p=find(s3r(:,1)==max(s3r(:,1))); 
    text(t(p),s3r(p,1),['(',num2str(t(p)),',',num2str(s3r(p,1)),')'],'color','b');
    p=find(s3r(:,1)==min(s3r(:,1))); 
    text(t(p),s3r(p,1),['(',num2str(t(p)),',',num2str(s3r(p,1)),')'],'color','b');
    text(s3r_infos.SettlingTime,-0.05,(num2str(s3r_infos.SettlingTime)),'color','b');
    axis([0 max(t) -inf inf],'auto y');
    
    yyaxis right; stairs(t,s3r(:,3)); xlabel('Time (seconds)'); ylabel('y_2: phi(t) (rad)');
    hold on; stairs(t,equ,'k--','LineWidth',0.05);
    p=find(s3r(:,3)==max(s3r(:,3))); 
    text(t(p),s3r(p,3),['(',num2str(t(p)),',',num2str(s3r(p,3)),')'],'color','r');
    p=find(s3r(:,3)==min(s3r(:,3))); 
    text(t(p),s3r(p,3),['(',num2str(t(p)),',',num2str(s3r(p,3)),')'],'color','r');
    text(s3r_infophi.SettlingTime,0,(num2str(s3r_infophi.SettlingTime)),'color','r');
    axis([0 max(t) -inf inf],'auto y');
    
    subplot(2,2,4); yyaxis left; stairs(t,s4r(:,1)); hold on; stairs(t,equ,'k--','LineWidth',0.05);
    title('Reponse of y(t) to Initial Condition: [0 0 0 -0.5]'); ylabel('y_1: s(t) (m)');
    p=find(s4r(:,1)==max(s4r(:,1))); axis([-inf inf -1.8 0.5])
    text(t(p),s4r(p,1),['(',num2str(t(p)),',',num2str(s4r(p,1)),')'],'color','b');
    p=find(s4r(:,1)==min(s4r(:,1))); 
    text(t(p),s4r(p,1),['(',num2str(t(p)),',',num2str(s4r(p,1)),')'],'color','b');
    text(s4r_infos.SettlingTime,-0.05,(num2str(s4r_infos.SettlingTime)),'color','b');
    axis([0 max(t) -inf inf],'auto y');
    
    yyaxis right; stairs(t,s4r(:,3)); xlabel('Time (seconds)'); ylabel('y_2: phi(t) (rad)');
    hold on; stairs(t,equ,'k--','LineWidth',0.05);
    p=find(s4r(:,3)==max(s4r(:,3))); 
    text(t(p),s4r(p,3),['(',num2str(t(p)),',',num2str(s4r(p,3)),')'],'color','r');
    p=find(s4r(:,3)==min(s4r(:,3))); 
    text(t(p),s4r(p,3),['(',num2str(t(p)),',',num2str(s4r(p,3)),')'],'color','r');
    text(s4r_infophi.SettlingTime,0,(num2str(s4r_infophi.SettlingTime)),'color','r');
    axis([0 max(t) -inf inf],'auto y');
end