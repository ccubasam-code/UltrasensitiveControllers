function [] = Func_AddLegend(x,y,j,a,b,scale3,var,Color3)
txt1 = strcat(num2str(scale3(1)),var(j));
txt2 = strcat(num2str(scale3(2)),var(j));
txt3 = strcat(num2str(scale3(3)),var(j));
% text(x,y,[txt1 txt2 txt3],'FontSize',12)
text(x,y+b,txt1,'FontSize',10)
text(x,y,txt2,'FontSize',10)
text(x,y-b,txt3,'FontSize',10)
plot(x-a,y+b,'o','Color',Color3(1,:),'MarkerFaceColor',Color3(1,:),'MarkerSize',6)
plot(x-a,y,'o','Color',Color3(2,:),'MarkerFaceColor',Color3(2,:),'MarkerSize',6)
plot(x-a,y-b,'o','Color',Color3(3,:),'MarkerFaceColor',Color3(3,:),'MarkerSize',6)
hold off
