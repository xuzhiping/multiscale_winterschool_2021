function [W]=CalcWork(rnnew,rnold,links,Stress);

[fp_inc]=calcplasticincrement(rnnew,rnold,links);
W=sum(sum(Stress.*fp_inc));


function [fp_inc]=calcplasticincrement(rnnew,rn,links);
seg= rn(links(:,2),1:3)-rn(links(:,1),1:3);
segnew=rnnew(links(:,2),1:3)-rnnew(links(:,1),1:3);
dx1=rnnew(links(:,2),1:3)-rn(links(:,1),1:3);
dA=cross(segnew+seg,dx1);
fp_inc=0.5.*(links(:,3:5)'*dA)
