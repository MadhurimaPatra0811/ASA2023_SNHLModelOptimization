% clearvars aa yy 
% aa=all_rt(:,:,4);
% yy=nanmean(aa);
% 
% a_st=min(yy)
% b_st=2*(max(yy)-a_st)/pi

aa=all_rt(:,:,5);
x=[0:5:100];
x=repmat(x,21,1);
x=x(:);
y=aa(:);
xx=0:5:100;
yy=nanmean(aa);
a_st=min(yy)
% b_st=2*(max(yy)-a_st)/pi
b_st=max(yy)-a_st