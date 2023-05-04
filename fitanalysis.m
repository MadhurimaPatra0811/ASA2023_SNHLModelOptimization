bc=coeffs(:,2).*coeffs(:,3);
a=coeffs(:,1);
b=coeffs(:,2);
c=coeffs(:,3);
d=coeffs(:,4);
k=0.01;
dy_ln = 2.*sqrt((bc-k)./k.*(c.^2)); %length of dynamic range
ll = d-dy_ln./2;
ul = d+dy_ln./2;
figure
plot(ul-ll);

% 
% 
dbs = 0:0.1:100;
y_hat=[]; %fitted curve
for ii=1:17
    aa=a(ii)+b(ii).*atan(c(ii).*(dbs-d(ii)));
    y_hat=[y_hat ; aa];
end


df = diff(y_hat,1,2);
sat_rate = a+b.*(pi/2);


ii=19;
figure
plot([0:5:100], nanmean(all_rt(:,:,ii)))
hold on
plot(dbs(1:end-1),df(20-ii,:),'r')


df2=diff(df,1,2);


% dbs = 0:0.1:100;
% y_hat=[]; %fitted curve
% for ii=1:17
%     %aa=a(ii)+b(ii).*atan(c(ii).*(dbs-d(ii)));
%     aa=a(ii)+b(ii)./(1+exp(-(dbs-c(ii))./d(ii)));
%     y_hat=[y_hat ; aa];
% end


% df = diff(y_hat,1,2);
% sat_rate = a+b;



