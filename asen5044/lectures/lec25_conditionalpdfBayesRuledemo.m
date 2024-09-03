%lec25_conditionalpdfBayesRuledemo.m

clc,clear,close all

%rng(995)

%% 
%%%joint pdf that is marginally uniform in x and conditionally uniform in y
%%% with y bounds depending quadratically on x
dx=0.05;
xrange = 0:dx:14;
yrange = 0:dx:13^2;

[X,Y]=meshgrid(xrange,yrange);

XY = [X(:),Y(:)];

pX = pdf('unif',X(:),1,11);
pY_X = pdf('unif',Y(:),X(:),(X(:)+1).^2);

pXYjoint = pX.*pY_X;
figure(),
surf(X,Y,reshape(pXYjoint,size(X)),'EdgeColor','none')
xlim([0 12])
ylim([0 12^2])
xlabel('x')
ylabel('y')
title('Joint pdf p(x,y) = U_x[1,11] U_{y|x}[x,(x+1)^2]')

%%show prior for x only
pxPrior = pdf('unif',xrange,1,11);
figure(),
plot(xrange,pxPrior,'r','LineWidth',3)
axis([0 14 0 0.2])
title('Prior p(x) = U[1,11]')
xlabel('x')
ylabel('p(x), prior')


%%% Draw sample from uniform p(x)
xtrue = 8.6%14*rand

%%% Draw 3 independent observations from p(y|xtrue)
yobs = xtrue + ((xtrue+1)^2 - xtrue)*rand(3,1)

%%find likelihood for first observation
yobs1 = yobs(1);
py1_x = pdf('unif',yobs1,xrange,(xrange+1).^2);
figure(),
plot(xrange,py1_x,'b','LineWidth',3)
xlim([0 14])
title(['Likelihood p(y_1=',num2str(yobs1),'|x) = U[x,(x+1)^2]'])
xlabel('x')
ylabel(['p(y_1=',num2str(yobs1),'|x), observation likelihood'])

%%show posterior for x after y1
pxPosty1_num = pxPrior.*py1_x;
pxPosty1_den = sum(pxPosty1_num)*dx;
pxPosty1 = pxPosty1_num./pxPosty1_den; 
figure(),
plot(xrange,pxPosty1,'m','LineWidth',3)
xlabel('x')
ylabel(['p(x|y_1=',num2str(yobs1),'), posterior'])
title(['Bayes'' Posterior p(x|y_1=',num2str(yobs1),')'])
%%compute conditional mean given y1
xhat_y1 = sum(xrange.*pxPosty1)*dx
varxhat_y1 = sum(xrange.^2 .* pxPosty1)*dx - xhat_y1^2


%%find likelihood for second observation
yobs2 = yobs(2);
py2_x = pdf('unif',yobs2,xrange,(xrange+1).^2);
figure(),
plot(xrange,py2_x,'b','LineWidth',3)
xlim([0 14])
title(['Likelihood p(y_2=',num2str(yobs2),'|x) = U[x,(x+1)^2]'])
xlabel('x')
ylabel(['p(y_2=',num2str(yobs1),'|x), observation likelihood'])

%%show posterior for x after y2
pxPrior2 = pxPosty1;
pxPosty2_num = pxPrior2.*py2_x;
pxPosty2_den = sum(pxPosty2_num)*dx;
pxPosty2 = pxPosty2_num./pxPosty2_den; 
figure(),
plot(xrange,pxPosty2,'m','LineWidth',3)
xlabel('x')
ylabel(['p(x|y_1,y_2=',num2str(yobs2),'), posterior'])
title(['Bayes'' Posterior p(x|y_1,y_2=',num2str(yobs2),')'])
%%compute conditional mean given y1 and y2
xhat_y12 = sum(xrange.*pxPosty2)*dx
varxhat_y12 = sum(xrange.^2 .* pxPosty2)*dx - xhat_y12^2


%%find likelihood for third observation
yobs3 = yobs(3);
py3_x = pdf('unif',yobs3,xrange,(xrange+1).^2);
figure(),
plot(xrange,py3_x,'b','LineWidth',3)
xlim([0 14])
title(['Likelihood p(y_2=',num2str(yobs3),'|x) = U[x,(x+1)^2]'])
xlabel('x')
ylabel(['p(y_2=',num2str(yobs3),'|x), observation likelihood'])

%%show posterior for x after y3
pxPrior3 = pxPosty2;
pxPosty3_num = pxPrior3.*py3_x;
pxPosty3_den = sum(pxPosty3_num)*dx;
pxPosty3 = pxPosty3_num./pxPosty3_den; 
figure(),
plot(xrange,pxPosty3,'m','LineWidth',3)
xlabel('x')
ylabel(['p(x|y_1,y_2,y_3=',num2str(yobs3),'), posterior'])
title(['Bayes'' Posterior p(x|y_1,y_2,y_3=',num2str(yobs3),')'])
%%compute conditional mean given y1 and y2 and y3
xhat_y123 = sum(xrange.*pxPosty3)*dx
varxhat_y3 = sum(xrange.^2 .* pxPosty3)*dx - xhat_y123^2