clear all;
global bigT;
global A;
global D;
global c1;
global c2;
global c3;
global b;
global v;
global y;
global u;
global steps;
global s;
global maxsize;
global maxvector;

global optstrat;
global optstratd;
global maxnstrats;

alctr = 0; %alctr counts loops through al
for al = 0.5:1.5:3.5 %al varies a parameter of interest (in this case, y)
    clear fit;
    clear fitd;
    clear fitq;
    clear fitdq;

    alctr = alctr+1;
maxnstrats = 2;

A = 0.05;
bigT =25;
z = 0.5;
b = 2/3;
c1 = 0.25;
c3 = 0.02;
y = 2;
u = 0.8;
s = 1;
steps = 2;
alpha = 0.5;
maxsize = 100;
v = 0.05;
maxvector = maxsize./steps;
ld = 0.9;
ls = 0.5;
mn_md = 0.5;
D = 0.00323;

Xs = linspace(steps,maxsize,maxvector);
Xd = linspace(steps,maxsize,maxvector);

   y = al;
   c2 = z.*c1;
     if u == 0
         u = 0.01;
     end
     
Dd = zeros(length(Xs),length(Xd));

fitd = ones(length(Xs),length(Xd),bigT);
Xdmat = repmat(Xd,length(Xs),1);
fitd(:,:,bigT) = Xdmat;


n = 0;
breaker = 0;
while breaker == 0;
    n = n+1
    if n > 1
        recoptstratd(:,:,:) = optstratd;
        recoptstrat(:,:,:) = optstrat;
    end
    
for xx = 1:length(Xs)
    for yy = 1:length(Xd)
   if yy < xx     
fit(xx,yy,bigT) = fitd(yy,xx,bigT);
   else
       fit(xx,yy,bigT) = alpha.*Xdmat(yy,xx);
   end
    end
end

for t = bigT-1:-1:1
    if n > 1
        Dd(:,:) = (optstratd(:,:,t)-1)./(maxnstrats-1);
    end
    fitq = fit(:,:,t+1);
    for littlex = 1:length(Xs)
        for littley = 1:length(Xd)
            if littlex > littley
                fitq(littlex,littley) = fitd(littley,littlex,t);
            end
        end
    end
    for X = 1:length(Xs)
        x = Xs(X);
        for Y = 1:length(Xd)
            yy = Xd(Y);
            r = x./yy;
            ctr = 0;
            for initialfights = 0:1/(maxnstrats-1):1
                ctr = ctr+1;
                pfights = initialfights;
    pfightd = Dd(X,Y);
    pswins = 1./2+1./2.*tanh(y.*(r-1));
    p(1) = u.*(1-pfights).*(1-pfightd)+(1-u);
    p(2) = u.*pfights.*pfightd.*(1-pswins).*(1-v);
    p(3) = u.*pfights.*pfightd.*pswins;
    p(4) = u.*(1-pfights).*pfightd;
    p(5) = u.*pfights.*(1-pfightd);
    p(6) = u.*pfights.*pfightd.*(1-pswins).*v;
    
    newx(1) = x+s.*A.*x.^b;
    newx(2) = x+s.*A.*(1-c2-c1-c3).*x.^b;
    newx(3) = x+s.*A.*(1-c1-c3).*x.^b;
    newx(4) = x+s.*A.*(1-c1).*x.^b;
    newx(5) = x+s.*A.*(1-c3).*x.^b;
    newx(6) = newx(2);
    
    newx = newx-s.*D.*x;

    
    newy(1) = yy+s.*A.*yy.^b;
    newy(2) = yy+s.*A.*(1-c1-c3).*yy.^b;
    newy(3) = yy+s.*A.*(1-c1-c2-c3).*yy.^b;
    newy(4) = yy+s.*A.*(1-c3).*yy.^b;
    newy(5) = yy+s.*A.*(1-c1).*yy.^b;
    newy(6) = newy(2);
    newy = newy -s.*D.*yy;
    
    newx = (newx < x).*x+(newx >= x).*newx;
    newy = (newy < yy).*yy+(newy >= yy).*newy;
    
    newxcat = newx./steps;
    newxfloor = floor(newxcat);
    newxprob = newxcat-newxfloor;
    newxceil = newxfloor+1;
    newxfloor = (newxfloor >= maxvector).*maxvector+(newxfloor < maxvector).*newxfloor;
    newxceil = (newxceil >= maxvector).*maxvector+(newxceil < maxvector).*newxceil;
    
    newycat = newy./steps;
    newyfloor = floor(newycat);
    newyprob = newycat-newyfloor;
    newyceil = newyfloor+1;
    newyfloor = (newyfloor >= maxvector).*maxvector+(newyfloor < maxvector).*newyfloor;
    newyceil = (newyceil >= maxvector).*maxvector+(newyceil < maxvector).*newyceil;
    
    
   for zz = 1:5
    efit(zz) = (1-newxprob(zz)).*(1-newyprob(zz)).*fitq(newxfloor(zz),newyfloor(zz))+(1-newxprob(zz)).*newyprob(zz).*fitq(newxfloor(zz),newyceil(zz))...
        + newxprob(zz).*(1-newyprob(zz)).*fitq(newxceil(zz),newyfloor(zz))+newxprob(zz).*newyprob(zz).*fitq(newxceil(zz),newyceil(zz));
   end
  efit(6) = ls.*((1-newxprob(6)).*fitd(round(newxfloor(6).*mn_md),newxfloor(6),t+1)+newxprob(6).*fitd(round(newxceil(6).*mn_md),newxceil(6),t+1));
  
    newfit(X,Y,ctr) = p*efit';
   
    
            end
        
        end
    end
    [As,B] = max(newfit,[],3);
    fit(:,:,t) = As;
    optstrat(:,:,t) = B;
    


end

for x = 1:length(Xs)
    for yy = 1:length(Xd)
   if yy >= x     
fitd(x,yy,bigT) = Xdmat(x,yy);
   else
       fitd(x,yy,bigT) = fit(yy,x,bigT);
   end
    end
end

%t = bigT-1;
for t = bigT-1:-1:1
    fitdq = fitd(:,:,t+1);
    for littlex = 1:length(Xs)
        for littley = 1:length(Xd)
            if littlex > littley
                fitdq(littlex,littley) = fit(littley,littlex,t+1);
            end
        end
    end
    for X = 1:length(Xs)
        x = Xs(X);
        for Y = 1:length(Xd)
            yy = Xd(Y);
            r = x./yy;
            ctr = 0;
            for initialfights = 0:1/(maxnstrats-1):1
                ctr = ctr+1;
    pfightd = initialfights;
    pfights = (optstrat(X,Y,t)-1)./(maxnstrats-1);
    pswins = 1./2+1./2.*tanh(y.*(r-1));
    p(1) = u.*(1-pfights).*(1-pfightd)+(1-u);
    p(2) = u.*pfights.*pfightd.*(1-pswins).*(1-v);
    p(3) = u.*pfights.*pfightd.*pswins;
    p(4) = u.*(1-pfights).*pfightd;
    p(5) = u.*pfights.*(1-pfightd);
    p(6) = u.*pfights.*pfightd.*(1-pswins).*v;
    
    newx(1) = x+s.*A.*x.^b;
    newx(2) = x+s.*A.*(1-c2-c1-c3).*x.^b;
    newx(3) = x+s.*A.*(1-c1-c3).*x.^b;
    newx(4) = x+s.*A.*(1-c1).*x.^b;
    newx(5) = x+s.*A.*(1-c3).*x.^b;
    newx(6) = newx(2);
    
    newx = newx-s.*D.*x;

    
    newy(1) = yy+s.*A.*yy.^b;
    newy(2) = yy+s.*A.*(1-c1-c3).*yy.^b;
    newy(3) = yy+s.*A.*(1-c1-c2-c3).*yy.^b;
    newy(4) = yy+s.*A.*(1-c3).*yy.^b;
    newy(5) = yy+s.*A.*(1-c1).*yy.^b;
    newy(6) = newy(2);
    newy = newy -s.*D.*yy;
    
    newx = (newx < x).*x+(newx >= x).*newx;
    newy = (newy < yy).*yy+(newy >= yy).*newy;
    
    newxcat = newx./steps;
    newxfloor = floor(newxcat);
    newxprob = newxcat-newxfloor;
    newxceil = newxfloor+1;
    newxfloor = (newxfloor >= maxvector).*maxvector+(newxfloor < maxvector).*newxfloor;
    newxceil = (newxceil >= maxvector).*maxvector+(newxceil < maxvector).*newxceil;
    
    newycat = newy./steps;
    newyfloor = floor(newycat);
    newyprob = newycat-newyfloor;
    newyceil = newyfloor+1;
    newyfloor = (newyfloor >= maxvector).*maxvector+(newyfloor < maxvector).*newyfloor;
    newyceil = (newyceil >= maxvector).*maxvector+(newyceil < maxvector).*newyceil;
    
    
    
   for zz = 1:5
    efit(zz) = (1-newxprob(zz)).*(1-newyprob(zz)).*fitdq(newxfloor(zz),newyfloor(zz))+(1-newxprob(zz)).*newyprob(zz).*fitdq(newxfloor(zz),newyceil(zz))...
        + newxprob(zz).*(1-newyprob(zz)).*fitdq(newxceil(zz),newyfloor(zz))+newxprob(zz).*newyprob(zz).*fitdq(newxceil(zz),newyceil(zz));
   end
   efit(6) = ld.*((1-newyprob(6)).*fitd(round(newxfloor(6).*mn_md),newxfloor(6),t+1)+newyprob(6).*fitd(round(newxceil(6).*mn_md),newyceil(6),t+1));
    newfitd(X,Y,ctr) = p*efit';
   
    
            end
        
        end
    end
    [Ad,Bd] = max(newfitd,[],3);
    fitd(:,:,t) = Ad;
    optstratd(:,:,t) = Bd;
    
end
if n > 1
for ii = 1:length(Xs)
    for jj = 1:length(Xd)
        if ii < jj
            if ii > jj/2
            oldstrat(ii,jj,:) = recoptstrat(ii,jj,:);
            newstrat(ii,jj,:) = optstrat(ii,jj,:);
            oldstratd(ii,jj,:) = recoptstratd(ii,jj,:);
            newstratd(ii,jj,:) = optstratd(ii,jj,:);
            end
        end
    end
end
end

if n > 1

doptstratd(n) = sum(sum(sum(abs(oldstratd-newstratd))));
doptstrat(n) = sum(sum(sum(abs(oldstrat-newstrat))));
end
if n > 1
if doptstratd(n) == 0
    if doptstrat(n) ==0;
        breaker = 1;
    end
end
if n == 50
    breaker = 1;
end

end

end
BB(:,:) = optstrat(:,:,20);
            BBB(:,:) = optstratd(:,:,20);
for ii = 1:length(Xs);
    for jj = 1:length(Xd)
        if ii < jj
            
        substrat(ii,jj) = (BB(ii,jj)-1)./(maxnstrats-1);
        domstrat(ii,jj) = (BBB(ii,jj)-1)./(maxnstrats-1);
        else
            substrat(ii,jj) = 0;
            domstrat(ii,jj) = 0;
        end
    end
end

 recsubstrat(:,:,alctr) = substrat;
 recdomstrat(:,:,alctr) = domstrat;
 nrec(alctr) = n;

[fevict, afights, histo1, histo2] = ForwardSimulation(optstrat,optstratd);
 
 finaldist(:,alctr) = histo1;
 finaldistrand(:,alctr) = histo2;
 finalfights(:,1:2,alctr) = afights';
 alrec(alctr) = al;
 finalevict(:,alctr) = fevict;

end
