function [fevict, afights, histo1, histo2] = ForwardSimulation(optstrat,optstratd)
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

for i = 1:10000
    breaker = 0;
while breaker < 2
    domvec(i) = nrmrnd(60,5,1,1);
    submult = nrmrnd(0.9,0.1,1,1);
    subvec(i) = domvec(i).*submult;
    if (domvec(i) > 0)
        if (domvec(i) < 100)
            if subvec(i) < domvec(i)
                if subvec(i) > 0
                breaker = 2;
                end
            end
        end
    end
end
end

rvec = subvec./domvec;
for xn = 1:2
subvecrec(:,1,xn) = subvec;
domvecrec(:,1,xn) = domvec;
end
for xn = 1:2
for t = 1:bigT-1;
    if t == 1
        subvec = subvecrec(:,1,xn);
        domvec = domvecrec(:,1,xn);
    end
    subdown = floor(subvec./steps);
subdown= (subdown<maxvector).*subdown+(subdown >= maxvector).*maxvector;    
    subup = subdown+1;
    subup = (subup<maxvector).*subup+(subup >=maxvector).*maxvector;
    subp = subvec-subdown;
    subsize(1,:) = subdown;
    subsize(2,:) = subup;
    
    domdown = floor(domvec./steps);
    domdown = (domdown<maxvector).*domdown+(domdown>=maxvector).*maxvector;
    domup = domdown+1;
    domup = (domup<maxvector).*domup+(domup >=maxvector).*maxvector;
    domp = domvec-domdown;
    domsize(1,:) = domdown;
    domsize(2,:) = domup;
    
  
  for i = 1:length(subvec)
      if rand()> subp
          xx = 1;
      else
          xx = 2;
      end
      if rand() > domp
          yy = 1;
      else
          yy = 2;
      end
        ssize = subvec(i);
        dsize = domvec(i);
        
            oss = optstrat(subsize(xx,i),domsize(yy,i),t);
            osd = optstratd(subsize(xx,i), domsize(yy,i),t);
            pfightsub(i) = (oss-1)./(maxnstrats - 1);
            pfightdom(i) = (osd-1)./(maxnstrats - 1);
             pswinssub(i) = 1./2+1./2.*tanh(y.*(rvec(i)-1));
             
             if xn == 2
                 isafight(i) = 0;
             elseif rand() < u
                 isafight(i) = 1;
             else
                 isafight(i) = 0;
             end
                 
              if rand()<pfightsub(i)
                  pfights(i) = 1;
              else
                  pfights(i) = 0;
              end
              if rand()< pfightdom(i)
                  pfightd(i) = 1;
              else
                  pfightd(i) = 0;
              end
              if rand()< pswinssub(i)
                  pswins(i) = 1;
              else
                  pswins(i) = 0;
              end
              if rand() < v
                  evicted(i) = 1;
              else
                  evicted(i) = 0;
              end
              
           if ssize > 1
    p(1) = isafight(i).*(1-pfights(i)).*(1-pfightd(i))+(1-isafight(i));
      
    p(2) = isafight(i).*pfights(i).*pfightd(i).*(1-pswins(i)).*(1-evicted(i));
    p(3) = isafight(i).*pfights(i).*pfightd(i).*pswins(i);
    p(4) = isafight(i).*(1-pfights(i)).*pfightd(i);
    p(5) = isafight(i).*pfights(i).*(1-pfightd(i));
        p(6) = isafight(i).*pfights(i).*pfightd(i).*(1-pswins(i)).*evicted(i);
           else
               p(1:5) = 0;
               p(6) = 1;
           end
        x = ssize;
    newx(1) = x+s.*A.*x.^b;
    newx(2) = x+s.*A.*(1-c2-c1-c3).*x.^b;
    newx(3) = x+s.*A.*(1-c1-c3).*x.^b;
    newx(4) = x+s.*A.*(1-c1).*x.^b;
    newx(5) = x+s.*A.*(1-c3).*x.^b;
    newx(6) = newx(2);
    newx = newx-s.*D.*x;

   yy = dsize; 
    newy(1) = yy+s.*A.*yy.^b;
    newy(2) = yy+s.*A.*(1-c1-c3).*yy.^b;
    newy(3) = yy+s.*A.*(1-c1-c2-c3).*yy.^b;
    newy(4) = yy+s.*A.*(1-c3).*yy.^b;
    newy(5) = yy+s.*A.*(1-c1).*yy.^b;
    newy(6) = newy(2);
    newy = newy -s.*D.*yy;
    
    newx = (newx < x).*x+(newx >= x).*newx;
    newy = (newy < yy).*yy+(newy >= yy).*newy;
    
    newx(6) = 1;
    newy(6) = 100;
    
    newx1(i) = 0;
    newy1(i) = 0;
    for zz = 1:6
        newx1(i) = newx1(i) + p(zz).*newx(zz);
        newy1(i) = newy1(i) + p(zz).*newy(zz);
    end
     prec(i,t,xn) = pfights(i).*pfightd(i);
  
  end
  

  for i = 1:length(subvec)
  if newx1(i) < newy1(i)
  subvec(i) = newx1(i);
  domvec(i) = newy1(i);
  else
      subvec(i) = newy1(i);
      domvec(i) = newx1(i);
  end
  end
    subvecrec(:,t+1,xn) = subvec;
  domvecrec(:,t+1,xn) = domvec;
  recvec = subvec./domvec;
  
 
end
end

freqfights = prec;
avfights = mean(freqfights,1);
afights(1,:) = avfights(1,:,1);
afights(2,:) = avfights(1,:,2);

rvecrec = subvecrec./domvecrec;
ed2 = linspace(0.5,1,50);
fevict = sum(rvecrec(:,bigT,1)<0.02);
figure(4)
h1 = histogram(rvecrec(:,bigT,1),ed2);
hold on
h2 = histogram(rvecrec(:,bigT,2),ed2);
hold off
h1.Normalization = 'probability';
h1.BinWidth = 0.01;
h2.Normalization = 'probability';
h2.BinWidth = 0.01;
xlabel('Size ratio, subordinate / dominant','FontSize',14);
ylabel('Frequency','FontSize',14);
axis('square');

histo1 = h1.Data;
histo2 = h2.Data;

end

