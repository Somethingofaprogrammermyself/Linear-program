
load lin_prog122; 
n = 4;






[xmax1 ii1] = min(b(1:n)./a1(1:n));
[xmax2 ii2] = min(b(n+1:2*n)./a1(n+1:2*n)); ii2 = ii2 +n;
xmax = ( b(ii1)*a2(ii2)-b(ii2)*a2(ii1) ) / ( a1(ii1)*a2(ii2)-a1(ii2)*a2(ii1) );
[xmin1 ii1] = max(b(2*n+1:3*n)./a1(2*n+1:3*n)); ii1 = ii1 + 2*n;
[xmin2 ii2] = max(b(3*n+1:4*n)./a1(3*n+1:4*n)); ii2 = ii2 + 3*n;
xmin = ( b(ii1)*a2(ii2)-b(ii2)*a2(ii1) ) / ( a1(ii1)*a2(ii2)-a1(ii2)*a2(ii1) );

xmax = ceil(xmax);
xmin = floor(xmin);


[ymax1 ii1] = min(b(1:n)./a2(1:n)); 
[ymax2 ii2] = min(b(2*n+1:3*n)./a2(2*n+1:3*n)); ii2 = ii2 + 2*n;
ymax = ( b(ii1)*a1(ii2)-b(ii2)*a1(ii1) ) / ( a2(ii1)*a1(ii2)-a2(ii2)*a1(ii1) );
[ymin1 ii1] = max(b(n+1:2*n)./a2(n+1:2*n)); ii1 = ii1 + n;
[ymin2 ii2] = max(b(3*n+1:4*n)./a2(3*n+1:4*n)); ii2 = ii2 + 3*n;
ymin = ( b(ii1)*a1(ii2)-b(ii2)*a1(ii1) ) / ( a2(ii1)*a1(ii2)-a2(ii2)*a1(ii1) );

ymax = ceil(ymax);
ymin = floor(ymin);


plot([xmin xmax],[0 0],'k')
hold on
plot([0 0],[ymin ymax],'k')
x = [xmin:0.01:xmax]';
for k = 1:4*n,
    plot(x,(b(k)-a1(k)*x)/a2(k));
end
for k = 1:n,
    fill([xmin xmax xmax],[(b(k)-a1(k)*xmin)/a2(k) (b(k)-a1(k)*xmax)/a2(k)  ymax],'r')
    fill([(b(k)-a2(k)*ymin)/a1(k) (b(k)-a2(k)*ymax)/a1(k)  xmax],[ymin ymax ymax],'r')
end
for k = n+1:2*n,
    fill([xmin xmax xmax],[(b(k)-a1(k)*xmin)/a2(k) (b(k)-a1(k)*xmax)/a2(k)  ymin],'g')
    fill([(b(k)-a2(k)*ymin)/a1(k) (b(k)-a2(k)*ymax)/a1(k)  xmax],[ymin ymax ymin],'g')
end
for k = 2*n+1:3*n,
    fill([xmin xmax xmin],[(b(k)-a1(k)*xmin)/a2(k) (b(k)-a1(k)*xmax)/a2(k)  ymax],'b')
    fill([(b(k)-a2(k)*ymin)/a1(k) (b(k)-a2(k)*ymax)/a1(k)  xmin],[ymin ymax ymax],'b')
end
for k = 3*n+1:4*n,
    fill([xmin xmax xmin],[(b(k)-a1(k)*xmin)/a2(k) (b(k)-a1(k)*xmax)/a2(k)  ymin],'y')
    fill([(b(k)-a2(k)*ymin)/a1(k) (b(k)-a2(k)*ymax)/a1(k)  xmin],[ymin ymax ymin],'y')
end


axis([xmin xmax ymin ymax])
t = xlabel('x');
set(t,'fontsize',24)
t = ylabel('y');
set(t,'fontsize',24)


x = [0;0];
mu_ip = 1000;

while mu_ip > 0.000001,
    mu_ip = 0.9 * mu_ip;
    
    for k = 1:10,
       
        sum= 0 ;
        for i=1:16
            Wang = [a1(i)^2 a1(i)*a2(i);a1(i)*a2(i) a2(i)^2];
            sum = sum + (1/((b(i)-a1(i)*x(1)-a2(i)*x(2))^2)*Wang);
        end
        Hessian = mu_ip*sum;
        sum1= 0;
        for i=1:16
            HEUH = [a1(i);a2(i)]; 
            sum1= sum1 + (1/((b(i)-a1(i)*x(1)-a2(i)*x(2))))*HEUH;
        end
        grad= f + mu_ip*sum1;
            
        
        direction = -Hessian \ grad;
        x1 = x + direction;
      
        while max((a1*x1(1) + a2*x1(2)-b)) <= 0 ,
            direction = 2 * direction;
            x1 = x + direction;
        end
       
        while max((a1*x1(1) + a2*x1(2)-b)) > 0,
            direction = direction * 0.9;
            x1 = x + direction;
        end
      
        tau = 2/(1+sqrt(5));
        w1 = (1-tau) * x + tau * x1;
        w = tau * x + (1-tau) * x1;
        while max(abs(x1-x)) > 0.00000000001,
           
            
            sum2= 0;
            for i=1:16
                sum2 = sum2 + log(b(i)-a1(i)*w(1)-a2(i)*w(2));
            end
            J_w = f(1)*w(1) + f(2)*w(2) - mu_ip*sum2;
            
             sum3= 0;
            for i=1:16
                sum3 = sum3 + log(b(i)-a1(i)*w1(1)-a2(i)*w1(2));
            end
            J_w1 = f(1)*w1(1) + f(2)*w1(2) - mu_ip*sum3; 
           
            
            if J_w < J_w1,
                x1 = w1;
                w1 = w;
                w = tau * x + (1-tau) * x1;
            else
                x = w;
                w = w1;
                w1 = (1-tau) * x + tau * x1;
            end
        end
        
    end 
    p = plot(x(1),x(2),'x');
    hold on
end 
set(p,'linewidth',3)
set(p,'markersize',15)
if x(2)>=0,
    t = text(-0.3,-0.1,'initial value');
else
    t = text(-0.3,0.1,'initial value');
end

