function [yout]=ivp_implicit()

a=0;%initial time
b=100*60;% final time which was 100 minutes converted into secs
N=21;%number of points in our case unlike other times we fixed N in our case
h=-1*(a-b)/(N); % calculating step size from N ,a&b
y(1)=1;%fixing value of hieght at t=o as 1;
tol=1e-8;% tolrance value for newtons method at each point/iterations
t=zeros(N+1,1);%initializing value of t as array of zeros of dimensions 12*1
exact=zeros(N+1,1);%initializing value of exact value as array of zeros of dimesions 12*1
t(1)=0;
for i=2:N+1
    t(i)=t(i-1)+h;%making mesh of time with h as step size
    yold=y;%storing value of y in yold for next iteration of newtons method
    f=getf(y,yold,h);
    while abs(f)>tol
        y=y-f/getdf(y,h);%newtons method
        f=getf(y,yold,h);
    end
    yout(1)=1;%cause value of y at t=0 was 1 
    yout(i)=y;%storing %value of y in array  and then storing that in yold for next iterations
  
end
 for j=1:N+1
exact(j)=getexact(t(j));%storing exact value of function in exact variable
RE(j)=abs(exact(j)-yout(j))/exact(j);
 end
 
 tiledlayout(2,1);
 nexttile
 plot(t,yout,"b");%ploting implicit solutions vs time
xlabel("Time");
ylabel('height of chemical in reactor left');
title('Plot of height of chemical in reactor left vs Time ');
hold on
 plot(t,exact,"g");%ploting exact solution vs time
 legend('Implicit','Exact solution')%labeling the graphs
 hold off
 nexttile
  plot(t,RE)%ploting relative error or deviation from solution vs time
 xlabel("Time");
ylabel('Deviation from exact solution');

end

function exact=getexact(t) % defining function to get exact value
d=2/1000;%constants 
D=250/1000;
g=9.8;
exact=((1)^(1/2)-((d/D)^(2))*t*(g*(1/2))^(1/2))^2;%exact solution expressions
end

function f=getf(y,yold,h) %defining function for newtons method and taking f in consideration to check if value is more than tolerance or not
d=2/1000;
D=250/1000;
g=9.8;
f=y-yold+h*(d/D)^2*(2*g*y)^(1/2);
end

function df=getdf(y,h)%defining derivative of function for newtons method
d=2/1000; 
D=250/1000;
g=9.8;
df=1+h*(d/D)^2*(g/2)^(1/2)*(y)^(-1/2);
end

