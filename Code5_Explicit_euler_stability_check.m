function y= Explicit_euler()
a=0;%initial time
b=100*60;% final time which was 100 minutes converted into secs
N=4;%number of points in our case unlike other times we fixed N in our case
h=-1*(a-b)/(N);% calculating step size from N ,a&b


t=zeros(N+1,1);%initializing value of t as array of zeros of dimensions 12*1
y=zeros(N+1,1);%initializing value of solutions as array of zeros of dimensions 12*1
t(1)=a;%setting value of starting t equal to a
y(1)=1;%fixing value of hieght at t=o as 1;
%euler's explicit method
for i=1:N 
    y(i+1)=y(i)+h*getf(y(i));%implementing the formula for euler's explicit method
    t(i+1)=t(i)+h;%increasing the value of time for each iteration
end

 for j=1:N+1
exact(j)=getexact(t(j));
RE(j)=abs(exact(j)-y(j))/exact(j);
 end
 tiledlayout(2,1)
 nexttile
 plot(t,y,"b")%ploting explicit solutions vs time
hold on;
 plot(t,exact,"g"); %ploting exact solution vs time
 xlabel("Time");%defining x axis of the plot
ylabel('height of chemical in reactor left');% defining y axis of the plot
title('Plot of height of chemical in reactor left vs Time ');% defining title of the plot
 legend('Explicit','Exact solution');%labeling the graphs
nexttile
plot(t,RE)%ploting relative error or deviation from solution vs time
 xlabel("Time");
ylabel('Deviation from exact solution');

end

function exact=getexact(t) % defining function to get exact value
d=2/1000;%constants 
D=250/1000;
g=9.8;
exact=((1)^(1/2)-((d/D)^(2))*t*(g*(1/2))^(1/2))^2;%expression for exact solution
end

function out=getf(y)%defining the function to be used in explicit method
d=2/1000;%constants 
D=250/1000;
g=9.8;
out=-1*(d/D)^2*(2*y*g)^(1/2);
end