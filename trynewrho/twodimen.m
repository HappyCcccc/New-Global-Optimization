

clear all; close all; clc;



%z =0.0;

% Compute the value of the quadratic function.

 

d = 2;

n = 5000;

Corner = zeros(1,d);% matrix for corner node

Width = zeros(1,d);% matrix for width and length

Center = zeros(2*n+1,d);%matrix for center node

area = ones(2*n+1,1); % matrix for retrangle area

f = zeros(1,1); % function value for center points

getrho = zeros(1,1);

a = zeros(1,d);

V = zeros(1,n);

M = zeros(1,n);

A = zeros(n,1);

B = zeros(n,1);

 

for i = 1:d

    Width(1,i)=1.0;

    Center(1,i)=0.5;

end

 

maxWidth = 0.0;

maxIndex = 1;

maxD = 0;

%MinArea = 1.0;

Center(maxIndex,:);

 
%A=zeros(1,n);
for i = 1:n

    

maxWidth = 0.0;

% split the retrangle

    %decide which axis to split on

for j = 1:d

        if(Width(maxIndex,j)>maxWidth)

            maxWidth = Width(maxIndex,j);

            maxD = j;

      %      maxIndex = j;

        end

end

 

%update corner node and width/length information for new retrangles

%NC1 = zeros(d);

 

Width(maxIndex, maxD)=1/3*maxWidth;

 

NC1 = zeros(1,d);

NC2 = zeros(1,d);

NW1 = zeros(1,d); 

NW2 = zeros(1,d);  

 

for j = 1:d

    if (j == maxD)

        NC1(j) = 1/3*maxWidth + Corner(maxIndex,j);

        NC2(j) = 2/3*maxWidth + Corner(maxIndex,j);

        NW1(j) = 1/3*maxWidth;

        NW2(j) = 1/3*maxWidth;

    else 

        NC1(j) = Corner(maxIndex,j);

        NC2(j) =  Corner(maxIndex,j);

        NW1(j) = Width(maxIndex,j);

        NW2(j) =  Width(maxIndex,j);

    end

end

 

%update corner and width information

Corner = [Corner;NC1;NC2];

Width = [Width;NW1;NW2];

[m,d] = size(Corner);



%calculate center nodes and retrangle area

Center = Corner + Width/2;

 

%{

for j = 1:m

    for k = 1:d

        Center(j,k)=Corner(j,k)+1/2*Width(j,k);

%        area(j,1)=area(j,1)*Width(j,k);

%        f(j) = funct(Center);

    end

end

%}

 

%get the smallest function value among all vertices

Mn = intmax('int64');

for j = 1:m

   % for k = 1:d

   %     a(k) = Center(j,k);

   % end

   a = Center(j,:);

        f(j) = funct1(a,d);

      

%         if((f(j)-0.5*Width(j,d)*abs(dfunct(a,d)))<Mn)
% 
%             Mn = f(j)-0.5*Width(j,d)*abs(dfunct(a,d));
% 
%         end

    

end

Mn = min(f);

M(1,i)=Mn;

% A(i,1) = (log(2*i+1)/(2*i+1)) * log(1/abs(Mn));
%A(i,1) = Mn;

B(i,1) = i;

% get the volume of smallest rectangle:travsal all rectangles

area = ones(2*n+1,1);

for j = 1:m

    for k = 1:d

        area(j,1)= area(j,1)*Width(j,k);

    end

end

 

 

Vn = intmax('int64');

for j = 1 : m

        if (area(j,1) < Vn)

            Vn = area(j,1);

        end

end

Vn;

V(1,i)=Vn;

 
%A(i,1) = (log(2*i+1)/(2*i+1)) * log(1/abs(Vn));
A(i,1) = log(1/abs(Vn));

% get GVn
GVn = (Vn * log(n)).^(2/d);
%GVn = (V(1,i)*log(n)).^(2/d);

%GVn =  (area(i,1)*log(n)).^(2/d);

%choose the retrangle which have the largest Rho

rho = intmin;

%figure (1)

%axis equal 

%axis ([0 1 0 1])

%hold on

for j = 1 : m

    a = Center(j,:);

      % getrho(j) = area(j,1).^(2/d)/(funct(a,d)-M(1,i)+d*GVn);
      getrho(j) = 1;
      grad = dfunct1(a,d);
      
       for k = 1:d
        getrho(j) = getrho(j) * (2*area(j,1)/(sqrt(max(0,f(j)+ grad(k,1)*0.5*Width(j,k)-Mn+GVn))+sqrt(max(0,f(j)- grad(k,1) *0.5*Width(j,k)-Mn)+GVn)));
       end
        if (getrho(j) > rho)

            rho = getrho(j);

            maxIndex = j;

      %      figure (1);

      %      plot(Center(maxIndex,1),Center(maxIndex,2),'k.');

      %      hold on

      %      Center(maxIndex,:);

            

      %      drawnow; 

       %     pause(0.00000000000001);

        end

end


 Center(maxIndex, :);
 

 %plot(Center(:,1), Center(:,2),'.');
 
 
 

end


f=f.';

%plot(Center(:,1),f,'.');

plot(Center(:,1), Center(:,2),'.')

%plot(B(:,1),A(:,1),'.');

%plot(A,'.');
 
%c= sqrt(2)/(8*sqrt(3))

 

 

function y = funct(c,d)

%n = size(c);

sum=0;

for i = 1:d

    sum = sum+ (c(i)-sqrt(0.5)).^2;

end

y = sum;

end

function y = dfunct(c,d)

%n = size(c);

sum=0;

for i = 1:d

    sum = sum+ 2*(c(i)-sqrt(0.5));

end

y = sum;

end

 

function y = funct1(c,d)

%n = size(c);

nc=zeros(d);

sum=20;

B=4.0;

for i = 1:d

    nc(i)=2*B*c(i)-B-0.1;

    sum = sum+ nc(i).^2-10*cos(2*pi*nc(i));

end

  y=sum;

%y = sum*0.001;

end

function y = dfunct1(c,d)

nc=zeros(d);

sum=zeros(d,1);

B=4.0;

for i = 1:d

    nc(i)=2*B*c(i)-B-0.1;

    sum(i) = 2*nc(i)+20*pi*sin(2*pi*nc(i));

end

  y=sum;
end

function y = dfunct2(c,d)

nc=zeros(d);

sum=0;

B=4.0;

for i = 1:d

    nc(i)=2*B*c(i)-B-0.1;

    sum = sum+ 2*nc(i)+20*pi*sin(2*pi*nc(i));

end

  y=sum;
end

 

function y = funct2(c,d)

%n = size(c);

nc=zeros(d);

sum=0;

B=5.1/(4*pi.^2);

C=5/pi;R=6;T=1/(8*pi);



for i = 1:d-1

    nc(i)=15*c(i)-5;

    nc(i+1)=15*c(i+1);

    

    sum = (1/51.95)*((nc(i+1)-B*nc(i).^2+C*nc(i)-R).^2+10*(1-T)*cos(nc(i))-44.81);

end

y = sum*0.001;

end

 

