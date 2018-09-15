syms f(x)
f(x,y) = sin(x)*sin(y);
df = diff(f,x)
diff(f,y)

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