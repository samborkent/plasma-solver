function [ o ] = myfastint2str( x )
maxvalue=max(x(:));
%maxvalue=intmax(class(x));%Alternative implementation based on class
required_digits=ceil(log(double(maxvalue+1))/log(10));
o=repmat(x(1)*0,size(x,1),required_digits);%initialize array of required size
for c=size(o,2):-1:1
   o(:,c)=mod(x,10);
   x=(x-o(:,c))/10;
end
o=char(o+'0');

if isempty(o) == 1
    o = '0';
end
end