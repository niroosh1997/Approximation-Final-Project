function y=f_source_mcl(x,C,N)
s=size(x,2);
y=zeros(1,s);
for i=1:N
    temp=x.^(i-1);
    y=y+C(i)*temp;
end
