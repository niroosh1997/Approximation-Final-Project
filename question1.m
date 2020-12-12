[num,txt,raw]=xlsread('COVID-19-geographic-disbtribution-worldwide');

countries=raw(:,7);
countries=unique(countries);
countries(strcmp(countries, 'countriesAndTerritories')) = [];


yIsraelCor=[];
xIsraelCor=[];
yUkCor=[];
xUkCor=[];
yUsaCor=[];
xUsaCor=[];

k=1;
j=1;
z=1;
for i=1:size(raw,1)
    if strcmp(raw(i,7),'Israel')
        yIsraelCor(k)=cell2mat(raw(i,5));
        xIsraelCor(k)=k;
        k=k+1;
    end
    if strcmp(raw(i,7),'United_Kingdom')
        yUkCor(j)=cell2mat(raw(i,5));
        yUkCor(j)=j;
        j=j+1;
    end
    
    if strcmp(raw(i,7),'United_States_of_America')
        yUsaCor(z)=cell2mat(raw(i,5));
        xUsaCor(z)=z;
        z=z+1;
    end
    
    
end

yIsraelCor=flip(yIsraelCor);
yUkCor=flip(yUkCor);
yUsaCor=flip(yUsaCor);

% figure
% %israelCor=smooth(israelCor,7);
% plot(israelCor);
% figure
% %ukCor=smooth(ukCor,7);
% plot(ukCor);
% figure
% %usaCor=smooth(usaCor,7);
% plot(usaCor);

y=yIsraelCor;
x=xIsraelCor;
N=size(y,2);

%delete
%  x=1:1:N;
%  y=log(1+x);


%---------- according to deteratives off a,b and compare gradient to 0.   https://en.wikipedia.org/wiki/Linear_least_squares#Example ------------
figure
plot(y);
a_close=(N*sum(x.*y)-sum(x)*sum(y))/(N*sum(x.^2)-sum(x)^2);
b_close=(sum(y)-a_close*sum(x))/N;

hold on
plot(f_source_line(1:1:N,a_close,b_close));



%---------- according to linear algebra. ----------------------------
%https://he.wikipedia.org/wiki/%D7%A9%D7%99%D7%98%D7%AA_%D7%94%D7%A8%D7%99%D7%91%D7%95%D7%A2%D7%99%D7%9D_%D7%94%D7%A4%D7%97%D7%95%D7%AA%D7%99%D7%9D
% we find mclorn function

newColumn=ones(N,1);
X=[newColumn];
k=63; %kelet
for i=1:k
    temp=x.^i;
    X=[temp' X];
end

Ans=((transpose(X)*X)\(transpose(X)))*y';
Ans=flip(Ans);


k=k+1;
hold on
output=f_source_mcl(1:1:N,Ans,k);
plot(output);


%------------------ pade approximation using mclorn series(Ans)---------
%https://math.stackexchange.com/questions/860293/how-to-compute-the-pade-approximation

%order of P and Q
n=floor(k/2); %kelet
m=k-n;

%delete
%Ans=[0,1, -1/2,1/3,-1/4];

A=specialRow(n+m+1,m+1);

for i=1:m
    tempRow1=specialRow(n+m+1,i)*(-1);
    tempRow2=zeros(1,n+m+1);
    for j=1:min([i (n+1)])
        tempRow2=tempRow2+specialRow(n+m+1,m+j)*Ans(i-j+1);
    end
    tempRow1=tempRow1+tempRow2;
    A=[A;tempRow1];
    
end



for i=1:n
    tempRow2=zeros(1,n+m+1);
    for j=1:n+1
        tempRow2=tempRow2+specialRow(n+m+1,m+j)*Ans(m+i-j+1);
    end
    A=[A;tempRow2];
    
end

b=specialRow(n+m+1,1);
b=b';


x=linsolve(A,b);


figure
plot(y);

P=f_source_mcl(1:1:250,x(1:m),m);
Q=f_source_mcl(1:1:250,x(m+1:m+n),n);
output=P./Q;
hold on
plot(output);
output=f_source_mcl(1:1:250,Ans,k);
plot(output);












