
function [ret]=predicate(XTrain,XTest)
[m2,n2]=size(XTest);%XTest=56*1024
ret=zeros(m2,1);
test=zeros(1,n2);
for i=1:1:m2
    test(1,:)=XTest(i,:);
    ret(i)=classifier(XTrain,test);
end

end

function [num]=classifier(XTrain,test)
kind=7;
[m1,n1]=size(XTrain);%m1=140,n1=1024
[m2,n2]=size(test);%m2=56,n2=1024

weight=testing(XTrain,test);
cla=zeros(kind,1);
nu=zeros(kind,1);

for m=1:1:7
    train=zeros(floor(m1/kind),n1);%20*1024
    w=zeros(floor(m1/kind),1);%20*1
    for i=1:1:(floor(m1/kind))
        train(i,:)=XTrain(i+(m-1)*(m1/kind),:);
        w(i)=weight(i+(m-1)*(m1/kind));
    end
    train=train';
    test=test';
    product=train*w;
    sum2=0;
    for i=1:1:m2
		I1=test(i);
        I2=product(i);
		I_MINUS=(I1-I2)^2;
		sum2=sum2+I_MINUS;
    end;
    cla(m)=sum2^(1/2);
    nu(m)=(m-1);
end
[cla,nu]=BubbleSort(cla,nu,kind);
num=nu(1);  
nu(1)
end

function [cla,nu]=BubbleSort(cla,nu,kind)
for i=1:1:kind
    for j=2:1:(kind+1-i)
        if cla(j-1)>cla(j)
            temp=cla(j-1);
            cla(j-1)=cla(j);
            cla(j)=temp;
            
            tem=nu(j-1);
            nu(j-1)=nu(j);
            nu(j)=tem;
        end
    end
end  
end

function [p]=testing(Xtrain,XTest)
R=1;
error=0.0001;
a=0;
[m1,n1]=size(Xtrain);
p=ones(m1,1);

while 0<1
    if a<=20
        a=a+1;
        tra=Xtrain';
        tea=XTest';
        [perish,p]=conjugate_g_test(tra,tea,R,p);
        if perish<error
            break;
        end
        R=R*2;
    end
    if(a>20)
        break;
    end
end
end