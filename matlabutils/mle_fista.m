function xFin = mle_fista(B,img,b,xInit)

for i = 1:size(xInit,1)
    xTmp = xInit;
    if rem(i-1,15) < 6
        xTmp(i) = xTmp(i) + 1;
        Ltmp(i) = norm(gradf(xTmp,B,b,img) - gradf(xInit,B,b,img));
    else
        xTmp(i) = xTmp(i) + 50;
        Ltmp(i) = norm(gradf(xTmp,B,b,img) - gradf(xInit,B,b,img))/50;
    end
end

L = max(Ltmp)/10;

eta = 1.2;
x = xInit;
y = x;
t = 1;

for iter = 1:100
    cgrad = gradf(y,B,b,img);
    i = 1;
    pLyk = pL(y,L,cgrad);
    Lnew = L;
    while any(f(pLyk,B,b,img) > qL(pLyk,y,Lnew,B,b,img))
        Lnew = eta^i*L;
        i = i + 1;
        pLyk = pL(y,Lnew,cgrad);
    end
    xnew = pLyk;
    if norm(xnew-x)/(norm(x)+norm(xnew)) < 5e-3 % change later
        break;
    end
    tnew = (1+sqrt(1+4*t^2))/2;
    y = xnew +(t-1)/tnew*(xnew-x);
    x = xnew;
    t = tnew;
end
xFin = xnew;

end

function c = f(x,B,b,img)
I = B*x+b;
c = sum(I - img.*log(I));
end

function cgrad = gradf(x,B,b,img)
I = B*x+b;
cgrad = B'*(1-img./I);
end

function x = pL(y,L,cgrad)
x = y - cgrad/L;
end

function q = qL(x,y,L,B,b,img)
q = f(y,B,b,img) + dot(x-y,gradf(y,B,b,img))+L*(x-y)'*(x-y)/2;
end
