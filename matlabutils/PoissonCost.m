function c = PoissonCost(m,img,B,b)

I = m*B+b;
I(I<b) = b;
c = sum(I - img.*log(I));

end