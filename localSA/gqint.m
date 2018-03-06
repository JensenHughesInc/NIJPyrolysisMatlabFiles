function yint = gqint( x0, y0 )
%Gaussian quadrature integration for the data from a to b
%   Find the Gaussian point order n = 2, three points method
a = min(x0);
b = max(x0);

xg = zeros(3,1);
ag = zeros(3,1);
fg = zeros(3,1);

xg(1) = -sqrt(3/5);
xg(2) = 0.0;
xg(3) = -xg(1);

ag(1) = 5/9;
ag(2) = 2 - 2*ag(1);
ag(3) = ag(1);

tg = (b-a)/2*xg + (b+a)/2;
bg = (b-a)/2*ag;

[m, n] = size(x0);
if(n ~= 1)
    k = dsearchn(x0', tg);
else
    k = dsearchn(x0, tg);
end


for i=1:3
    if((x0(k(i))-tg(i))>0)
        fg(i) = (y0(k(i))-y0(k(i)-1))/(x0(k(i))-x0(k(i)-1))*(tg(i)-x0(k(i)-1))+y0(k(i)-1);
    elseif((x0(k(i))-tg(i))<0)
        fg(i) = (y0(k(i)+1)-y0(k(i)))/(x0(k(i)+1)-x0(k(i)))*(tg(i)-x0(k(i)))+y0(k(i));
    else
        fg(i) = y0(k(i));
    end
end

% Compute Gaussian quadrature integration
yint = bg'*fg;

end

