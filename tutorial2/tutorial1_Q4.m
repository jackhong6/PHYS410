for ii = 1:200
    x(ii) = ii/100 -1;
    y(ii) = ii/200;
end

for ii = 1:200
    for jj = 1:200
        z(ii,jj) = x(ii)^2 - y(jj)^2;
    end
end

x2 = [1:200]/100 - 1;
y2 = [1:200]/200;

[xx,yy] = ndgrid(x2,y2);
z2 = xx.^2 - yy.^2;

surf(xx,yy,z2)
