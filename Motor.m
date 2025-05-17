function FM=Motor(Throt)
if Throt>=0 && Throt<=1
    f=-10.89*Throt^3+21.848*Throt^2-2.0347*Throt+0.037762;
    M=0.08741*Throt^3-0.2267*Throt^2-0.008974*Throt-0.002028;
elseif Throt<0
    f=0.037762;
    M=-0.002028;
else
    f=-10.89*(1)^3+21.848*(1)^2-2.0347*(1)+0.037762;
    M=0.08741*(1)^3-0.2267*(1)^2-0.008974*(1)-0.002028;
end
FM=[f, M];
end