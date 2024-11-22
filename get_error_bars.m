function wsi = get_error_bars(jumpers,total)                        %finds error bars
        x = jumpers;
        n = total;
        alpha = .05;
        phat =  (x./n);
        z=sqrt(2).*erfcinv(alpha);
        den=1+(z^2./n);xc=(phat+(z^2)./(2*n))./den;
        halfwidth=(z*sqrt((phat.*(1-phat)./n)+(z^2./(4*(n.^2)))))./den;
        wsi=[xc(:) xc(:)]+[-halfwidth(:) halfwidth(:)];
%        wsi = wsi .* 100;
    end
