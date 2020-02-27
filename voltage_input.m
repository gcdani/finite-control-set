function output = voltage_input( voltage, w, t, theta, option )
    wt = w * t;
    
    if option == 0
        vas = sqrt(2/3) * voltage * cos(wt);
        vbs = sqrt(2/3) * voltage * cos(wt - 2 * pi / 3);
        vcs = sqrt(2/3) * voltage * cos(wt + 2 * pi / 3);
    
        Ks = [  cos(theta) cos(theta - 2 * pi / 3) cos(theta + 2 * pi / 3);
            sin(theta) sin(theta - 2 * pi / 3) sin(theta + 2 * pi / 3)];
        Ks = (2/3) .* Ks;
    
        vqd = Ks * [vas;vbs;vcs];
        vqd(3) = 0;
        vqd(4) = 0;
    
    elseif option == 1
       f = w / (2 * pi );
       time22pi = 1 / f;
       if t > 0
           K = ceil( t / time22pi );
       else
           K = 1;
       end
%        if wt > 2 * pi
%            wt
%        end
       
       wt = wt - (K - 1) * 2 * pi;
       voltage = voltage * sqrt(3 / 2);
        
       if wt >= - pi / 6 && wt <= pi / 6
           mode = 1;
       elseif wt > pi / 6 && wt <= pi / 2
           mode = 2;
       elseif wt > pi / 2 && wt <= 5 * pi / 6
           mode = 3;
       elseif wt > 5 * pi / 6 && wt <= 7 * pi / 6
           mode = 4;
       elseif wt > 7 * pi / 6 && wt <= 3 * pi / 2
           mode = 5;
       elseif wt > 3 * pi / 2 && wt <= 11 * pi / 6
           mode = 6;
       else
           mode = 1;
       end
       
       vqds = (2 / 3) * voltage * exp(1i * (mode - 1) * pi / 3);
       vqd = [real(vqds); -imag(vqds); 0; 0];
       
    end
    
    output = vqd;

end

