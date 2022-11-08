## Copyright (C) 2022 jmiguel
## Author: jmiguel <jmiguel@jmiguel-nitro>
## Created: 2022-11-04

function [mtz_trayec] = multigrafica (x1, x2, num_exp, exp, datos, final_ruido)
  x_tray = [x1; x2];
  
  ipp = 1;
  ii = 1;
    for i = ipp:ipp + 1 
      mtz_trayec(i, :) = x_tray(ii, :);
      ii = ii + 1;
    endfor
    mtz_trayec;
    
  figure(6);
    subplot(2,2,exp)
      hold on;
      grid on;

      plot( (0:datos), mtz_trayec(1,1:datos+1),"marker", "v", "markerEdgeColor", "k", ... 
           "markersize", 2, "linewidth", 1, "color","r");
           
      plot( (0:datos), mtz_trayec(2,1:datos+1),"marker", "o", "markerEdgeColor", "k", ... 
            "markersize", 2, "linewidth", 1, "color","b");
    if exp>2
      xlabel("Tiempo k")
    end
      yticks([0])
    if (exp == 1) 
      ylabel({"x_{1}, x_{2}"})
    end
    if (exp == 3) 
      ylabel({"x_{1}, x_{2}"})
    end   
      set(0, "defaultlinelinewidth", 2)
      %set(gca, 'FontSize', 17)
      line([final_ruido final_ruido], [0 1], "linestyle", "--", "linewidth", 2,  "color", "black")
      xlim ([0, datos])
%title(['Temperature is ',num2str(exp),' C'])
      title(["Exp: ", num2str(exp)])
    ##  legend("Estado x_{1}","Estado x_{2}")
      
  
endfunction
