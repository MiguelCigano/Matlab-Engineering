function [u,wk]=michirp(wini,wfin,N,T)
  for i=0:(N-1),
      wk(i+1)=wini+(i-1)/N*(wfin-wini); 
      w(i+1)=sat(10*i/N)*sat((N-i)/(0.1*N));    
      u(i+1,1)=i*T; 
      u(i+1,2)=w(i+1)*sin(wk(i+1)*i*T);
  end
end