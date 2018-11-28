function [t_nt] = newton_opt(df,ddf)
  %Dynamics of magnetization in primary rotating coordinate system

  t = [5e-4,1e-3,5e-3,1e-2,2e-2,4e-2,6e-2,8e-2,1e-1,1.5e-1]; %Initial Value
  N_max = 100;
  t_nt = zeros(size(t));

  for i = 1:size(t,2)
    t_old = t(i);

    for n = 1:N_max
      t_new = t_old - df(t_old)/ddf(t_old);
      if abs(t_new - t_old) < 1e-8  %Convergence determination
        break
      end
      t_old = t_new;
    end

    if n==N_max
      t_nt(i) = 0; %Convergence
    else
      if t_old < 0
        t_nt(i) = 0; %Inappropriate
      elseif df(t_old - 1e-4)>0 && df(t_old + 1e-4)<0
        t_nt(i) = t_old; %Maximal Value
      else
        t_nt(i) = 0;  %Minimum value
      end
    end
  end
end
