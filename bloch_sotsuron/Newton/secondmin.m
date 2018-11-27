function [tsl_opt] = secondmin(ta)
  %Ask for the second smallest ta(The smallest value is 0)
  %極大値をとる最小のtslを求めたい。別のプログラムで極小値や発散などを0にしているので必然的に0の次に大きい値が求める値へ
  t = sort(ta);

  for i = 1:size(ta,2)
    if min(ta) < t(i)
      tsl_opt = t(i);
      break;
    end
  end
end
