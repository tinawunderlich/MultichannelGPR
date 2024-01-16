function r=empcorrcoeff(A,B)

% r=empcorrcoeff(A,B)
% 
% empirical correlation coefficient
% Dr. Tina Wunderlich, CAU Kiel, 2020, tina.wunderlich@ifg.uni-kiel.de

r=sum((A-mean(A)).*(B-mean(B)))/sqrt(sum((A-mean(A)).^2)*sum((B-mean(B)).^2));