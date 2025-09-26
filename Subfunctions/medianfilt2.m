function A=medianfilt2(A,n)

% 2d-median filter (Tina Wunderlich, CAU Kiel, September 2025)

if length(n)==1
    n=[n n];
end
nr=(n(1)-1)/2; % half width for rows
nc=(n(2)-1)/2; % half width for columns

[r,c]=size(A);

% padding with zeros
B=[zeros(nr,c+2*nc); zeros(r,nc) A zeros(r,nc); zeros(nr,c+2*nc)];

for i=nc+1:nc+c
    for j=nr+1:nr+r
        A(j-nr,i-nc)=median(B(j-nr:j+nr,i-nc:i+nc),'all','omitnan');
    end
end