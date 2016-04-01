function [WM] = image_to_network(I,pad,MAXVAL)
% function [Ncut] = graphcuts(I)
% Input: I image
% pad: spatial connectivity; eg. 3
% MAXVAL: maximum image value
% Output: Ncut: Binary map 0 or 1 corresponding to image segmentation
I = double(I); [H,W] = size(I);
% Find weights between nodes I1 and I2, w = exp(a*abs(I1-I2));
% Set a to have a weight of 0.01 for diff = MAXVAL
a = log(0.01)/MAXVAL;
x = [0:MAXVAL/100:MAXVAL]';
y = exp(a*x);

ws = 2*pad + 1;
if(ws <= 3)
    ws = 3;
end

%Build the weight matrix
WM = zeros(H*W,H*W); countWM = 0;
for kk = 1:W
    for jj = 1:H
        mask = zeros(H,W);
        cs = kk-pad; ce = kk+pad; rs = jj-pad; re = jj+pad;
        if(cs<1)
            cs = 1;
        end
        if(ce>W)
            ce = W;
        end
        if(rs<1)
            rs = 1;
        end
        if(re>H)
            re = H;
        end
        mask(rs:re,cs:ce) = 1;
        idx = find(mask==1);
        p = abs(I(idx) - I(jj,kk)); p = exp(a*p);
        countWM = countWM + 1;
        WM(countWM,idx) = p(:)';
    end
end

% Weight between a node and iteself is 0
WM(1:(H*W+1):(H*W)^2)=0;

