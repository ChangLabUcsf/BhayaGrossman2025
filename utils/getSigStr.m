function [str, sig] = getSigStr(P, type, numcomps)

if nargin<2, type = 1; end

str = '';
switch type
    case 1
        if P<0.001
            str = '***';
        elseif P<0.01
            str = '**';
        elseif P<0.05
            str = '*';
        else
            str = 'n.s.';
        end
    case 2
        if P<0.001
            str = 'p<.001';
        elseif P<0.01
            str = ['p=' num2str(P, 2)];
        elseif P<0.05
            str = ['p=' num2str(P, 2)];
        else            
            str = ['p=' num2str(P, 2)];
        end
    case 3 % bonferroni corrected threshold for multiple comparisons
        if P<(0.001/numcomps)
            str = '***';
        elseif P<(0.01/numcomps)
            str = '**';
        elseif P<(0.05/numcomps)
            str = '*';
        else
            str = 'n.s.';
        end
end
sig = P<0.05;
end