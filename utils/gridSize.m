function [rows, cols] = gridSize(elements)
    divs = divisors(elements);
    rows = divs(ceil(length(divs)/2));
    cols = elements/rows;
    
    % recursively search for better solution
    if (rows == 1 && cols > 10) || (rows > 10 && cols ==1)
        [rows, cols] = gridSize(elements+1);
    end
end