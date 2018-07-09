function d = bhattaDist(p,q)
    d = -log(BC(p,q));
end

function bc = BC(p,q)
    bc = sum( p .* q);
end