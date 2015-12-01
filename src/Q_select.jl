function Q_selectn_hl(J)
    K = 3;
    return (J^2-K^2)/(J*(2*J+1))
end

function Q_selectn_lh(J)
    K = 3;
    return ((J+1)^2-K^2)/((J+1)*(2*J+1))
end