function check_sum(p, sol)
    sub_sum = zeros(p.num_layers)
    ri = p.r_int
    for i = 1:p.num_layers
        index1 = (18*p.num_freq + p.n_vib) * (i-1) + 1
        index2 = (18*p.num_freq + p.n_vib) * i
        sub_sum[i] = sum(sol[index1:index2]) * ri[i]
    end
    return sum((sub_sum))
end

function pump_dist(p, sol)
    pump_layer = zeros(p.num_layers)
    for i in 1:p.num_layers
        pump_layer[i] = pump_total(p, sol, i)
    end
    return pump_layer
end
