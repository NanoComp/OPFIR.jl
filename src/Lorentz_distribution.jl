function lorentz_dist(ν, Δ_f_NT, f_pump)
    return  1/π .* Δ_f_NT ./ ((ν - f_pump).^2 + Δ_f_NT^2);
end

function f_NT_ampl(ν, Δ_f_NT, f_pump)
    SHB = Δ_f_NT^2 ./ ((ν - f_pump).^2 + Δ_f_NT^2)
    SHB = SHB/sum(SHB)
    return SHB
end
