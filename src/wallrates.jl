function fn_1(θ, MFP, r₀, x)
    return 2*MFP^2*sin(θ).*real(acos((r₀^2-MFP^2*sin(θ).^2-x^2)/2/MFP./sin(θ)/x));
end

function α_max(θ, MFP, r₀, x)
    return  real(acos((r₀^2-MFP^2*sin(θ).^2-x^2)/2/MFP./sin(θ)/x));
end

function f1(θ, MFP, r₀, x)
    f0 = (sqrt(((r₀-x)./tan(θ[1])).^2+(-x*cos(θ[2]*α_max(θ[1], MFP, r₀, x))
    + sqrt(x^2*cos(θ[2]*α_max(θ[1], MFP, r₀, x))^2-x^2+r₀^2))^2)).*sin(θ[1]);
    return f0*α_max(θ[1], MFP, r₀, x);
end

function g(θ, MFP, r₀, x)
    return sin(θ[1])*α_max(θ[1], MFP, r₀, x);
end

function fn_2(v, mass, kBT, l)
    return l./v * 4*π*(mass/2/π/kBT)^(3/2).*(v).^2.*exp(-mass*(v).^2/2/kBT);
end

function WallRate(radius, pressure, r_int, ntotal, M, T, NA, v, σ_GKC)
    r₀ = radius/100; # cm to m
    mass = M/1000/NA; # kg per molecule
    τ_GKC = 1/(ntotal*v*σ_GKC*(1e-10)^2); # sec
    MFP = 0.732*T/pressure/σ_GKC/100; # m
    kB = 1.3806488*1e-23;
    kBT = kB * T;

    index = find(r_int .< r₀ - MFP);
    if isempty(index) == true
        index = 0;
    end
# println(MFP)
    surf_fraction = zeros(1,length(r_int));
    θ_min = 0;
    θ_max = 0;

    for k in index[end]+1:length(r_int)
        x = r_int[k]; # m
        fn_0(θ) = fn_1(θ, MFP, r₀, x);
        θ_min = real(asin((r₀-x)/MFP));
        θ_max = π/2;
        tmp0, = quadgk(fn_0,θ_min,θ_max);
        surf_fraction[k] = 2*real(tmp0)/(4*π*MFP^2);
    end
# println(surf_fraction)
    distance = zeros(1,length(r_int));
    distance_num = zeros(1,length(r_int));
    distance_den = zeros(1,length(r_int));

    for k in index[end]+1:length(r_int)
        x = r_int[k];

        θ_min = real(asin((r₀-x)/MFP));
        θ_max = π/2;
        α_min = 0;

        f1_0(θ) = f1(θ, MFP, r₀, x);
        g_0(θ) = g(θ, MFP, r₀, x);

        tmp0 = hcubature(f1_0,[θ_min; α_min], [θ_max; 1])
        distance_num[k] = real(tmp0[1])
        tmp0 = hcubature(g_0,[θ_min; α_min], [θ_max; 1])
        distance_den[k] = real(tmp0[1])
        distance[k] = distance_num[k]/distance_den[k]
    end

    τ = zeros(1,length(r_int))
    for k in index[end]+1:length(r_int)
        l = distance[k]
        fn_0(v) = fn_2(v, mass, kBT, distance[k])
        τ[k], = quadgk(fn_0, 0.0, Inf)
    end
    kwall = zeros(size(r_int))
    kwall[index[end]+1:end] = surf_fraction[index[end]+1:end]./τ[index[end]+1:end]/1e6
    return kwall
end
