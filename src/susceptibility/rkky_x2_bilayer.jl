const X1_INDEX = 1
const X2_INDEX = 3

function x2_layer_channels(chi::AbstractMatrix)
    size(chi) == (4, 4) || throw(DimensionMismatch("chi must be 4x4"))
    chi11 = chi[X1_INDEX, X1_INDEX]
    chi22 = chi[X2_INDEX, X2_INDEX]
    chi12 = chi[X1_INDEX, X2_INDEX]
    chi21 = chi[X2_INDEX, X1_INDEX]
    return (;
        intra = (chi11 + chi22) / 2,
        inter = (chi12 + chi21) / 2,
        intra_layer_residual = abs(chi11 - chi22),
        inter_layer_residual = abs(chi12 - chi21),
    )
end

function x2_layer_channels_q_index(
    qx_index::Integer,
    qy_index::Integer,
    cache,
    mu;
    temperature = 0.001,
    tol = 1e-8,
)
    Nk = cache.Nk
    chi11 = 0.0 + 0.0im
    chi22 = 0.0 + 0.0im
    chi12 = 0.0 + 0.0im
    chi21 = 0.0 + 0.0im

    for ix in 1:Nk, iy in 1:Nk
        ixq = mod1(ix + qx_index, Nk)
        iyq = mod1(iy + qy_index, Nk)
        for n in 1:4, m in 1:4
            kernel = susceptibility_kernel(
                cache.energies[ix, iy, n],
                cache.energies[ixq, iyq, m],
                mu,
                temperature;
                tol,
            )
            u1 = cache.vectors[ix, iy, X1_INDEX, n]
            u2 = cache.vectors[ix, iy, X2_INDEX, n]
            v1 = cache.vectors[ixq, iyq, X1_INDEX, m]
            v2 = cache.vectors[ixq, iyq, X2_INDEX, m]
            chi11 += kernel * abs2(u1) * abs2(v1)
            chi22 += kernel * abs2(u2) * abs2(v2)
            chi12 += kernel * conj(u1) * v1 * conj(v2) * u2
            chi21 += kernel * conj(u2) * v2 * conj(v1) * u1
        end
    end

    scale = inv(Nk^2)
    chi11 *= scale
    chi22 *= scale
    chi12 *= scale
    chi21 *= scale
    return (;
        intra = (chi11 + chi22) / 2,
        inter = (chi12 + chi21) / 2,
        intra_layer_residual = abs(chi11 - chi22),
        inter_layer_residual = abs(chi12 - chi21),
    )
end

function compute_x2_bilayer_rkky(
    U;
    Nk::Integer = 100,
    Nq::Integer = Nk,
    n_target = 1.5,
    temperature = 0.01,
    JK = 1.0,
    tol = 1e-8,
)
    Nk > 0 || throw(ArgumentError("Nk must be positive"))
    Nq == Nk || throw(ArgumentError("Nq must equal Nk for indexed cached momenta"))
    temperature > 0 || throw(ArgumentError("temperature must be positive"))
    isfinite(JK) || throw(ArgumentError("JK must be finite"))

    mu = find_chemical_potential(n_target, U; Nk, temperature)
    filling = filling_number(mu, U; Nk, temperature)
    cache = diagonalized_kmesh(Nk, U)
    contributions = zeros(ComplexF64, Nq^2, 5)
    intra_residuals = zeros(Float64, Nq^2)
    inter_residuals = zeros(Float64, Nq^2)

    Threads.@threads for linear_index in 1:(Nq^2)
        iqx = div(linear_index - 1, Nq)
        iqy = mod(linear_index - 1, Nq)
        qx = 2pi * iqx / Nq
        qy = 2pi * iqy / Nq
        channels = x2_layer_channels_q_index(
            iqx,
            iqy,
            cache,
            mu;
            temperature,
            tol,
        )
        contributions[linear_index, 1] = channels.inter
        contributions[linear_index, 2] = exp(im * qx) * channels.intra
        contributions[linear_index, 3] = exp(im * (qx + qy)) * channels.intra
        contributions[linear_index, 4] = exp(im * 2qx) * channels.intra
        contributions[linear_index, 5] = exp(im * qx) * channels.inter
        intra_residuals[linear_index] = channels.intra_layer_residual
        inter_residuals[linear_index] = channels.inter_layer_residual
    end

    fourier = vec(sum(contributions; dims = 1)) ./ Nq^2
    bare = -real.(fourier)
    physical = JK^2 .* bare
    return (;
        U = Float64(U),
        mu,
        filling,
        filling_error = filling - n_target,
        bare_J0 = bare[1],
        bare_J1 = bare[2],
        bare_J2 = bare[3],
        bare_J3 = bare[4],
        bare_J1p = bare[5],
        J0 = physical[1],
        J1 = physical[2],
        J2 = physical[3],
        J3 = physical[4],
        J1p = physical[5],
        max_intra_layer_residual = maximum(intra_residuals),
        max_inter_layer_residual = maximum(inter_residuals),
        max_fourier_imaginary_residual = maximum(abs, imag.(fourier)),
    )
end

function x2_rkky_metadata(
    U_values;
    Nk::Integer,
    Nq::Integer,
    n_target,
    temperature,
    JK,
)
    return (;
        schema_version = 1,
        basis = "x1,z1,x2,z2",
        itinerant_orbital = "x2-y2",
        local_moment_orbital = "z2",
        fourier_convention = "chi(R)=sum_q exp(+i q.R) chi(q)/Nq^2",
        susceptibility_unit = "eV^-1",
        exchange_unit = "eV",
        sign_convention = "H=sum_ij J_ij S_i.S_j; positive J is antiferromagnetic",
        channel_intra = "(chi[1,1]+chi[3,3])/2",
        channel_inter = "(chi[1,3]+chi[3,1])/2",
        U_values = Float64.(collect(U_values)),
        Nk = Int(Nk),
        Nq = Int(Nq),
        n_target = Float64(n_target),
        temperature = Float64(temperature),
        JK = Float64(JK),
    )
end

function validate_x2_rkky_resume(expected, actual)
    expected == actual || throw(ArgumentError("RKKY resume metadata mismatch"))
    return nothing
end

function save_x2_rkky_result_atomic(filepath::AbstractString, result)
    mkpath(dirname(filepath))
    temporary = filepath * ".tmp"
    isfile(temporary) && rm(temporary; force = true)
    try
        jldsave(temporary; result)
        mv(temporary, filepath; force = true)
    finally
        isfile(temporary) && rm(temporary; force = true)
    end
    return filepath
end

load_x2_rkky_result(filepath::AbstractString) = load(filepath, "result")

function _empty_x2_rkky_result(metadata)
    return (;
        metadata,
        completed_U_values = Float64[],
        mu = Float64[],
        filling = Float64[],
        filling_error = Float64[],
        bare_J0 = Float64[],
        bare_J1 = Float64[],
        bare_J2 = Float64[],
        bare_J3 = Float64[],
        bare_J1p = Float64[],
        J0 = Float64[],
        J1 = Float64[],
        J2 = Float64[],
        J3 = Float64[],
        J1p = Float64[],
        max_intra_layer_residual = Float64[],
        max_inter_layer_residual = Float64[],
        max_fourier_imaginary_residual = Float64[],
    )
end

function _validate_x2_result(result, requested_U_values)
    fields = propertynames(result)
    lengths = [
        length(getproperty(result, field)) for field in fields if
        field != :metadata && getproperty(result, field) isa AbstractVector
    ]
    isempty(lengths) || all(==(first(lengths)), lengths) ||
        throw(ArgumentError("RKKY result arrays have inconsistent lengths"))
    n_done = length(result.completed_U_values)
    n_done <= length(requested_U_values) ||
        throw(ArgumentError("RKKY result is longer than the requested scan"))
    result.completed_U_values == requested_U_values[1:n_done] ||
        throw(ArgumentError("completed U values are not a requested prefix"))
    for field in fields
        value = getproperty(result, field)
        if field != :metadata && value isa AbstractVector{<:Number}
            all(isfinite, value) || throw(ArgumentError("non-finite RKKY result in $field"))
        end
    end
    return nothing
end

function scan_x2_bilayer_rkky(
    U_values;
    Nk::Integer = 100,
    Nq::Integer = Nk,
    n_target = 1.5,
    temperature = 0.01,
    JK = 1.0,
    filepath::AbstractString,
    point_calculator = compute_x2_bilayer_rkky,
)
    requested = Float64.(collect(U_values))
    metadata = x2_rkky_metadata(
        requested;
        Nk,
        Nq,
        n_target,
        temperature,
        JK,
    )
    result = if isfile(filepath)
        loaded = load_x2_rkky_result(filepath)
        validate_x2_rkky_resume(metadata, loaded.metadata)
        loaded
    else
        _empty_x2_rkky_result(metadata)
    end
    _validate_x2_result(result, requested)

    for U in requested[(length(result.completed_U_values) + 1):end]
        point = point_calculator(U; Nk, Nq, n_target, temperature, JK)
        additions = Dict{Symbol, Any}(:completed_U_values => U)
        for field in propertynames(result)
            if field != :metadata && field != :completed_U_values
                additions[field] = getproperty(point, field)
            end
        end
        result = (; (
            field => field == :metadata ? result.metadata :
                [getproperty(result, field); additions[field]]
            for field in propertynames(result)
        )...)
        _validate_x2_result(result, requested)
        save_x2_rkky_result_atomic(filepath, result)
    end
    return result
end
