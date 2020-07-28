using LambertW
# TODO: Add documentation
function lamb_smooth_minimum(
    l::AbstractArray;
    lower_bound::FT = 0.1,
    frac_upper_bound::FT = 1.5,
) where {FT}

    leng = size(l)
    xmin = minimum(l)

    lambda0 = max(
        xmin * lower_bound /
        real(LambertW.lambertw(FT(2) / MathConstants.e)),
        frac_upper_bound,
    )

    i = 1
    num = 0
    den = 0
    while (tuple(i) < leng)
        num += l[i] * exp(-(l[i] - xmin) / lambda0)
        den += exp(-(l[i] - xmin) / lambda0)
        i += 1
    end
    smin = num / den

    return smin
end
