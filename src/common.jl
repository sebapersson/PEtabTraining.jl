
function _makechunks(x::AbstractVector, n::Integer)
    c = length(x) ÷ n
    return [x[(1 + c * k):(k == n - 1 ? end : c * k + c)] for k in 0:(n - 1)]
end
