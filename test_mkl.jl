function eigvals(x)
    if ndims(x) > 1  
        x = reshape(x, :)  
    end
    dim = Int32(sqrt(size(x, 1)))
    c_eigvals = ccall((:eigvals, "./lib/libjulia.so"), Ptr{Cdouble}, (Ptr{Cdouble}, Cint), x, dim)
    eigenvalues = unsafe_wrap(Array, c_eigvals, (dim,), own=true)
    return copy(eigenvalues)
end

function eigvecs(x)
    if ndims(x) > 1  
        x = reshape(x, :)  
    end
    dim = Int32(sqrt(size(x, 1)))
    c_eigvecs = ccall((:eigvecs, "./lib/libjulia.so"), Ptr{Cdouble}, (Ptr{Cdouble}, Cint), x, dim)
    eigenvectors = unsafe_wrap(Array, c_eigvecs, (dim, dim), own=true)
    return copy(eigenvectors)
end

function matmul(A, B)
    m, k = size(A)
    _, n = size(B)
    if size(A, 2) != size(B, 1)
        error("The dimensions of A and B don't match for multiplication!")
    end
    A = reshape(A, :)
    B = reshape(B, :)
    c_result = ccall((:matmul, "./lib/libjulia.so"), Ptr{Cdouble}, (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cint, Cint), A, B, m, n, k)
    C = unsafe_wrap(Array, c_result, (m, n), own=true)  
    return copy(C)
end

function cholesky(A)
    if ndims(A) > 1  
        A = reshape(A, :)  
    end
    dim = Int32(sqrt(size(A, 1)))
    try
        c_result = ccall((:cholesky, "./lib/libjulia.so"), Ptr{Cdouble}, (Ptr{Cdouble}, Cint), A, dim)
        B = unsafe_wrap(Array, c_result, (dim, dim), own=true)
        return copy(B)
    catch e
        println("Error: ", e)
        return nothing
    end
end


function inverse(A)
    if ndims(A) > 1  
        A = reshape(A, :)  
    end
    dim = Int32(sqrt(size(A, 1)))
    c_result = ccall((:inverse, "./lib/libjulia.so"), Ptr{Cdouble}, (Ptr{Cdouble}, Cint), A, dim)
    B = unsafe_wrap(Array, c_result, (dim, dim), own=true)
    return copy(B)
end

x = Float64.([1 -3 3; 3 -5 3; 6 -6 4])
I = Float64.([1 0 0; 0 1 0; 0 0 1])

eigvalues = eigvals(x)
eigvectors = eigvecs(x)

for i in 1:3
    println("lambda ", i, " is ", eigvalues[i])
    println("corresponding eigenvector: ", eigvectors[:,i])
end

C = matmul(x, I)
println(C)

x2 = Float64.([4 12 -16; 12 37 -43; -16 -43 98])
D = cholesky(x2)
println(D)


E = inverse(I)
println(E)
