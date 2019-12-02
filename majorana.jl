include("./tensor.jl")
using LinearAlgebra
#module majorana
# using tensor
#export majorana

function majorana(N::Int)
    #MAJORANA creates a representation of N majorana fermions
    #   Creates N Majorana fermions in terms of spin chain variables (Pauli
    #   sigma matrices). The only catch is that they satisfy the Dirac algebra
    #   {Xi_i, Xi_j} = delta_{ij}. This is basically Jordan-Wignering.
    #
    #   N.B. the Hilbert space is 2^(N/2) dimensional
    #
    #   Outputs:
    #   Xi is an 2^(N/2) x 2^(N/2) x N array. 
    #   You access the ith fermion as Xi(:,:,i)
    #
    #   Inputs:
    #   N is the number of majorana fermions. (Twice the number of qubits)
    #       thus, N should be even

    if N % 2 != 0
        error("N must be even")
    end

    #SU(2) Matrices
    X = [0 1; 1 0]

    Y = [0 -1im; 1im 0]

    Z = [1 0; 0 -1]

    # output of majoranas here
    dim = Int(2^(N/2))
    Xi = zeros(Complex, (dim, dim, N))

    # growing chain of Xs, start with 0
    Xs = Matrix{Complex}(I, dim, dim)

    for i in 1:N
        # the Y or Z at the end
        if i % 2 == 1
            xi = tensor(Z, (i+1)/2 , N/2)
        else
            xi = tensor(Y, i/2 , N/2)
        end

        Xi[:,:,i] = Xs*xi

        # build an increment chains of X's
        if i % 2 == 0
            Xs = Xs*tensor(X, i/2, N/2)
        end
    end
    # right now, it's normalized {Xi_i, Xi_j} = 2*delta_{ij}
    # but Kitaev's normalization doens't have the 2
    Xi = Xi / sqrt(2)


    # this code checks the commutation relation {Xi_i, Xi_j} = delta_{ij}

    #zeromat = zeros(2^(N/2));
    #Imat = eye(2^(N/2));

    #for i = 1:N
    #    for j = 1:N
    #        if i != j
    #            @assert min(min(acom(Xi(:,:,i), Xi(:,:,j)) == zeromat))
    #        else
    #            @assert min(min(abs(acom(Xi(:,:,i), Xi(:,:,j)) - Imat))) < 1e-10
    #        end
    #    end
    #end
    return Xi
end
#end
