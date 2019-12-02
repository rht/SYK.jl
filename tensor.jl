#module tensor
#export tensor
using LinearAlgebra

function tensor(a, i, N, p=2)
    #TENSOR takes matrix a and implements it as the ith qubit
    #   Embeds the matrix a in a p^N dimensional Hilbert space of N qudits
    #   in the ith position in the tensor product.
    #   It is assumed a is a pxp 1-qudit matrix.
    #
    #   a is a pxp Unitary matrix to operate on a single qudit
    #   i is the position in the tensor product, indexing from 1
    #   N is the number of qudits in the Hilbert space
    #   varargin is an optional arg to specify the degree of freedom per site. 
    #   by default it's qubits, i.e. varargin{1} = 2.
    dim1 = Int(p^(i-1))
    dim2 = Int(p^(N-i))

    A = kron(kron(Matrix{Complex}(I, dim1, dim1), a), Matrix{Complex}(I, dim2, dim2));
    return A
end
#end
