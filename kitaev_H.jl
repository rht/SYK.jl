#module kitaev_H
#export kitaev_H

function kitaev_H(J, Xi, N)
    #KITAEV_H returns Hamiltonian for Kitaev's chaotic Hamiltonain.
    #   Specifically, this builds a nonlocal system of N majorana fermions with
    #   random interactions chosen from a Gaussian with variance J^2.
    #   
    #   Let Xi_i be a Majorana, where i runs over 1..N, N = 2*log2(length(Xi)),
    #   and {Xi_i, Xi_j} = delta_ij. The Hamiltonian is given by
    #
    #   H =  \Sum_{i<j<k<l} J_{ijkl} Xi_i Xi_j Xi_k Xi_l
    #
    #   N.B. HL comes with a negative sign.
    #   N.B.2. The Hilbert space is 2^(N/2) dimensional
    #
    #   Outputs:
    #   H is the Hamiltonian
    #
    #   Inputs:
    #   N is the number of majorana fermions. (Twice the number of qubits)
    #   J is the spread of the coupling distribution
    #   Xi are the fermions

    if N % 2 != 0
        error("N must be even")
    end

    if N < 4
        error("N must be greater than 4")
    end

    # build Hamiltonian
    H_length = Int(2^(N/2))
    H = zeros(H_length, H_length)

    var_J_jklm = factorial(3)*J^2/( (N-1)*(N-2)*(N-3) )
    sigma = sqrt(var_J_jklm)

# what he writes

#     for j = 1:N
#         
#         for k = 1:N
#             
#             if k == j
#                 continue
#             end
#             
#             for l = 1:N
#                 
#                 if l == j || l == k
#                     continue
#                 end
#                 
#                 for m = 1:N
#                     
#                     if m == j || m == k || m == l
#                         continue
#                     end
#                     
#                     J_jklm = normrnd(0, sigma);
#                     H = H + 1/factorial(4) * J_jklm * ...
#                                    Xi(:,:,j)*Xi(:,:,k)*Xi(:,:,l)*Xi(:,:,m);
#                 end
#             end
#         end
#     end

    # this is faster

    for j = 1:N
        for k = (j+1):N
            for l = (k+1):N
                for m = (l+1):N
                    J_jklm = randn() * sigma
                    H = H + J_jklm * Xi[:,:,j]*Xi[:,:,k]*Xi[:,:,l]*Xi[:,:,m]
                end
            end
        end
    end
    return H
end
#end
