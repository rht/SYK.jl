using Statistics
include("./majorana.jl")
include("./kitaev_H.jl")
#using majorana
#using kitaev_H

function remove_matrix_elements_randomly(H, dim, num_remove)
    for x = 1:num_remove
        i, j = rand(1:dim, 2)
        H[i, j] = 0
    end
end

function get_spectral_form_factor(energies, beta)
    # Z is partition function
    Z = sum([exp(-beta * E) for E in energies])
    return real(conj(Z) * Z)
end

# This code generates, diagonalizes and stores SYK eigenvalues
# size parameters
N = 12  # set number of Majorana fermions
numH = 5  # set number of samples
Htype = "SYK"  # type of the Hamiltonian used
matrix_size = Int(2^(N/2))

# SYK parameter
J = 1

# get Majorana fermions
Xi = majorana(N)


ts = exp10.(range(-2,stop=2.5,length=50))
gss = zeros(numH, length(ts))
for i = 1:numH
    # diagonalize H once
    H = kitaev_H(J, Xi, N) # creates the Hamiltonian
    remove_matrix_elements_randomly(H, matrix_size, 20)  # make it sparse
    VV, DD = eigen(H)  # diagonalizes the Hamiltonian

    gs = []
    for (j, t) in enumerate(ts)
        # begin calculating spectral form factor g = ZxZ*
        # for the definition of spectral form factor, see
        # https://physics.stackexchange.com/questions/481786/what-is-the-spectral-form-factor
        g = get_spectral_form_factor(VV, 1im * t)
        gss[i, j] = g
        # end spectral form factor
    end
end

# average across hamiltonians
final_gs = mean(gss, dims=1)
println(ts)
println(final_gs[:])
