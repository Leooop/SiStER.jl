# [p, vx, vy]=SiStER_reshape_solver_output(S,Kc,Nx,Ny)
# reshapes the output of the Stokes solver into pressure and velocity
# arrays
# B.Z. Klein, rewritten from J.-A. Olive, Spring 2013

function reshape_solver_output(S,Kc,Nx,Ny)

    INDp  = 1:3:length(S)
    INDvx = INDp .+ 1
    INDvy = INDp .+ 2

    p  = reshape(S[INDp], Ny, Nx).*Kc
    vx = reshape(S[INDvx], Ny, Nx)
    vy = reshape(S[INDvy], Ny, Nx)

    return p, vx, vy

end
