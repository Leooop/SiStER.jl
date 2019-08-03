# [T, rhs, Lii, Ljj, Lvv]=SiStER_thermal_solver_sparse_CFD(x,y,Told,rho,cp,k,dt,BCtherm,H)
# implicit Solver for thermal diffusion
# rho cp dT/dt = div (k grad T) + H
# J.-A. Olive, November 2014
# B.Z. Klein, added sparse matrix filling, End of 2014
# Added internal heat gen - howellsm 2/2017
# Updated to fully conservative (centered finite difference)
# for variable thermal diffusivity
# by S. Howell, now accepts k on shear nodes

function thermal_solver_sparse_CFD(x,y,Told,rho,cp,k,dt,BCtherm,H)

    Nx=length(x)
    Ny=length(y)
    dx=diff(x)
    dy=diff(y)

    Li=Int[]
    Lj=Int[]
    Lv=Float64[]
    rhs=zeros(Float64, Nx*Ny)

    n = 1

    for i=1:Ny
        for j=1:Nx

            gin=i+Ny*(j-1) # global index

            if i==1 # top boundary

                if BCtherm.top[1]==1 # Dirichlet

                    push!(Li, gin)
                    push!(Lj, gin)
                    push!(Lv, 1.0)

                    n = n+1

                    rhs[gin] = BCtherm.top[2]

                #L[gin,gin]=1
                #rhs[gin]=BCtherm.top[2]

                elseif BCtherm.top[1]==2 # Neumann

                    push!(Li, gin)
                    push!(Lj, gin)
                    push!(Lv, -1.0)

                    n = n+1

                    push!(Li, gin)
                    push!(Lj, gin+1)
                    push!(Lv, 1.0)

                    n = n+1

                    rhs[gin]=BCtherm.top[2]*dy[i]

    #             L(gin,gin]=-1
    #             L(gin,gin+1]=1
    #             rhs(gin]=BCtherm.top(2]*dy(i]

                end

            elseif i==Ny # bottom boundary


                if BCtherm.bot[1]==1 # Dirichlet

                    push!(Li, gin)
                    push!(Lj, gin)
                    push!(Lv, 1.0)

                    n = n+1

                    rhs[gin] = BCtherm.bot[2]

    #             L[gin,gin]=1
    #             rhs[gin]=BCtherm.bot[2]

                elseif BCtherm.bot[1]==2 # Neumann

                    push!(Li, gin)
                    push!(Lj, gin)
                    push!(Lv, 1.0)

                    n = n+1

                    push!(Li, gin)
                    push!(Lj, gin-1)
                    push!(Lv, -1.0)

                    n = n+1

                    rhs[gin] = BCtherm.bot[2]*dy[i-1]

    #             L[gin,gin]=1
    #             L[gin,gin-1]=-1
    #             rhs[gin]=BCtherm.bot(2]*dy(i-1]

                end

            elseif j==1 # left boundary


                if BCtherm.left[1]==2 # Neumann

                    push!(Li, gin)
                    push!(Lj, gin)
                    push!(Lv, -1.0)

                    n = n+1

                    push!(Li, gin)
                    push!(Lj, gin+Ny)
                    push!(Lv, 1.0)

                    n = n+1

                    rhs[gin] = BCtherm.left[2]*dx[j]

    #             L[gin,gin]=-1
    #             L[gin,gin+Ny]=1
    #             rhs[gin]=BCtherm.left[2]*dx[j]

                elseif BCtherm.left[1]==1 # Dirichlet

                    push!(Li, gin)
                    push!(Lj, gin)
                    push!(Lv, 1.0)

                    n = n+1

                    rhs[gin] = BCtherm.left[2]

    #             L[gin,gin]=1
    #             rhs[gin]=BCtherm.left[2]

                end


            elseif j==Nx # right boundary


                if BCtherm.right[1]==2 # Neumann

                    push!(Li, gin)
                    push!(Lj, gin)
                    push!(Lv, 1.0)

                    n = n+1

                    push!(Li, gin)
                    push!(Lj, gin-Ny)
                    push!(Lv, -1.0)

                    n = n+1

                    rhs[gin] = BCtherm.right[2]*dx[j-1]

    #             L[gin,gin]=1
    #             L[gin,gin-Ny]=-1
    #             rhs[gin]=BCtherm.right[2]*dx[j-1]

                elseif BCtherm.right[1]==1 # Dirichlet

                    push!(Li, gin)
                    push!(Lj, gin)
                    push!(Lv, 1.0)

                    n = n+1

                    rhs[gin] = BCtherm.right[2]

    #             L[gin,gin]=1
    #             rhs[gin]=BCtherm.right[2]

                end



            else


    #         #   internal nodes

                ddx=dx[j-1]+dx[j]
                ddy=dy[i-1]+dy[i]

                # Follows Gerya lettering scheme, P139-140
                kA    = (2 * k[i,j-1] * k[i,j]   )/( k[i,j-1] + k[i,j]   )
                kB    = (2 * k[i,j]   * k[i,j+1] )/( k[i,j]   + k[i,j+1] )
                kC    = (2 * k[i-1,j] * k[i,j]   )/( k[i-1,j] + k[i,j]   )
                kD    = (2 * k[i,j]   * k[i+1,j] )/( k[i,j]   + k[i+1,j] )

                push!(Li, gin)
                push!(Lj, gin)
                push!(Lv, rho[i,j]*cp[i,j]+2*dt*kB/(dx[j]*ddx) + 2*dt*kA/(dx[j-1]*ddx) + 2*dt*kD/(ddy*dy[i]) + 2*dt*kC/(ddy*dy[i-1]))


                n = n+1

    #             Lv(n) = rho(i,j)*cp(i,j)+2*dt*kB/(dx(j)*ddx) + 2*dt*kA/(dx(j-1)*ddx) + ...
    #                 2*dt*kD/(ddy*dy(i)) + 2*dt*kC/(ddy*dy(i-1))

                push!(Li, gin)
                push!(Lj, gin+Ny)
                push!(Lv, -2*dt*kB/(dx[j]*ddx))

                n = n+1

    #             L[gin,gin+Ny)=-2*dt*k[i,j)/[dx[j)*ddx)
    #
                push!(Li, gin)
                push!(Lj, gin-Ny)
                push!(Lv, -2*dt*kA/(dx[j-1]*ddx))


                n = n+1
    #
    # #             L[gin,gin-Ny]=-2*dt*k[i,j-1]/[dx[j-1]*ddx]
    #
                push!(Li, gin)
                push!(Lj, gin+1)
                push!(Lv, -2*dt*kD/(dy[i]*ddy))


                n = n+1
    #
    # #             L[gin,gin+1]=-2*dt*k[i,j]/[dy[i]*ddy]
    #
                push!(Li, gin)
                push!(Lj, gin-1)
                push!(Lv, -2*dt*kC/(dy[i-1]*ddy))

                n = n+1

    #             L[gin,gin-1]=-2*dt*k[i-1,j]/[dy[i-1]*ddy]

                rhs[gin]=rho[i,j]*cp[i,j]*Told[i,j] + H[i,j]*dt

            end

        end
    end

    nn = n-1

    Lii = Li[1:nn]
    Ljj = Lj[1:nn]
    Lvv = Lv[1:nn]

    L = sparse(Lii, Ljj, Lvv)

    tvec=L\rhs

    T=reshape(tvec,Ny,Nx)

    #return T, rhs, Lii, Ljj, Lvv
    return T, L, rhs

end
