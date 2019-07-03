# seed_markers_uniformly(x,y,dx,dy,N)
# seeds N markers in the cell or quadrant whose upper-left node
# has coordinates (x,y) and width, height (dx,dy)

function seed_markers_uniformly(x,y,dx,dy,N)

    # create a subgrid

    fact=0.4 # randomization factor

    nx=Int(ceil(sqrt(N*dx/dy)))
    ny=Int(ceil(N/nx))
    actual_N=nx*ny

    ddx=range(x,x+dx,length = nx+1)
    ddy=range(y,y+dy,length = ny+1)

    #[DX, DY]=meshgrid(ddx,ddy); matlab version
    DX = repeat(ddx',outer = size(ddy))
    DY = repeat(ddy,outer = size(ddx'))

    xm=zeros(Float64,actual_N)
    ym=copy(xm)

    k=0
    for i=1:ny
        for j=1:nx

            k=k+1
            xsub=DX[i,j]
            ysub=DY[i,j]
            dxsub=dx/nx
            dysub=dy/ny

            # randomize marker position around sub-cell center
            xm[k]=xsub+(dxsub/2)*(1+fact*2*(rand()-0.5))
            ym[k]=ysub+(dysub/2)*(1+fact*2*(rand()-0.5))

        end
    end


    if actual_N>N # only keep N random markers out of actual_N markers

        idx=randperm(actual_N)
        idx=idx[1:N]
        xm=xm[idx]
        ym=ym[idx]

    end

    return xm, ym

end
