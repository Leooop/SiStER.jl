# [dt_m]=SiStER_set_timestep(dx,dy,vx,vy,PARAMS)
# sets the advection time step
# J.-A. Olive, November 2014

function set_timestep(dx,dy,vx,vy,PARAMS)
    dt_m=PARAMS.fracCFL*0.5*min(minimum(dx),minimum(dy))./max(maximum(abs.(vx)),maximum(abs.(vy)))
end
