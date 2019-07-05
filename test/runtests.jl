using Test

# Put tests here

function test_locate_markers_in_grid()
    xm = 10 .*rand(100)
    ym = 10 .*rand(100)
    x = [0.0, 1.0, 2.0 , 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
    y = [0.0, 1.0, 2.0 , 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
    dx = diff(x)
    dy = diff(y)
    qd,icn,jcn = locate_markers_in_grid(xm,ym,x,y,dx,dy)
    return (all(abs.(qd) .<= 4)) & (maximum(icn) == length(x)-1) & (maximum(jcn) == length(y)-1)
end

@testset "initialize file tests" begin
    @test test_locate_markers_in_grid()
end
