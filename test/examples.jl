######################################
# End to end identification examples #
######################################
using Test: @testset, @test
using Apolo

@testset "Uniaxial Extension Example" begin

    # FIXME: Fix absolute versus relative paths
    include("../examples/uniaxial_extension/iden_uniaxial_extension.jl")

end
