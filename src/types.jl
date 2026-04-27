"""
        OBC
empty structure used to exploit julia automatic dispacth. It selects the function to use when Open Boundary Conditions applies
"""
struct OBC end
"""
        PBC
empty structure used to exploit julia automatic dispacth. It selects the function to use when Periodic Boundary Conditions applies
"""
struct PBC end


"""
    Check_flag

bitflag to select which `juobs.Corr` field *not* to check.

## Elements:
    - no_flag   = 0x00
    - no_gamma  = 0x01
    - no_obs    = 0x02
    - no_kappa  = 0x04
    - no_mu     = 0x08
    - no_y0     = 0x10
    - no_theta1 = 0x20
    - no_theta2 = 0x40
    - no_thetas = no_theta1 | no_theta2
"""
@bitflag Check_flag begin
    no_flag   = 0x00
    no_gamma  = 0x01
    no_obs    = 0x02
    no_kappa  = 0x04
    no_mu     = 0x08
    no_y0     = 0x10
    no_theta1 = 0x20
    no_theta2 = 0x40
end

no_thetas = no_theta1 | no_theta2
