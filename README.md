# GTOC11Utils
## Installation
```julia
(@v1.6)> add https://github.com/buaa1205/GTOC11Utils
```

## Use
Note: λ can be typed in the Julia REPL in VSCode as `\lambda<TAB>` and ṙ can be typed as `r\dot<TAB>`.
```julia
using GTOC11Utils


# Define state vectors
chaser = [-1.078389739950731, -2.4354584930816596, -0.0, 3.5201927832810287, -1.5586961514320443, -0.0]

target = [-1.0797379621096455, -2.5442185136451623, -0.14643380352709984, 3.1756439174678857, -1.3357628170506850, -0.27747062891731783]


# Solve the optimization problem 
out = low_thrust_transfer(chaser, target)


# Inspect returned optimal values in the `u` property
out.u.t     # 1.6238469777682614
out.u.λ.r.x   # -2.4348070681482498e11
out.u.λ.ṙ.y   # -2.2321843456764603e11


# Run the problem with the optimal values to get the full solution
sol = run_solution(out)


# Plot the solution
chase_plot(sol; show_accel=true)    # Plot the trajectories in 3D
acceleration_plot(sol)              # Plot the accelerations vs. time




```
