using DifferentialBases: differential_basis, double_pendulum

println("Double Pendulum")
ideal, derivatives, R = double_pendulum()
diff_basis = differential_basis(ideal, derivatives, R, true)
println(diff_basis)

