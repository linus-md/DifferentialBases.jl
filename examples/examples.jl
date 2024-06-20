using DifferentialBases: differential_basis, akzo_nobel 

println("akzo_nobel")
ideal, derivatives, R = akzo_nobel()
diff_basis = differential_basis(ideal, derivatives, R, true)
println(diff_basis)
