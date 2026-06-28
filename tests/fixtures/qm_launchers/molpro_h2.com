***,MATRIX Molpro launcher smoke H2
memory,64,m
geometry={
H 0.00000000 0.00000000 0.00000000
H 0.00000000 0.00000000 0.74000000
}
basis=sto-3g
hf
