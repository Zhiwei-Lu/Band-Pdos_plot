from band import draw_band_structure
pos_file='POSCAR'
eigenval_filename='EIGENVAL'
dos_file='DOSCAR'
kpoint_filename='KPOINTS2'
yaxis_max=2
yaxis_min=-2  #调整y轴范围
draw_band_structure(eigenval_filename, dos_file, pos_file, kpoint_filename,yaxis_max,yaxis_min)
