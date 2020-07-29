# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 17:46:56 2020

@author: luzhiwei
"""

# 画能带

# from band import draw_band_structure
# pos_file='POSCAR'
# eigenval_filename='EIGENVAL'
# dos_file='DOSCAR'
# kpoint_filename='KPOINTS'
# yaxis_max=5
# yaxis_min=-3  #调整y轴范围
# draw_band_structure(eigenval_filename, dos_file, pos_file, kpoint_filename,yaxis_max,yaxis_min)


# # 画原子pdos
from band import draw_pdos
pos_file='POSCAR'
eigenval_filename='EIGENVAL'
dos_file='DOSCAR'
kpoint_filename='KPOINTS'
yaxis_max=5
yaxis_min=-3  #调整y轴范围
draw_pdos(eigenval_filename, dos_file, pos_file, kpoint_filename,yaxis_max,yaxis_min)

# 画band-pdos
# from band import draw_band_pdos
# pos_file='POSCAR'
# eigenval_filename='EIGENVAL'
# dos_file='DOSCAR'
# kpoint_filename='KPOINTS'
# yaxis_max=5
# yaxis_min=-3  #调整y轴范围
# draw_band_pdos(eigenval_filename, dos_file, pos_file, kpoint_filename,yaxis_max,yaxis_min)
