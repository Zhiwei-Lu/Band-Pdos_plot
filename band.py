import matplotlib.pyplot as plt
import numpy as np
import re
from matplotlib.pyplot import MultipleLocator
from math import sqrt, pow, sin, cos, pi
import sys
import math 
def draw_band_structure(eigenval_filename, dos_file, pos_file, kpoint_filename,yaxis_max,yaxis_min):
    yaxis_max=yaxis_max
    yaxis_min=yaxis_min
    engval,kpoints,ispin_type=read_eigenval(eigenval_filename)
    E_fermi = get_fermi_from_doscar(dos_file)
    engval=engval-E_fermi
    latt_rec=read_recip(pos_file)
    latt,latt_vec,atom,S,sys_name=read_poscar(pos_file)
    hsp,hsp_label,node=read_high_sym_point(kpoint_filename)
    temp=[]
    for ii in range(1,kpoints.shape[0]):
        temp.append(np.linalg.norm((kpoints[ii-1]-kpoints[ii])*latt_rec))
    temp.insert(0,0)
    temp=np.cumsum(temp)
    temp=temp.reshape(kpoints.shape[0],1)
    kpoints=np.append(kpoints,temp,axis=1)
    aa=[]
    for ii in range(1,int(hsp.shape[0]/2)+1):
        aa.append(kpoints[node*(ii-1),:])
        aa.append(kpoints[node*ii-1,:])    
    new_array = np.array(aa)
    x_scale=[]
    x_scale=np.array(list(set([tuple(x) for x in new_array])))
    coordskeys2=[]
    coordskeys3=[]
    new_array2=list(new_array)
    ii=0
    for aa in new_array:
        bb=np.sum(aa)
        ii=ii+1
        if bb not in coordskeys2:            
            coordskeys2.append(bb)
            coordskeys3.append(ii)
        coordskeys5=[x-1 for x in coordskeys3]   #编号
    coordskeys5=np.array(coordskeys5)
    coordskeys6=[]
    for ii in range(int(len(coordskeys5))):
        coordskeys6.append(tuple(hsp_label[coordskeys5[ii]])) 
    coef = np.array(coordskeys6).flatten()
    kk=[]
    for bb in coef:
        kk.append( re.findall(r'\b\w', bb))
    hsp_label2=np.array(kk).flatten()          #k点符号
    if ispin_type==1:
        plt.figure(dpi=300)
        energy=np.zeros([engval.shape[0],int(engval.shape[1]/2)])
        for ii in range(energy.shape[0]): 
            for jj in range(energy.shape[1]):
                energy[ii][jj]=engval[ii][int(2*jj)]
    #    print(coordskeys2,coordskeys5)        
        plot_band(kpoints,x_scale,hsp_label2,E_fermi,sys_name,energy,'r',yaxis_max,yaxis_min)
    
    else:
        plt.figure(dpi=300)
        plt.subplot(121)
        sys_name='spin_up'
        energy_up=np.zeros([engval.shape[0],int(engval.shape[1]/2)])
        for ii in range(energy_up.shape[0]): 
            for jj in range(energy_up.shape[1]):
                energy_up[ii][jj]=engval[ii][int(2*jj)][0]            
        plot_band(kpoints,x_scale,hsp_label2,E_fermi,sys_name,energy_up,'b',yaxis_max,yaxis_min)
        plt.subplot(122)
        sys_name='spin_down'
        energy_down=np.zeros([engval.shape[0],int(engval.shape[1]/2)])
        for ii in range(energy_down.shape[0]): 
            for jj in range(energy_down.shape[1]):
                energy_down[ii][jj]=engval[ii][int(2*jj)][1]
        plot_band(kpoints,x_scale,hsp_label2,E_fermi,sys_name,energy_down,'r',yaxis_max,yaxis_min)
    
    plt.savefig('Band.jpg') 
def get_fermi_from_doscar(dos_file):
    f=open(dos_file)
    for ii in range(6):
        line=f.readline()
    s=np.array(line.split())
    E_fermi=float(s[3])
    return E_fermi    
    
def plot_band(kpoints,x_scale,hsp_label2,E_fermi,sys_name,energy,color,yaxis_max,yaxis_min):          
    plt.plot(kpoints[:,3],energy,color=color)
    for ii in range(int(x_scale.shape[0])):
        plt.axvline(x=x_scale[ii,3], color='r', linestyle='--')          
    coordskeys7=np.array(sorted(list(x_scale[:,3])))
    plt.xticks(coordskeys7,hsp_label2)  
    plt.axhline(0, color='r', linestyle='--') 
    x_max=coordskeys7[int(coordskeys7.shape[0]-1)]
    x_min=coordskeys7[0]
    y_max=yaxis_max
    y_min=yaxis_min
    y_major_locator=MultipleLocator(1)
    plt.gca().yaxis.set_major_locator(y_major_locator)
    plt.xlim(x_min,x_max)
    plt.ylim(y_min,y_max)
    plt.title(sys_name,y=1)

def read_doscar(dos_file):    
    f=open(dos_file)
    line=f.readline()
    s1=np.array(line.split())
    for ii in range(5):
        line=f.readline()
    s=np.array(line.split())
    NEDOS=int(s[2])
    n_element=int(s1[0])
    line=f.readline()
    line=f.readline()
    line=f.readline()
    split=line.split()
    if len(split)==3:
        ispin_type=1
    else:
        ispin_type=2
    sum_dos=np.zeros([NEDOS,len(split)])
    f.close()
    f=open(dos_file)
    for ii in range(6):
        line=f.readline()
    for ii in range(NEDOS):
        line=f.readline()
        split=line.split()
        sum_dos[ii][:]=[float(split[jj]) for jj in range(len(split))]
    line=f.readline()
    line=f.readline()
    split=line.split()
    if not line:
        have_pdos=False
    else:
        have_pdos=True
    if have_pdos:
        n_pdos = len(split);
    else:
        n_pdos = 0;
    f.close()
    if n_pdos!=0:
        f=open(dos_file)
        for nn in range(NEDOS+6):
            line=f.readline()
        p_dos=np.zeros([n_element,NEDOS,n_pdos])
        temp_dos=[]  
        for ii in range(n_element):
            line=f.readline()
            temp_dos=[] 
            for jj in range(NEDOS):
                line=f.readline()
                split=line.split()
                p_dos[ii][jj][:]=[float(split[kk]) for kk in range(len(split))]
    else:
        p_dos=0
    
    return sum_dos,p_dos    
   
def read_eigenval(eigenval_filename):
    f=open(eigenval_filename)
    line=f.readline()
    split=line.split()
    ispin_type=int(split[-1])
    for ii in range(5):
        line=f.readline()
    split=line.split()
    nbands=int(split[-1])
    nkpoints=int(split[1])    
    kpoints=[]
    if ispin_type==1:
        engval=np.zeros([nkpoints,nbands*2])
#        print(engval.size)
        for ii in range(nkpoints):
            line=f.readline()
            line=f.readline()
            split=line.split()
            kpoints.append([float(split[j]) for j in range(3)])
            for jj in range(nbands):
                line = f.readline()
                split=line.split()
                engval[ii][2*jj]=[float(split[kk]) for kk in range(3)][1]
                engval[ii][2*jj+1]=[float(split[kk]) for kk in range(3)][-1]
        kpoints=np.array(kpoints)
#        engval=engval-E_fermi
    else:        
        engval= np.zeros([nkpoints,nbands*2,2])
        for ii in range(nkpoints):
            line=f.readline()
            line=f.readline()
            split=line.split()
            kpoints.append([float(split[j]) for j in range(3)])
            for jj in range(nbands):
                line = f.readline()
                split=line.split()
                engval[ii][2*jj][0]=float(split[1])
                engval[ii][2*jj+1][0]=float(split[3])
                engval[ii][2*jj][1]=float(split[2])
                engval[ii][2*jj+1][1]=float(split[4])                    
        kpoints=np.array(kpoints)    
#        engval=engval-E_fermi
    f.close()
    return engval,kpoints,ispin_type    
    
def read_high_sym_point(kpoint_filename):
    f=open(kpoint_filename) 
    hsp_1=[]
    hsp_label=[]
    k=1;
    kk=1
    line1=f.readline()
    line2=f.readline()
    mat=re.findall('\d+',line2)
    node=int(mat[0])
#    print(node)
    mat2=[]
    for line in f.readlines()[2:]:
        aa=line[0:line.rfind('!')]
        hsp_1.append(aa.split())       
        bb=line[line.rfind('!')+1:]
        hsp_label.append(bb.split())
    while [] in hsp_1 and hsp_label:
        hsp_1.remove([])
        hsp_label.remove([])
    hsp=[]
    for ii in range(len(hsp_1)):
        a1=hsp_1[ii]
        for jj in range(3):
            a1[jj]=np.float(a1[jj])
        hsp.append(a1)
    hsp=np.array(hsp)
    f.close()
    return hsp,hsp_label,node    

 
def read_poscar(pos_file):
    try:
        latt_vec=np.zeros([3,3])    
        poscar=open(pos_file)
        line=poscar.readline()
        sys_name=line
        latt= float(poscar.readline())
        for i in range(3):
        	line = poscar.readline().split()
        	for j in range(3):
        		latt_vec[i,j] = latt*float(line[j])
        atom_ele=poscar.readline().split()
        atom_num_old=poscar.readline().split()
        atom_num=[]
        for n in atom_num_old:
            atom_num.append(int(n)) 
        atom=atom_ele+atom_num
        atom_sum=(np.array(atom_num)).sum()
        line2=poscar.readline()
        S=np.zeros([atom_sum,3])
        for ii in range(atom_sum):
        	S1 = poscar.readline().split()
        	for jj in range(3):
        		S[ii,jj] =np.float(S1[jj])
    except IOError:
    	print ('POSCAR is not exist or file open failed! \n')
    poscar.close()
    return latt,latt_vec,atom,S,sys_name    
    
def read_recip(pos_file):    
    pos_file='POSCAR'
    read_poscar(pos_file)
    latt,latt_vec,atom,S,sys_name=read_poscar(pos_file)
    latt_rec = 2*math.pi*np.mat(latt_vec).I.T
    return latt_rec    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    #pos_file='POSCAR'
    #eigenval_filename='EIGENVAL'
    #dos_file='DOSCAR'
    #kpoint_filename='KPOINTS2'
    #engval,kpoints=read_eigenval(eigenval_filename)
    #E_fermi = get_fermi_from_doscar(dos_file)
    #if len(engval.shape)==2:
    #    energy_up=engval
    #    latt_rec=read_recip(pos_file)
    #    latt,latt_vec,atom,S,sys_name=read_poscar(pos_file)
    #    temp=[]
    #    for ii in range(1,kpoints.shape[0]):
    #        temp.append(np.linalg.norm((kpoints[ii-1]-kpoints[ii])*latt_rec))
    #    temp.insert(0,0)
    #    temp=np.cumsum(temp)
    #    temp=temp.reshape(kpoints.shape[0],1)
    #    kpoints=np.append(kpoints,temp,axis=1)
    ##    print(kpoints)
    #    hsp,hsp_label,node=read_high_sym_point(kpoint_filename)
    #    energy=np.zeros([engval.shape[0],int(engval.shape[1]/2)])
    #    for ii in range(energy.shape[0]): 
    #        for jj in range(energy.shape[1]):
    #            energy[ii][jj]=engval[ii][int(2*jj)]    
    #    pl
