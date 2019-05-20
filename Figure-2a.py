#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 14 18:17:45 2018

@author: seyed
"""
import numpy as np
np.set_printoptions(threshold=np.nan)
import math
from scipy.stats import linregress
import matplotlib.pyplot as plt
import pandas as pd
import random



cell_genos=np.linspace(1,2,10)
st_dev_cell=0.25

st_dev_group=0.25

st_dev_cluster=0.25
number_of_stdevs=5 #this name is a bit confusing, but this is the number of different st_devs to be run in the iteration
replicates=10

#use these lists for running multiple iterations of cell and group stdevs
group_stdev_iterator=np.linspace(.0001,st_dev_cell,number_of_stdevs)
cell_stdev_iterator=np.linspace(.0001,st_dev_group,number_of_stdevs)
cluster_stdv_iterator=np.linspace(.0001,st_dev_cluster,number_of_stdevs)

group_sizes=[]
st_dev_cells=[]
st_dev_groups=[]
slope_cell=[]
slope_group_volume=[]
slope_group_radius=[]
slope_group_settling=[]
stdev_list = [ [] for i in range(number_of_stdevs) ]	
group_sizes_list = [ [] for i in range(len(cell_genos)) ]
#print(group_sizes_list)
time_step=6 #number of generation
gs=32  #is the average number of cells per group from 2 to the assigned value
mean9=[]
mean2=[]
x_1=[-1,34]
y_1=[1,1]

fig = plt.figure()


for ijk in [1,2]:
    stdev_list = [ [] for i in range(number_of_stdevs) ]	
    group_sizes_list = [ [] for i in range(len(cell_genos)) ]
    if ijk == 1:
        st_dev_cluster=0.000001
    if ijk == 2:
        st_dev_cluster=0.25
    plots=np.array([[],[],[],[],[],[],[],[],[],[]])
    plotg=np.array([[],[],[],[],[],[],[],[],[],[]])
    for ii in range(0,len(stdev_list)):
        st_dev_cell=cell_stdev_iterator[ii]
        #st_dev_group=group_stdev_iterator[ii]
        #st_dev_cluster=cluster_stdv_iterator[ii]
        for a in np.arange(2,gs+2,2):
            trigger = 0 
            cluster_size=a
            cluster_std=a*st_dev_cluster
            for b in range(0,10):
                cell_x=np.array([])
                cell_y=np.array([])
                pop=np.zeros((len(cell_genos)*cluster_size,10))
                np.set_printoptions(threshold=np.nan)
                for i in range(0,np.shape(pop)[0]):
                    pop[i][0]=i #Number each cell
                    pop[i][1]=cell_genos[math.floor(i/cluster_size)]
                    pop[i][2]=np.random.normal(pop[i][1],st_dev_cell)
                    pop[i][4]=math.floor(i/cluster_size)
                curr_gen=pop
                    
                for j in range(1,time_step+1): #looping through generations
                    cluster_IDs=np.zeros(2**j)
                    k=j-1
        
                    pre_gen_CI=np.zeros(2**k)
                    for i in range(len(cluster_IDs)):
                        cluster_IDs[i]=len(cluster_IDs)-1+i
                    for i in range(2**(j-1)):
                        pre_gen_CI[i]=len(pre_gen_CI)-1+i
                    nex_gen=np.zeros([1,10])
                    per_gen=len(cell_genos)*cluster_size*2**j
                    cluster_variance_factor=np.random.normal(1,st_dev_group)
                    cell_max=int(max(curr_gen[:,0]))+1
                    cluster_counter=int(max(curr_gen[:,4]))+1
                    y=cell_max
                    next_gen=np.zeros([1,10])     
                    num_clusters=len(cell_genos)*2**j #total number of clusters for each generation
                    par_clusters_ID=[[] for i in range(len(cell_genos))]
                                         
                    for iii in range(num_clusters): 
                        cluster_variance_factor=np.random.normal(1,st_dev_group)
                        extra=np.mod(cluster_counter+iii,len(cell_genos))
                        xx=(cluster_counter+iii-extra)/len(cell_genos)
                        for k in range(len(cluster_IDs)):
                            if xx == cluster_IDs[k]:
                                par_clu_co=pre_gen_CI[int(k/2)] #parents cluster counter
                        par_cluster_id=par_clu_co*10+extra 
                            
                        par=np.zeros([1,10])
                        cs=int(np.random.normal(cluster_size+0.455,cluster_std))
                            #cs=cluster_size
                            #print(cs)
                        if cs <= 0:
                            cs=1
                        cluster=np.zeros([cs,10])
                        for ij in range(len(curr_gen)):
                            if par_cluster_id == curr_gen[ij][4]:
                                par=np.vstack([par,curr_gen[ij]])
                        part_can=par[1:] #parent candidates
                        aux0=np.arange(len(part_can))
                        for n in range(len(cluster)):
                                
                            aux1=random.choice(aux0)
                            parent=part_can[aux1]
                            if len(aux0) > 1:
                                aux0=np.setdiff1d(aux0,np.array([aux1]))
                            else:
                                aux0=np.arange(len(part_can))
                            cluster[n][0]=y+n
                            cluster[n][1]=cell_genos[np.mod(cluster_counter+iii,10)]
                            cluster[n][2]=np.random.normal(cluster[n][1],st_dev_cell)*cluster_variance_factor
                            cluster[n][3]=parent[0]
                            cluster[n][4]=cluster_counter+iii
                            cluster[n][5]=parent[4]
                            cluster[n][6]=0
                            cluster[n][7]=0
                            cluster[n][8]=0
                            cluster[n][9]=parent[2]
                                
                        y=len(cluster)+y
                        cell_max=int(max(curr_gen[:,0]))+1+len(cluster)
                        nex_gen=np.vstack([nex_gen,cluster])
                    curr_gen=nex_gen[1:]
                    pop=np.vstack([pop,curr_gen])
                    ng=nex_gen[1:]
               
                np.savetxt("full-population-test.csv", pop, delimiter=",")
                cell_x=pop[:,9]
                cell_y=pop[:,2]
                cell_x=cell_x[len(cell_genos)*cluster_size+1:]
                cell_y=cell_y[len(cell_genos)*cluster_size+1:]
                mean9.append(np.mean(cell_x))
                mean2.append(np.mean(cell_y))
                df=pd.DataFrame(pop)
                size_by_ID=df.groupby(4)[2].sum()
                parent_by_ID=df.groupby(4)[5].mean()
                joined=pd.concat([size_by_ID,parent_by_ID], axis=1, ignore_index=True)    
                parent_size=[]
                for i in range(0,len(joined[0])):
                    j=joined[1][i]
                    parent_size.append(joined[0][j])
                offspring_size=joined[0]
                parent_size_cleaned=list(parent_size[len(cell_genos):])
                offspring_size_cleaned=list(offspring_size[len(cell_genos):])
                tempratio=(linregress(parent_size_cleaned,offspring_size_cleaned)[0]) / (linregress(cell_x,cell_y)[0])   
                #if tempratio < .9 or tempratio > 2.5:
                    #print('tempratio', tempratio, ' st_dev_cell=', st_dev_cell, ' cluster size=', cs)
                stdev_list[ii].append(tempratio)
                group_sizes_list[ii].append(a)        
    if ijk == 1:
        np.savetxt("A_group_sizes_list_low_stdv.csv", np.array(group_sizes_list), fmt='%5s', delimiter=",")
        np.savetxt("A_stdev_list_list_low_stdv.csv", np.array(stdev_list), fmt='%5s', delimiter=",")
        ax = fig.add_subplot(1,2,ijk)
        #plt.subplot(1,2,ijk)
        ax.plot(x_1,y_1,'--',color='k')
        ax.set_title("$\sigma''=10^{-4}$")
        #plt.plot(x_1,y_1,'--',color='k')
        cmap = plt.get_cmap('rainbow')
        colors = [cmap(i) for i in np.linspace(0, 1, len(stdev_list))]
        #plt.xlim(1,gs+1)
        ax.set_xlim(1,gs+1)
        #ax.set_aspect('equal')
        for i, color in enumerate(colors, start=0):
            ax.scatter(group_sizes_list[i],stdev_list[i], color=color, alpha=.5)
        #ax.set_xlabel('Mean number of cells per group')
        ax.set_ylabel('Ratio of group to cellular-level heritability for size', fontsize='10')
        ax.set_ylim(0,1.8)
    else:
        np.savetxt("A_group_sizes_list_high_stdv.csv", group_sizes_list, fmt='%5s', delimiter=",")
        np.savetxt("A_stdev_list_list_high_stdv.csv", stdev_list, fmt='%5s', delimiter=",")
        ax = fig.add_subplot(1,2,ijk)
        ax.set_title("$\sigma''=0.25$")
        ax.set_xlim(1,gs+1)
        ax.plot(x_1,y_1,'--',color='k')
        cmap = plt.get_cmap('rainbow')
        colors = [cmap(i) for i in np.linspace(0, 1, len(stdev_list))]      
        for i, color in enumerate(colors, start=0):
            ax.scatter(group_sizes_list[i],stdev_list[i], color=color, alpha=.5)
        #ax.set_aspect('equal')
        ax.set_ylim(0,1.8)

        #plt.xlabel('Mean number of cells per group')
        #plt.ylabel('Ratio of group to cellular-level heritability for size')
fig.text(0.5, 0.04, 'Mean number of cells per group', ha='center', va='center', fontsize='10')
plt.show()
plt.savefig("A_Ratio of group to cell heritability iterator=cell variance, group sd=%s.png" %(st_dev_group), dpi=300, orientation='landscape')
plt.savefig("A_Ratio of group to cell heritability iterator=cell variance, group sd=%s.pdf" %(st_dev_group), dpi=300, orientation='landscape')
#print(max(pop[:,0]))
print("Col 9", np.mean(mean9))
print("Col 2", np.mean(mean2))