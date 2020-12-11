import sys 
import csv 
import math 
import random 
from textwrap import dedent
import matplotlib.pyplot as plt 

'''
    Owais A. Shahzada 
    04/01/2019 
    Objective: K-Means CLustering on Gene Expression Data 
    Description: Using RNA Sequences samples after genome assembly and identifying which 
    genes are differentialy expressed you can group the genes which express similarly using 
    K-means clustering. This can identify which genes are "on" or "off" 

'''


def euclidean_distance(gene, centroid):
    '''
        @Parameters gene, centroid
        
        Each gene is assigned to a cluster based on how 
        close the gene is to cetnroid using the Euclidian 
        Distance formula. 
        
        return euclidean_distance : Float 
        
    '''
    # centroid = [mean_con, mean_h9, mean_h18]
    mean_con,mean_h9,mean_h18 = centroid
    gene_con = gene.con
    gene_h9 = gene.h9
    gene_h18 = gene.h18
    euclidean_distance = (gene_con - mean_con)**2 + (gene_h9 - mean_h9)**2 + (gene_h18 - mean_h18)**2
    return euclidean_distance

#centroids = {0: [x,y,z], 1:[x,y,z]...K:[]}
def reassign_genes_to_cluster(gene_list, centroid_dict, new_cluster):
    '''
        @Parameters gene_list, centroid_dict, new_cluster

        Each genes distance will be checked using the
        Euclidian Distance formula and resassigned to 
        the closest centroid. 
        
        return new_cluster : dict 
    '''
    
    for gene in gene_list:
        # Start with first gene from gene_list and first centroid from centroids_dict
        smallest_distance = euclidean_distance(gene, centroid_dict[0])
        centroid_index = 0
        for centroid_number in centroid_dict:
            distance = euclidean_distance(gene, centroid_dict[centroid_number])
            if distance < smallest_distance:
                smallest_distance = distance
                centroid_index = centroid_number
        new_cluster[centroid_index][gene.id] = gene 
    return new_cluster

def create_empty_cluster(cluster,k):
    for cluster_number in range(0,k):
            cluster[cluster_number] = {}
    return cluster

def centroids(k, clusters):
    '''
        @Parameters k, clusters 
        
        Add all cons in each cluster, and divide by total number of genes in the cluster 
        Add all h9s in each cluster, and divide by total number of genes in the cluster 
        Add all h18s in each cluster, and divide by total number of genes in the cluster 
        
        return centoids: Dictionary 
    '''
    centroids = {}
    for centroid_numbers in range(0,k):
        centroids[centroid_numbers] = []
    
    for cluster_number in clusters:
        length_list = len(clusters[cluster_number])
        if length_list == 0:
            length_list = 1 
        mean_con = sum(clusters[cluster_number][genes].con for genes in clusters[cluster_number])/length_list
        mean_h9 =  sum(clusters[cluster_number][genes].h9 for genes in clusters[cluster_number])/length_list
        mean_h18 = sum(clusters[cluster_number][genes].h18 for genes in clusters[cluster_number])/length_list
        
        centroids[cluster_number].append(mean_con)
        centroids[cluster_number].append(mean_h9)
        centroids[cluster_number].append(mean_h18)
    return centroids
            
class Gene:
    def __init__(self,id,con,h9,h18):
        self.id = id 
        self.con = con
        self.h9 = h9
        self.h18 = h18

def main(k,should_print):
    with open(sys.argv[1], 'r') as csv_file:
        k_input = k
        csv_reader = csv.reader(csv_file, delimiter=',')
        next(csv_reader)

        clusters = {}
        clusters = create_empty_cluster(clusters, k_input)
        
        gene_list = []

        for row in csv_reader:
            gene_id = row[0]
            con = row[3]
            h9 = row[4]
            h18 = row[5]
            genes = Gene(gene_id, float(con), float(h9), float(h18))
            gene_list.append(genes)
            random_num = random.randrange(0,k_input)
            clusters[random_num][gene_id] = genes
            '''
            clusters = {
                0: {
                    'gene4_id': Gene4,
                    'gene5_id': Gene5,
                    'gene6_id': Gene6
                },
                1: {
                    'gene8_id': Gene8,
                    'gene9_id': Gene9
                }
            }
            '''

        # Starting point
        centroids_dict = centroids(k_input,clusters)
        old_cluster = clusters
        iterations = 0
        while True:
            iterations += 1 
            new_cluster = {}
            new_cluster = create_empty_cluster(new_cluster,k_input)
            new_cluster = reassign_genes_to_cluster(gene_list, centroids_dict,new_cluster)
            if new_cluster == old_cluster:
                break
            else:
                centroids_dict = centroids(k_input, new_cluster)
                old_cluster = new_cluster
        if should_print is True:
            print(dedent(f"""
            K = {k_input}
            Completed in {iterations} iterations.
            """))       
        return centroids_dict, new_cluster 

        
def calc_sse():
    '''
        @Parameters None 
        
        Minimize the total within-cluster sum of squared error 
        and visualize it the graph.
        
        return None 
    
    '''
    
    all_sses = {}
    for k in range(2, 11):
        centroids_dict, new_cluster = main(k, False)
        sse = 0
        for cluster_number in new_cluster:
            for gene, gene_val in new_cluster[cluster_number].items():
                mean_con, mean_h9, mean_h18 = centroids_dict[cluster_number]
                sse += pow((gene_val.con - mean_con), 2)
                sse += pow((gene_val.h9 - mean_h9), 2)
                sse += pow((gene_val.h18 - mean_h18), 2)
        all_sses[k] = sse 
    fig, ax = plt.subplots()
    x_values = list(all_sses.keys())
    y_values = list(all_sses.values())
    ax.plot(x_values, y_values)
    ax.set(xlabel='Number of Clusters', ylabel='Total Sum Squared Error')
    # fig.savefig('shahzada_hw7.png')
    plt.show()

def random_color():
    return (random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1))

def plot_best_kmeans():
    centroids_dict, new_cluster = main(5, False)
    ax = plt.axes(projection='3d')
    for centroid_number in centroids_dict:
        color = random_color()
        max_con, max_h9, max_h18 = centroids_dict[centroid_number] 
        ax.scatter3D(max_con, max_h9, max_h18, color='black', marker='D')
        for gene_id, gene in new_cluster[centroid_number].items():
            ax.scatter3D(gene.con,gene.h9,gene.h18, color=color)
    plt.show()


if __name__ == '__main__':
    k_input = int(sys.argv[2])
    main(k_input, True)
    calc_sse()
    plot_best_kmeans() # Best K-mean 
    exit()
