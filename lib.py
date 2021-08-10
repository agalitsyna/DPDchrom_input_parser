import pandas as pd
import numpy as np
import scipy

### Small class to structure the code
class ImagingExperiment():
    def __init__(self, dist_matrix, metadata):
        self.dist_matrix = dist_matrix.copy()
        self.l = len(dist_matrix)
        self.metadata = dict(metadata)
        
    def convert_to_contacts(self, threshold=150):
        self.contacts = np.array(np.where(self.dist_matrix<=threshold)).T
        
    def annotate_contacts(self, segments_table, one_based=True):
        self.contacts_dataframe = pd.DataFrame([ 
            np.concatenate([segments_table.loc[x[0]+1 if one_based else 0, :].values, 
                            segments_table.loc[x[1]+1 if one_based else 0, :].values, [1]]) 
            for x in self.contacts], 
            columns=['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'count']).astype({
                    'chrom1':str, 'start1':int, 'end1':int, 
                    'chrom2':str, 'start2':int, 'end2':int, 
                    'count':int})
        
    def annotate_distances(self, segments_table, one_based=True):
        self.dist_dataframe = pd.DataFrame([ 
            np.concatenate( [segments_table.loc[i+1 if one_based else 0, :].values, 
                             segments_table.loc[j+1 if one_based else 0, :].values, 
                             [self.dist_matrix[i, j] ]] ) for i in range(self.l) for j in range(i, self.l)  ], 
            columns=['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'distance']).astype({
                    'chrom1':str, 'start1':int, 'end1':int, 
                    'chrom2':str, 'start2':int, 'end2':int, 
                    'distance':float})
        
    def save_contacts(self, filename_mask):
        self.contacts_dataframe.to_csv(filename_mask.format(**self.metadata), header=True, index=False, sep='\t')
        
    def save_distances(self, filename_mask):
        self.dist_dataframe.to_csv(filename_mask.format(**self.metadata), header=True, index=False, sep='\t')