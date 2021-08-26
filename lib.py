import pandas as pd
import numpy as np
import scipy

### Small class to structure the code
class ImagingExperiment():
    def __init__(self, dist_matrix, metadata):
        """
        ImagingExperiment class to store and convert between representations of multi-regions microscopy results.
        data_matrix: 2d numpy array, NxN
                    distance matrix, where the size N corresponds to number of marks studied. 
        """
        self.dist_matrix = dist_matrix.copy()
        self.l = len(dist_matrix)
        self.metadata = dict(metadata)
        
    def convert_to_contacts(self, threshold=150):
        """
        Convert microscopy distances to contacts using specified threshold.
        """
        self.contacts = np.array(np.where(self.dist_matrix<=threshold)).T
        
    def annotate_contacts(self, segments_table, one_based=True):
        """
        Annotate contacts by table of segments. 
        segments_table: pd.DataFrame with columns: ['chrom', 'start', 'end'] (assumed but not checked)
                        The numeration of the index should correspond to the contatcing pairs.
        one_based: if 1, it assumes that segments_table index starts from 1. This adds +1 to the numbers of contacts. 
        """
        shift = 1 if one_based else 0
        self.contacts_dataframe = pd.DataFrame([ 
            np.concatenate([segments_table.loc[x[0]+shift, :].values, 
                            segments_table.loc[x[1]+shift, :].values, [1]]) 
            for x in self.contacts], 
            columns=['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'count']).astype({
                    'chrom1':str, 'start1':int, 'end1':int, 
                    'chrom2':str, 'start2':int, 'end2':int, 
                    'count':int})
        
    def annotate_distances(self, segments_table, one_based=True):
        """
        Annotate distances by table of segments. 
        segments_table: pd.DataFrame with columns: ['chrom', 'start', 'end'] (assumed but not checked)
                        The numeration of the index should correspond to the contatcing pairs.
        one_based: if 1, it assumes that segments_table index starts from 1. This adds +1 to the numbers of contacts. 
        """
        shift = 1 if one_based else 0
        self.dist_dataframe = pd.DataFrame([ 
            np.concatenate( [segments_table.loc[i+shift, :].values, 
                             segments_table.loc[j+shift, :].values, 
                             [self.dist_matrix[i, j] ]] ) for i in range(self.l) for j in range(i, self.l)  ], 
            columns=['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'distance']).astype({
                    'chrom1':str, 'start1':int, 'end1':int, 
                    'chrom2':str, 'start2':int, 'end2':int, 
                    'distance':float})
        
    def save_contacts(self, filename_mask):
        """
        Store contacts in local tsv file with given filename. 
        filename_mask: string
                        Either a file name or name mask, where the metadata will be inserted.
        """
        self.contacts_dataframe.to_csv(filename_mask.format(**self.metadata), header=True, index=False, sep='\t')
        
    def save_distances(self, filename_mask):
        """
        Store distances in local tsv file with given filename. 
        filename_mask: string
                        Either a file name or name mask, where the metadata will be inserted.
        """
        self.dist_dataframe.to_csv(filename_mask.format(**self.metadata), header=True, index=False, sep='\t')