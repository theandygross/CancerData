#! /usr/bin/env python

import pickle as pickle

import pandas as pd
import numpy as np

from Stats.Scipy import anova

from Data.Containers import Dataset
from Helpers.LinAlg import frame_svd, extract_pc
from Helpers.Pandas import true_index, screen_feature
from Stats.Scipy import pearson_pandas

from Data.Intermediate import read_data


def exp_change(s):
    """
    Calculates an anova for the change in expression across a variable
    on the second level of a MultiIndex. (eg. tumor/normal).
    """
    return anova(pd.Series(s.index.get_level_values(1), s.index), s)


def extract_pc_filtered(df, pc_threshold=.2, filter_down=True,
                          tumor_code='01', normal_code='11'):
    """
    First pre-filters for patients with no tumor/normal change.
    Then normalizes by normals.
    """
    if (normal_code in df.columns.levels[1]) and filter_down:
        tt = df.xs('11', axis=1, level=1)
        rr = df.apply(exp_change, 1).sort('p')
        m, s = tt.mean(1), tt.std(1)
        df_n = df.xs(tumor_code, axis=1, level=1)
        df_n = ((df_n.T - m) / s).T
        df_n = df_n.ix[true_index(rr.p < .05)]
    else:  # No matched normals
        df_n = df.xs(tumor_code, axis=1, level=1)
        df_n = ((df_n.T - df_n.mean(1)) / df_n.std(1)).T
    pc = extract_pc(df_n, pc_threshold, standardize=False)
    return pc


def extract_geneset_pcs(df, gene_sets, filter_down=True):
    """
    Extract PCs for all gene sets.
    """
    pc = {p: extract_pc_filtered(df.ix[g].dropna(), .2, filter_down) for p, g,
          in gene_sets.iteritems()}
    pc = pd.DataFrame({p: s for p, s in pc.iteritems() if s}).T
    pat_vec = pd.DataFrame(pc.pat_vec.to_dict()).T
    return pc.gene_vec, pc.pct_var, pat_vec


def extract_features(df, low_expr=-1, var_cutoff=.25):
    """
    Extract features from a data-frame of expression data.

    Binary markers are used when expression is only present
    (having more than 0.5 reads per million) in a moderate
    fraction of the cohort (between 20 cases and half of the
    cohort).

    Real-valued gene and miRNA expression levels were used
    for differentially expressed features not assigned as
    binary markers.

    Both binary and real valued markers are filtered based
    on differential expression with the adjacent normal
    tissue that is available for some samples.
    """
    df_n = df.xs('01', level=1, axis=1)
    binary = df_n > low_expr
    binary = binary[binary.sum(1).isin(range(20, df.shape[1] / 2))]
    rr = df.ix[binary.index].apply(exp_change, 1)
    binary = binary.ix[true_index(rr.p < .05)]
    
    real = df_n.ix[df_n.index.diff(binary.index)]
    singles = real[((real.max(1) - real.min(1)) > 1)]
    singles = singles[(singles.std(1) > var_cutoff)]
    ch = df.ix[singles.index].apply(exp_change, 1)
    singles = df_n.ix[true_index(ch.p < .01)]
    return binary, singles, real


class RealDataset(Dataset):
    """
    Inherits from Dataset class.  Adds some added processing for real valued
    data.
    """
    def __init__(self, run, cancer, data_type, patients=None, drop_pc1=False,
                 create_real_features=True, create_meta_features=True,
                 filter_down=True):
        """
        """
        Dataset.__init__(self, cancer.path, data_type, compressed=True)
        self.df = read_data(run.data_path, cancer.name, data_type,
                               tissue_code='All')
        if patients is not None:
            self.df = self.df.ix[:, patients].dropna(axis=1, how='all')
            self.patients = patients
        else:
            self.patients = self.df.xs('01', 1, 1).columns
        
        self.global_vars = pd.DataFrame(index=self.patients)
        self.features = {}
        self.global_loadings = pd.DataFrame(index=self.df.index)
        self._calc_global_pcs(drop_pc1)
        
        if create_real_features is True:
            self._get_real_features()
        
        if create_meta_features is True:
            self._get_meta_features(run.gene_sets, filter_down)
            
        self.features = pd.concat(self.features)
            
    def _get_real_features(self):
        binary, singles, real = extract_features(self.df)
        background_df = real.ix[real.index.diff(singles.index)].dropna()
        background = extract_pc(background_df, 0)
        ss = screen_feature(background['pat_vec'], pearson_pandas, singles)
        singles = singles.ix[ss.p > 10e-5]
        
        singles = ((singles.T - singles.mean(1)) / singles.std(1)).T
        U, S, pc = frame_svd(singles)
        
        self.features['binary'] = binary
        self.features['real'] = singles
        self.global_vars['background'] = background['pat_vec']
        self.global_vars['filtered_pc1'] = pc[0]
        self.global_vars['filtered_pc2'] = pc[1]
        self.global_loadings['background'] = background['gene_vec']
        self.global_loadings['filtered_pc1'] = U[0]
        self.global_loadings['filtered_pc2'] = U[1]
        
    def _get_meta_features(self, gene_sets, filter_down):
        gs = extract_geneset_pcs(self.df, gene_sets, filter_down)
        self.loadings, self.pct_var, pathways = gs
        if hasattr(self.global_vars, 'background'):
            r = screen_feature(self.global_vars.background, pearson_pandas,
                               pathways)
            pathways = pathways.ix[r.p > 10e-5]
        pathways = ((pathways.T - pathways.mean(1)) / pathways.std(1)).T
        U, S, pc = frame_svd(pathways)
        
        self.pathways = pathways
        self.features['pathways'] = pathways
        self.global_vars['pathway_pc1'] = pc[0]
        self.global_vars['pathway_pc2'] = pc[1]
        self.global_loadings['pathway_pc1'] = U[0]
        self.global_loadings['pathway_pc2'] = U[1]

    def _calc_global_pcs(self, drop_pc1=False):
        """
        Normalize data and calculate principal components. If drop_pc1 is
        set to True, also reconstructs the normalized data without the
        first PC. 
        """
        df = self.df.xs('01', axis=1, level=1)
        norm = ((df.T - df.mean(1)) / df.std(1)).T
        U, S, vH = frame_svd(norm)
        self.global_vars['pc1'] = vH[0]
        self.global_vars['pc2'] = vH[1]
        self.global_loadings['pc1'] = U[0]
        self.global_loadings['pc2'] = U[1]        
        if drop_pc1 is True:
            S_n = S.copy()
            S_n[0] = 0
            norm = U.dot(pd.DataFrame(np.diag(S_n)).dot(vH.T))
            
        return norm
            
        
def initialize_real(cancer_type, report_path, data_type, patients=None,
                    drop_pc1=False, create_real_features=True,
                    create_meta_features=True, filter_down=False,
                    save=True):
    """
    Initialize real-valued data for down-stream analysis.
    """
    run = pickle.load(open(report_path + '/RunObject.p', 'rb'))
    cancer = run.load_cancer(cancer_type)
    
    if data_type is 'miRNASeq':
        create_meta_features = False
    
    data = RealDataset(run, cancer, data_type, patients, drop_pc1,
                       create_real_features, create_meta_features, filter_down)

    if save is True:
        data.save()
    return data
