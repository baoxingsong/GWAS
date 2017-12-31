"""
Created: 25/12/2017  00:40:39
Author: Baoxing Song
Email: songbaoxing168@163.com
This source code is partially adopted from https://github.com/bvilhjal/mixmogam/
And the R source code of emma is referred
"""

import numpy as np
from numpy import linalg
from scipy import stats
from scipy import optimize
import warnings


class phenotype_data:
    """
    A class that encapsulates phenotype values and provides basic functionality for these.
    """

    def __init__(self, phen_dict=None, phenotype_names=[], phen_ids=None):
        if phen_dict:
            self.phen_dict = phen_dict
            self.phen_ids = phen_dict.keys()
            for pid in phen_dict:
                self.phen_dict[pid]['transformation'] = None
                self.phen_dict[pid]['raw_values'] = []
        else:
            if phen_ids:
                self.phen_ids = phen_ids
            else:
                self.phen_ids = range(len(phenotype_names))
            self.phen_dict = {}
            for i, pid in enumerate(self.phen_ids):
                self.phen_dict[pid] = {'name': phenotype_names[i], 'ecotypes': [], 'values': [], 'transformation': None,
                                       'raw_values': []}

    def get_values(self, pid):
        return self.phen_dict[pid]['values']


class _SnpsData_(object):
    """
    An abstract superclass.
    """

    def __init__(self, snps, positions, baseScale=None, accessions=None, arrayIds=None, chromosome=None,
                 alignment_positions=None, id=None, marker_types=None, alphabet=None, associated_positions=None):
        self.snps = snps  # list[position_index][accession_index]
        self.positions = positions  # list[position_index]
        if accessions:
            self.accessions = accessions  # list[accession_index]
            # self._convert_to_tg_ecotypes_()
        if arrayIds:
            self.arrayIds = arrayIds  # list[accession_index]
        self.chromosome = chromosome
        self.alignment_positions = alignment_positions
        self.id = id
        self.marker_types = marker_types  # Where do these markers come frome, what type are they?  Useful for later analysis.
        self.alphabet = None
        self.missingVal = None
        self.associated_positions = associated_positions


class SNPsData(_SnpsData_):
    """
    An alternative to the old SnpsData class, where this uses scipy to speed things up when possible.

    """
    alphabet = [-1, 0, 1, 2, 3]  # Here -1 is thought to be missing data value.
    def __init__(self, snps, positions, accessions=None, arrayIds=None, chromosome=None,
                 alignment_positions=None, id=None, marker_types=None, missing_val=-1):
        self.snps = snps
        self.positions = positions
        self.accessions = accessions
        self.arrayIds = arrayIds
        self.chromosome = chromosome
        self.alignment_positions = alignment_positions
        self.marker_types = marker_types  # Where do these markers come frome, what type are they?  Useful for later analysis.
        self.id = id
        self.missingVal = missing_val


class SNPsDataSet:
    """
    A class that encompasses multiple _SnpsData_ chromosomes objects (chromosomes), and can deal with them as a whole.

    This object should eventually replace the snpsdata lists..
    """
    snpsDataList = None
    chromosomes = None
    accessions = None
    def __init__(self, snpsds, chromosomes, id=None, call_method=None, data_format=None):
        self.snpsDataList = snpsds
        self.chromosomes = chromosomes
        self.accessions = self.snpsDataList[0].accessions
        self.array_ids = self.snpsDataList[0].arrayIds
        self.id = id
        self.missing_val = snpsds[0].missingVal
        self.call_method = call_method
        self.data_format = data_format  # binary, diploid_ints, floats, int
        if not id and snpsds[0].id:
            self.id = id
        for i in range(1, len(self.chromosomes)):
            if self.accessions != self.snpsDataList[i].accessions:
                raise Exception("Accessions (or order) are different between SNPs datas")
        self.is_binary = list(snpsds[0].snps[0]).count(0) or list(snpsds[0].snps[0]).count(1)

    def get_snps(self, random_fraction=None, cache=False):
        if cache:
            try:
                return self.snps
            except Exception:
                pass
        snplist = []
        if random_fraction:
            import random
            for snpsd in self.snpsDataList:
                for snp in snpsd.snps:
                    if random.random() < random_fraction:
                        snplist.append(snp)
        else:
            for snpsd in self.snpsDataList:
                snplist.extend(snpsd.snps)
                if cache:
                    self.snps = snplist
        return snplist


def scale_k(k, verbose=False):
    c = np.sum((np.eye(len(k)) - (1.0 / len(k)) * np.ones(k.shape)) * np.array(k))
    scalar = (len(k) - 1) / c
    if verbose:
        print ('Kinship scaled by: %0.4f' % scalar)
    k = scalar * k
    return k



class LinearMixedModel:
    """
    A class for linear mixed models
    """
    def __init__(self, Y=None, dtype='single'):
        """
        The fixed effects should be a list of fixed effect lists (SNPs)
        """

        self.n = len(Y)

        self.Y = np.matrix(Y, dtype=dtype)
        residuals = self.Y - np.mean(self.Y)
        self.y_var = residuals * residuals.T
        
        self.Y.shape = (self.n, 1) 
        self.X = np.matrix(np.ones((self.n, 1), dtype=dtype))  # The intercept
        self.p = 1
        self.beta_est = None

        # A list of random effect type, and the cov matrix.
        self.random_effects = [('normal', np.matrix(np.identity(self.n)))]  # The first random effect is the IID error.

    def add_factor(self, x, lin_depend_thres=1e-4):
        """
        Adds an explanatory variable to the X matrix.
        """
        # Checking whether this new cofactor in linearly independent.
        new_x = np.array(x)
        new_x.shape = len(x)
        (beta, rss, rank, sigma) = linalg.lstsq(self.X, new_x)
        if float(rss) < lin_depend_thres:
            warnings.warn(
                'A factor was found to be linearly dependent on the factors already in the X matrix.  Hence skipping it!')
            return False
        new_x.shape = (self.n, 1)
        self.X = np.hstack([self.X, new_x])
        self.p += 1
        return True

    def clear_factors(self, include_intercept=True):
        """
        remove all factors
        """
        self.p = 0
        if include_intercept:
            self.X = np.mat(np.repeat(1, self.n), dtype='single').T  # The intercept
            self.p = 1

    def add_random_effect(self, cov_matrix=None, effect_type='normal'):
        if effect_type != 'normal':
            raise Exception('Currently, only Normal random effects are allowed.')
        self.random_effects.append((effect_type, scale_k(cov_matrix)))

    def _get_eigen_L_(self, K=None, dtype='single'):
        if K is None:
            K = self.random_effects[1][1]
        if np.__version__ < '0.8':
            K = np.mat(K, dtype=dtype)
        evals, evecs = linalg.eigh(K)
        evals = np.array(evals, dtype=dtype)
        return {'values':evals, 'vectors':np.mat(evecs, dtype=dtype).T}

    def _get_eigen_R_(self, X=None, K=None, hat_matrix=None, dtype='single'):
        if X is None:
            X = self.X
        q = X.shape[1]
        if not hat_matrix:
            X_squared_inverse = linalg.pinv(X.T * X)  # (X.T*X).I
            hat_matrix = X * X_squared_inverse * X.T
        if K is None:
            K = self.random_effects[1][1]
        S = np.mat(np.identity(self.n)) - hat_matrix  # S=I-X(X'X)^{-1}X'
        M = np.mat(S * (K + self.random_effects[0][1]) * S, dtype='double')
        if np.__version__ < '0.8':
            M = np.mat(M, dtype=dtype)
        evals, evecs = linalg.eigh(M)  # eigen of S(K+I)S
        eig_values = np.array(evals[q:], dtype=dtype) - 1  # # Don't know what minute one, the emma R code did this
        return {'values':eig_values, 'vectors':(np.mat(evecs, dtype=dtype).T[q:])}

    def _rell_(self, delta, eig_vals, sq_etas):
        num_eig_vals = len(eig_vals)
        c_1 = 0.5 * num_eig_vals * (np.log(num_eig_vals / (2.0 * np.pi)) - 1)
        v = eig_vals + delta
        res = c_1 - 0.5 * (num_eig_vals * np.log(np.sum(sq_etas.flatten() / v)) + np.sum(np.log(v)))
        return res  # log-likelihoods (eq. 7 from paper)


    def _redll_(self, delta, eig_vals, sq_etas):
        num_eig_vals = len(eig_vals)
        v1 = eig_vals + delta
        v2 = sq_etas.flatten() / v1
        res = (num_eig_vals * np.sum(v2 / v1) / np.sum(v2) - np.sum(1.0 / v1))
        return res  # diffrentiated log-likelihoods (*2) (eq. 9 from paper)


    def _ll_(self, delta, eig_vals, eig_vals_L, sq_etas):
        n = self.n
        c_1 = 0.5 * n * (np.log(n / (2.0 * np.pi)) - 1)
        v1 = eig_vals + delta
        v2 = eig_vals_L + delta
        res = c_1 - 0.5 * (n * np.log(np.sum(sq_etas.flatten() / v1)) + np.sum(np.log(v2)))
        return res  # log-likelihoods (eq. 6 from paper)


    def _dll_(self, delta, eig_vals, eig_vals_L, sq_etas):
        num_eig_vals = len(eig_vals)
        v1 = eig_vals + delta
        v2 = sq_etas.flatten() / v1
        v3 = eig_vals_L + delta
        res = (self.n * np.sum(v2 / v1) / np.sum(v2) - np.sum(1.0 / v3))
        return res  # diffrentiated log-likelihoods (*2) (eq. 8 from paper)


    def get_estimates(self, eig_L=None, xs=None, ngrids=50, llim= -10, ulim=10, esp=1e-6,
                return_pvalue=False, return_f_stat=False, method='REML', verbose=False,
                dtype='single', eig_R=None, rss_0=None):
        """
        Get ML/REML estimates for the effect sizes, as well as the random effect contributions.
        Using the EMMA algorithm (Kang et al., Genetics, 2008).
        
        Methods available are 'REML', and 'ML'        
        """
        if verbose:
            print ('Retrieving %s variance estimates' % method)
        if xs is not None:
            X = np.hstack([self.X, xs])
        else:
            X = self.X
        K=self.random_effects[1][1]
        if not (eig_R and xs is not None):
            eig_R = self._get_eigen_R_(X=X, K=K)
        q = X.shape[1]  # number of fixed effects
        n = self.n  # number of individuls
        p = n - q
        m = ngrids + 1

        etas = np.array(eig_R['vectors'] * self.Y, dtype=dtype)
        sq_etas = etas * etas
        log_deltas = (np.arange(m, dtype=dtype) / ngrids) * (ulim - llim) + llim  # a list of deltas to search
        deltas = np.exp(log_deltas)
        assert len(deltas) == m, 'Number of deltas is incorrect.'
        eig_vals = np.array(eig_R['values'], dtype=dtype)
        assert len(eig_vals) == p, 'Number of eigenvalues is incorrect.'

        lambdas = np.reshape(np.repeat(eig_vals, m), (p, m)) + np.reshape(np.repeat(deltas, p), (m, p)).T
        s1 = np.sum(sq_etas / lambdas, axis=0)
        if method == 'REML':
            if verbose: print ('Calculating grid REML values')
            s2 = np.sum(np.log(lambdas), axis=0)
            lls = 0.5 * (p * (np.log((p) / (2.0 * np.pi)) - 1 - np.log(s1)) - s2)
            s3 = np.sum(sq_etas / (lambdas * lambdas), axis=0)
            s4 = np.sum(1 / lambdas, axis=0)
            dlls = 0.5 * (p * s3 / s1 - s4)
        elif method == 'ML':
            if verbose: print ('Calculating grid ML values')
            # Xis < -matrix(eig.L$values, n, m) + matrix(delta, n, m, byrow=TRUE)
            eig_vals_L = np.array(eig_L['values'], dtype=dtype)
            xis = np.reshape(np.repeat(eig_vals_L, m), (n, m)) + \
                np.reshape(np.repeat(deltas, n), (m, n)).T
            # LL <- 0.5*(n*(log(n/(2*pi))-1-log(colSums(Etasq/Lambdas)))-colSums(log(Xis)))    
            # dLL <- 0.5*delta*(n*colSums(Etasq/(Lambdas*Lambdas))/colSums(Etasq/Lambdas)-colSums(1/Xis))    

            s2 = np.sum(np.log(xis), axis=0)
            lls = 0.5 * (n * (np.log((n) / (2.0 * np.pi)) - 1 - np.log(s1)) - s2)
            s3 = np.sum(sq_etas / (lambdas * lambdas), axis=0)
            s4 = np.sum(1 / xis, axis=0)
            dlls = 0.5 * (n * s3 / s1 - s4)

        max_ll_i = np.argmax(lls)
        max_ll = lls[max_ll_i]

        last_dll = dlls[0]
        last_ll = lls[0]
        zero_intervals = []
        for i in range(1, len(dlls)):
            if dlls[i] < 0 and last_dll > 0:
                zero_intervals.append(((lls[i] + last_ll) * 0.5, i))
            last_ll = lls[i]
            last_dll = dlls[i]

        if len(zero_intervals) > 0:
            if verbose: print ('Found a zero interval... now performing Newton-Rhapson alg.')
            opt_ll, opt_i = max(zero_intervals)
            opt_delta = 0.5 * (deltas[opt_i - 1] + deltas[opt_i])
            # Newton-Raphson
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    if method == 'REML':
                        new_opt_delta = optimize.newton(self._redll_, opt_delta, args=(eig_vals, sq_etas), tol=esp,
                                                        maxiter=100)
                    elif method == 'ML':
                        new_opt_delta = optimize.newton(self._dll_, opt_delta, args=(eig_vals, eig_vals_L, sq_etas),
                                                        tol=esp, maxiter=100)
            except Exception:
                if verbose:
                    print ('Problems with Newton-Raphson method.')
                    print ("Using the maximum grid value instead.")
                    print ('opt_i:', opt_i)
                    print ('Grid optimal delta:', opt_delta)
                new_opt_delta = opt_delta
            # Validating the delta
            if opt_i > 1 and deltas[opt_i - 1] - esp < new_opt_delta < deltas[opt_i] + esp:
                opt_delta = new_opt_delta
                opt_ll = self._rell_(opt_delta, eig_vals, sq_etas)
            # Cheking lower boundary
            elif opt_i == 1 and 0.0 < new_opt_delta < deltas[opt_i] + esp:
                opt_delta = new_opt_delta
                opt_ll = self._rell_(opt_delta, eig_vals, sq_etas)
            # Cheking upper boundary
            elif opt_i == len(deltas) - 1 and new_opt_delta > deltas[opt_i - 1] - esp \
                        and not np.isinf(new_opt_delta):
                opt_delta = new_opt_delta
                opt_ll = self._rell_(opt_delta, eig_vals, sq_etas)
            else:
                if verbose:
                    print ('Local maximum outside of suggested possible areas?')
                    print ('opt_i:', opt_i)
                    print ('Grid optimal delta:', opt_delta)
                    print ("Newton's optimal delta:", new_opt_delta)
                    print ('Using the maximum grid value instead.')

            if verbose: print ('Done with Newton-Rahpson')
            if method == 'REML':
                opt_ll = self._rell_(opt_delta, eig_vals, sq_etas)
            elif method == 'ML':
                opt_ll = self._ll_(opt_delta, eig_vals, eig_vals_L, sq_etas)

            if opt_ll < max_ll:
                opt_delta = deltas[max_ll_i]
        else:
            if verbose: print ('No zero-interval was found, taking the maximum grid value.')
            opt_delta = deltas[max_ll_i]
            opt_ll = max_ll

        if verbose: print ('Finishing up.. calculating H_sqrt_inv.')
        l = sq_etas / (eig_vals + opt_delta)
        opt_vg = np.sum(l) / p  # vg
        opt_ve = opt_vg * opt_delta  # ve

        H_sqrt_inv = np.mat(np.diag(1.0 / np.sqrt(eig_L['values'] + opt_delta)), dtype=dtype) * eig_L['vectors']
        # V = opt_vg * K + opt_ve * sp.eye(len(K))
        # H_sqrt = cholesky(V).T
        # H_sqrt_inv = H_sqrt.I
        X_t = H_sqrt_inv * X
        Y_t = H_sqrt_inv * self.Y
        (beta_est, mahalanobis_rss, rank, sigma) = linalg.lstsq(X_t, Y_t)
        x_beta = X * beta_est
        residuals = self.Y - x_beta
        rss = residuals.T * residuals
        # x_beta_var = sp.var(x_beta, ddof=1)
        # var_perc = x_beta_var / self.y_var
        res_dict = {'max_ll':opt_ll, 'delta':opt_delta, 'beta':beta_est, 've':opt_ve, 'vg':opt_vg,
            'rss':rss, 'mahalanobis_rss':mahalanobis_rss, 'H_sqrt_inv':H_sqrt_inv,
            'pseudo_heritability':1.0 / (1 + opt_delta)}
        if return_pvalue:
            h0_X = H_sqrt_inv * self.X
            (h0_betas, h0_rss, h0_rank, h0_s) = linalg.lstsq(h0_X, Y_t)
            f_stat = (h0_rss / mahalanobis_rss - 1) * p / xs.shape[1]
            res_dict['var_perc'] = 1.0 - mahalanobis_rss / h0_rss
            res_dict['f_stat'] = float(f_stat)
            p_val = stats.f.sf(f_stat, (xs.shape[1]), p)
            res_dict['p_val'] = float(p_val)
        return res_dict  # , lls, dlls, sp.log(deltas)

    def variance_explained_by_significant_loci(self, sig_snps, K=None, ngrids=50, llim=-4, ulim=10, esp=1e-6, eig_L=None, includedSnpThreadshold = 1e-4, lin_depend_thres=1e-4):

        assert len(self.random_effects) == 2, "Expedited REMLE only works when we have exactly two random effects."
        if K is None:
            K = self.random_effects[1][1]
        if eig_L is None:
            eig_L = self._get_eigen_L_(K)

        # Run the model....

        for snp in sig_snps: #
            res = self.get_estimates(eig_L=eig_L, xs=np.matrix(snp).T, ngrids=ngrids, llim=llim, ulim=ulim,
                                     esp=esp, return_pvalue=True, return_f_stat=True)
            if res['p_val'] < includedSnpThreadshold: # if the variants could explain significantly proportion variance
                self.add_factor(snp, lin_depend_thres=lin_depend_thres) # linearly dependence will be checked

        res = self.get_estimates(eig_L=eig_L, xs=None, ngrids=ngrids, llim=llim, ulim=ulim,
                                 esp=esp, return_pvalue=False, return_f_stat=True)

        explained_variance_proportion = (self.y_var-res['rss'])/(self.y_var)
        return explained_variance_proportion


def parse_plink_tped_file(file_prefix, imputation_type='simple', return_kinship=False):
    """
    Requires a .tped file in 12 format.

    - Converts (on-the-fly) to a integer format.
    - Imputes missing data.
    """
    tped_filename = file_prefix + '.tped'
    tped_pickled_filename = tped_filename + '.imputed.pickled'
    tfam_filename = file_prefix + '.tfam'
    tfam_pickled_filename = tfam_filename + '.pickled'

    individs = []
    sex_list = []
    with open(tfam_filename) as f:
        for line in f:
            l = list(map(str.strip, line.split()))
            individs.append(l[1])
            sex_list.append(int(l[4]))
    num_individs = len(individs)

    #    k_mat = sp.zeros((num_individs, num_individs))

    chrom_pos_snp_dict = {}
    with open(tped_filename) as f:
        cur_chrom = -1
        for line_i, line in enumerate(f):
            # if line_i % 1000 == 0:
            #     print(line_i)
            l = list(map(str.strip, line.split()))
            chrom = int(l[0])
            if chrom != cur_chrom:
                chrom_pos_snp_dict[chrom] = {'positions': [], 'snps': []}
                cur_chrom = chrom
            chrom_pos_snp_dict[chrom]['positions'].append(int(l[3]))
            snp = np.zeros(num_individs, dtype='int8')
            j = 0
            w_missing = False
            for i in range(4, 2 * num_individs + 4, 2):
                nt1 = int(l[i])
                nt2 = int(l[i + 1])
                if nt1 == 0 or nt2 == 0:
                    snp[j] = 3
                    w_missing = True
                elif nt1 == 2 and nt2 == 2:
                    snp[j] = 2
                elif nt1 != 1 or nt2 != 1:
                    snp[j] = 1
                #                    #Calculating K
                #                    for ind_i in range(j):
                #                        if snp[j] != 3 and snp[ind_i] != 3:
                #                            k_mat[ind_i, j] = int(snp[j] == snp[ind_i]) + 0.5 * int(sp.absolute(snp[j] - snp[ind_i]) == 1)
                #                            k_mat[ind_i, j] += 1
                j += 1
            #                print k_mat

            bin_counts = np.bincount(snp)
            if w_missing:

                if imputation_type == 'simple': # assign mean value to missing value
                    mean = (bin_counts[1] + 2 * bin_counts[2]) / (bin_counts[0] + bin_counts[1] + bin_counts[2])
                    snp[snp == 3] = round(mean)
                if imputation_type == 'simple2':
                    snp[snp == 3] = np.argmax(bin_counts[:-1])

            chrom_pos_snp_dict[chrom]['snps'].append(snp)

    chromosomes = sorted(chrom_pos_snp_dict.keys())
    snpsds = []
    for chrom in chromosomes:
        snps = chrom_pos_snp_dict[chrom]['snps']
        positions = chrom_pos_snp_dict[chrom]['positions']
        snpsds.append(SNPsData(snps, positions, accessions=individs, chromosome=chrom))
    sd = SNPsDataSet(snpsds, chromosomes, data_format='diploid_int')

    return sd



import sys
import os

#read genotype data
snp = parse_plink_tped_file("./snp", imputation_type='simple') # simple means assign mean value to missing value
indel = parse_plink_tped_file("./indel3", imputation_type='simple')
orf = parse_plink_tped_file("./orfn", imputation_type='simple')

#read phenotype data
f = open("./snp.tfam", 'rU')
#d = {'name':sys.argv[4], 'ecotypes':[], 'values':[]}
d = {'name':"LD", 'ecotypes':[], 'values':[]}
for line in f:
    l = line.split()
    d['ecotypes'].append(l[0])
    d['values'].append(float(l[5]))
pid = 1
phen_dict = {}
phen_dict[pid] = d
f.close()

phend = phenotype_data(phen_dict=phen_dict, phen_ids=phen_dict.keys())
# get the genotypic variants id to line number map. line number could be queryed from the genotypic variants
snp_dic = {}
f1 = open("./snp.tped", 'rU')
i = 0
for line in f1:
    l = line.split()
    snp_dic[l[1]] = i # snp_id snp_line_number
    i = i+1
f1.close()

indel_dic = {}
f2 = open("./indel3.tped", 'rU')
i = 0
for line in f2:
    l = line.split()
    indel_dic[l[1]] = i
    i = i+1
f2.close()

orf_dic = {}
f3 = open("./orfn.tped", 'rU')
i = 0
for line in f3:
    l = line.split()
    orf_dic[l[1]] = i
    i = i + 1
f3.close()

# get a list of significant genotypic variants and transform genotypic variant id to line number
sig_snp_dic = {}
snp_indexes = []
f4 = open("./snp_summary.txt", 'rU')
for line in f4:
    l = line.split()
    if len(l) > 1:
        if l[1] in snp_dic:
            sig_snp_dic[l[1]] = 1
            snp_indexes.append(snp_dic[l[1]]) # snp line_number
        else:
            print (l[1])
f4.close()

indel_indexes = []
f5 = open("./indel_summary.txt", 'rU')
for line in f5:
    l = line.split()
    if len(l) > 1:
        if l[1] in indel_dic:
            indel_indexes.append(indel_dic[l[1]])
        else:
            print(l[1])
f5.close()

orf_indexes = []
f6 = open("./orf_summary.txt", 'rU')
for line in f6:
    l = line.split()
    if len(l) > 1:
        if l[1] in orf_dic:
            orf_indexes.append(orf_dic[l[1]])
        else:
            print(l[1])
f6.close()

# remove significant snps from the complete list
# reconstruct kinship matrix
# and read kinship matrix
f = open("significant_snps",'w')
for snp_id in sig_snp_dic:
    f.write(snp_id + "\t" + snp_id+"\n")
f.close()
os.system('plink -tfile snp --exclude significant_snps --recode transpose --out snp_significant')
os.system('emmax-kin-intel64 -v -d 10 snp_significant')
k = np.loadtxt("snp_significant.aBN.kinf")
k = np.mat(k)

# get the genotypic variants value of all the significant variants
xsSnp = []
xsIndel = []
xsOrf = []

xsSnpIndel = []
xsSnpOrf = []
xsIndelOrf = []

xsSnpIndelOrf = []

for snp_index in snp_indexes:
    xsSnp.append(snp.get_snps()[snp_index])
    xsSnpIndel.append(snp.get_snps()[snp_index])
    xsSnpIndelOrf.append(snp.get_snps()[snp_index])
    xsSnpOrf.append(snp.get_snps()[snp_index])

for indel_index in indel_indexes:
    xsIndel.append(indel.get_snps()[indel_index])
    xsSnpIndel.append(indel.get_snps()[indel_index])
    xsSnpIndelOrf.append(indel.get_snps()[indel_index])
    xsIndelOrf.append(indel.get_snps()[indel_index])

for orf_index in orf_indexes:
    xsOrf.append(orf.get_snps()[orf_index])
    xsSnpIndelOrf.append(orf.get_snps()[orf_index])
    xsSnpOrf.append(orf.get_snps()[orf_index])
    xsIndelOrf.append(orf.get_snps()[orf_index])

lmm = LinearMixedModel(phend.get_values(1))
lmm.add_random_effect(k)


orig_stdout = sys.stdout
f = open('mixmodel_emma.var_prec', 'w')
sys.stdout = f

var_perc = lmm.variance_explained_by_significant_loci( xsSnp, K=k, includedSnpThreadshold=1e-4, lin_depend_thres=1e-4)
print ("SNP: ", var_perc)
lmm.clear_factors()

var_perc = lmm.variance_explained_by_significant_loci( xsIndel, K=k, includedSnpThreadshold=1e-4, lin_depend_thres=1e-4)
print ("INDEL: ", var_perc)
lmm.clear_factors()

var_perc = lmm.variance_explained_by_significant_loci( xsOrf, K=k, includedSnpThreadshold=1e-4, lin_depend_thres=1e-4)
print ("ORF: ", var_perc)
lmm.clear_factors()

var_perc = lmm.variance_explained_by_significant_loci( xsSnpIndel, K=k, includedSnpThreadshold=1e-4, lin_depend_thres=1e-4)
print ("SNPINDEL: ", var_perc)
lmm.clear_factors()

var_perc = lmm.variance_explained_by_significant_loci( xsSnpOrf, K=k, includedSnpThreadshold=1e-4, lin_depend_thres=1e-4)
print ("SNPORF: ", var_perc)
lmm.clear_factors()

var_perc = lmm.variance_explained_by_significant_loci( xsIndelOrf, K=k, includedSnpThreadshold=1e-4, lin_depend_thres=1e-4)
print ("INDELORF: ", var_perc)
lmm.clear_factors()

var_perc = lmm.variance_explained_by_significant_loci( xsSnpIndelOrf, K=k, includedSnpThreadshold=1e-4, lin_depend_thres=1e-4)
print ("SNPINDELORF: ", var_perc)
lmm.clear_factors()

sys.stdout = orig_stdout
f.close()
