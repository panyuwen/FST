## T = -ln(1 - Fst)
## PBS = (T[A&B] + T[A&C] - T[B&C]) / 2.0
## population-specific branch length for A for each SNP

import argparse
import pandas
import numpy as np
from concurrent.futures import ProcessPoolExecutor
import collections

parser = argparse.ArgumentParser()
parser.add_argument("--fst", type=str, help="Variant_FST.txt.gz", required=True)
parser.add_argument("--pop", type=str, required=True, \
					help="<target pop name> <ref1 pop name> <ref2 pop name>, 3 columns, no header.")
parser.add_argument("--out", type=str, help="/path/for/output", required=False, default='./')
args = parser.parse_args()

fstfile = args.fst
popfile = args.pop
outpath = args.out

fst = pandas.read_csv(fstfile,sep='\t',index_col='ID')
pop = pandas.read_csv(popfile,sep='\s+',header=None)

def pbs_one_combination(input):
	target,ref1,ref2 = input
	pbs = pandas.DataFrame(columns=[target+":"+ref1, target+":"+ref2, ref1+":"+ref2], index=list(fst.index))
	for pair in [target+":"+ref1, target+":"+ref2, ref1+":"+ref2]:
		if pair in list(fst.columns):
			pbs[pair] = list(fst[pair])
		elif ':'.join(list(reversed(pair.split(':')))) in list(fst.columns):
			pbs[pair] = list(fst[':'.join(list(reversed(pair.split(':'))))])
		else:
			print('No Fst values found for combination of '+pair+'.')
			exit()
	pbs['t1'] = -1.0 * np.log(1.0 - pbs[target+":"+ref1])
	pbs['t2'] = -1.0 * np.log(1.0 - pbs[target+":"+ref2])
	pbs['t3'] = -1.0 * np.log(1.0 - pbs[ref1+":"+ref2])
	pbs['PBS_'+target+'_'+ref1+'_'+ref2] = (pbs['t1'] + pbs['t2'] - pbs['t3']) / 2.0

	return 'PBS_'+target+'_'+ref1+'_'+ref2, pbs['PBS_'+target+'_'+ref1+'_'+ref2]

with ProcessPoolExecutor(max_workers = min(20,pop.shape[0])) as pool:
	res = list(pool.map(pbs_one_combination, pop.values.tolist()))

df = collections.defaultdict(list)
for p, r in res:
	df[p] = r
df = pandas.DataFrame(df)
df.to_csv(outpath+'/PBS.txt.gz',sep='\t',compression='gzip')

df = df.replace([np.inf, -np.inf], np.nan).dropna()
df.to_csv(outpath+'/PBS.clean.txt.gz',sep='\t',compression='gzip')
