# -------------------------------------------
# Haonan Tong
# -------------------------------------------
import pandas as pd
import numpy as np
db = pd.read_csv('DB6.csv',index_col = 'id')

# Expression for TFs 
for atp in range(1,7):
	db_tmp = db.loc[(db['ethylene']==1) & (db['ATP']==atp) & (db['TF']==1), 't0':'t6']
	db_tmp = np.log2(db_tmp)
	db_tmp.to_csv('data/log2Expr_C2H4_TF_AT'+str(atp)+'.csv')

# Expression for nTFs 
for atp in range(1,7):
	db_tmp = db.loc[(db['ethylene']==1) & (db['ATP']==atp) & (db['TF']==0), 't0':'t6']
	db_tmp = np.log2(db_tmp)
 	db_tmp.to_csv('data/log2Expr_C2H4_nTF_AT'+str(atp)+'.csv')