import pandas as pd
import numpy as np;
import matplotlib.pyplot as plt

db = pd.read_csv('DB6.csv',index_col = 'id')
db = db[db['ethylene']==1]

# Plot All trajectories in ATP 
for atp in range(2,6):
	data = db.loc[(db['ATP']==atp) | ((db['ATP']<atp) & db['TF']==1),'t0':'t6']
	# print len(data)
	# data = ( data - data.mean() ) / data.std()
	# data = data.cumsum()
	# plt.figure()
	data = data.T
	data = ( data - data.mean() ) / data.std()
	# print data
	gh = data.plot(title='Trajectories for Genes And Corresponding Potential Regulators, \n ATP = ' 
		+ str(atp), lw=2)
	gh.set_xlabel("Time (hour)")
	gh.set_ylabel("Normalized Expression Level")
	gh.set_xticklabels(['0','0.25','0.5','1','4','12','24'])
	gh.legend(loc = 'best')

	# fig = plt.gcf()
	# fig.set_size_inches(18.5,10.5)
	# fig.savefig('Img/Trajectories_atp_'+str(atp)+'.pdf',bbox_inches='tight')


# Plot All trajectories in ATP But Highleights Targets 
for atp in range(2,6):
	data_TF = db.loc[( ( db['ATP']==atp ) & ( db['TF']==1 ) )
	 | ( (db['ATP']<atp) & (db['TF']==1) ),'t0':'t6']
	# print len(data_TF)
	data_nTF = db.loc[( (db['ATP']==atp) & (db['TF']==0) ),'t0':'t6']

	# print len(data_nTF)

	# print len(data)

	# Normalize Data Set
	data_TF = data_TF.T
	data_TF = ( data_TF - data_TF.mean() ) / data_TF.std()
	data_nTF = data_nTF.T
	data_nTF = ( data_nTF - data_nTF.mean() ) / data_nTF.std()

	# Plot Data
	gh = data_nTF.plot(title='Trajectories for Genes And Corresponding Potential Regulators, \n ATP = ' 
		+ str(atp), lw=4)
	data_TF.plot(lw=2,ax=gh)

	gh.set_xlabel("Time (hour)")
	gh.set_ylabel("Normalized Expression Level")
	gh.set_xticklabels(['0','0.25','0.5','1','4','12','24'])
	gh.legend(loc = 'best')

	# Output Figures
	# fig = plt.gcf()
	# fig.set_size_inches(18.5,10.5)
	# fig.savefig('Img/Trajectories_Highlighted_atp_'+str(atp)+'.pdf',bbox_inches='tight')


# Focus on AT1G04310
atp = 3
data_TF = db.loc[( ( db['ATP']==atp ) & ( db['TF']==1 ) )
	 | ( (db['ATP']<atp) & (db['TF']==1) ),'t0':'t6']

data_nTF = db.loc[( (db['ATP']==atp) & (db['TF']==0) ),'t0':'t6']

# Normalize Data Set
data_TF = data_TF.T
data_TF = ( data_TF - data_TF.mean() ) / data_TF.std()
data_nTF = data_nTF.T
data_nTF = ( data_nTF - data_nTF.mean() ) / data_nTF.std()

# Plot
gh = data_TF.plot(title='Trajectories for AT1G04310 And Corresponding Potential Regulators'
	, lw=2 , legend=False , style = 'lightgrey')
data_nTF.AT1G04310.plot(lw=2,ax=gh, legend = True)

gh.set_xlabel("Time (hour)")
gh.set_ylabel("Normalized Expression Level")
gh.set_xticklabels(['0','0.25','0.5','1','4','12','24'])
# gh.legend(loc = 'best')

# Output Figures
# fig = plt.gcf()
# fig.set_size_inches(18.5,10.5)
# fig.savefig('Img/Trajectories_AT1G04310.pdf',bbox_inches='tight')


# Focus on AT1G04310 and AT5G13910, Interp and Pearson Correlation
gh = data_TF.plot(title='Trajectories for AT1G04310 And Corresponding Potential Regulators'
	, lw=2 , legend=False , style = 'lightgrey')
data_nTF.AT1G04310.plot(lw=4,ax=gh, legend = True)
data_TF.AT5G13910.plot(lw=4,ax=gh, legend = True)

gh.set_xlabel("Time (hour)")
gh.set_ylabel("Normalized Expression Level")
gh.set_xticklabels(['0','0.25','0.5','1','4','12','24'])
# gh.legend(loc = 'best')

# Output Figures
fig = plt.gcf()
fig.set_size_inches(18.5,10.5)
fig.savefig('Img/Trajectories_AT1G04310_Interp_PC.pdf',bbox_inches='tight')


# Focus on AT1G04310 and AT5G13910, Interp and 4 highest Pearson Correlation
gh = data_TF.plot(title='Trajectories for AT1G04310 And Corresponding Potential Regulators'
	, lw=2 , legend=False , style = 'lightgrey')
data_nTF.AT1G04310.plot(lw=4,ax=gh, legend = True)
data_TF.AT5G13910.plot(lw=4,ax=gh, legend = True)
data_TF.AT4G16750.plot(lw=4,ax=gh, legend = True)
data_TF.AT1G25560.plot(lw=4,ax=gh, legend = True)
data_TF.AT5G47230.plot(lw=4,ax=gh, legend = True)

gh.set_xlabel("Time (hour)")
gh.set_ylabel("Normalized Expression Level")
gh.set_xticklabels(['0','0.25','0.5','1','4','12','24'])
# gh.legend(loc = 'best')

# Output Figures
fig = plt.gcf()
fig.set_size_inches(18.5,10.5)
fig.savefig('Img/Trajectories_AT1G04310_Interp_PC_Highest4.pdf',bbox_inches='tight')



# Focus on AT1G04310 and AT3G23240, and 7 time point Pearson Correlation
gh = data_TF.plot(title='Trajectories for AT1G04310 And Corresponding Potential Regulators'
	, lw=2 , legend=False , style = 'lightgrey')
data_nTF.AT1G04310.plot(lw=4,ax=gh, legend = True)
data_TF.AT3G23240.plot(lw=4,ax=gh, legend = True)

gh.set_xlabel("Time (hour)")
gh.set_ylabel("Normalized Expression Level")
gh.set_xticklabels(['0','0.25','0.5','1','4','12','24'])
# gh.legend(loc = 'best')

# Output Figures
fig = plt.gcf()
fig.set_size_inches(18.5,10.5)
fig.savefig('Img/Trajectories_AT1G04310_PC.pdf',bbox_inches='tight')

# Focus on AT1G04310, and 7 time point 4 highest Pearson Correlation
gh = data_TF.plot(title='Trajectories for AT1G04310 And Corresponding Potential Regulators'
	, lw=2 , legend=False , style = 'lightgrey')
data_nTF.AT1G04310.plot(lw=4,ax=gh, legend = True)
data_TF.AT3G23240.plot(lw=4,ax=gh, legend = True)
data_TF.AT1G28370.plot(lw=4,ax=gh, legend = True)
data_TF.AT5G13910.plot(lw=4,ax=gh, legend = True)
data_TF.AT5G47220.plot(lw=4,ax=gh, legend = True)

gh.set_xlabel("Time (hour)")
gh.set_ylabel("Normalized Expression Level")
gh.set_xticklabels(['0','0.25','0.5','1','4','12','24'])
# gh.legend(loc = 'best')

# Output Figures
fig = plt.gcf()
fig.set_size_inches(18.5,10.5)
fig.savefig('Img/Trajectories_AT1G04310_PC_Highest4.pdf',bbox_inches='tight')


# Focus on AT1G04310 DBN, 1 Parent with Highest Score
gh = data_TF.plot(title='Trajectories for AT1G04310 And Corresponding Potential Regulators'
	, lw=2 , legend=False , style = 'lightgrey')
data_nTF.AT1G04310.plot(lw=4,ax=gh, legend = True)
data_TF.AT1G74840.plot(lw=4,ax=gh, legend = True)
data_TF.AT3G61630.plot(lw=4,ax=gh, legend = True)
data_TF.AT1G64380.plot(lw=4,ax=gh, legend = True)
data_TF.AT4G11140.plot(lw=4,ax=gh, legend = True)

gh.set_xlabel("Time (hour)")
gh.set_ylabel("Normalized Expression Level")
gh.set_xticklabels(['0','0.25','0.5','1','4','12','24'])
# gh.legend(loc = 'best')

# Output Figures
fig = plt.gcf()
fig.set_size_inches(18.5,10.5)
fig.savefig('Img/Trajectories_AT1G04310_DBN_Highest4.pdf',bbox_inches='tight')


