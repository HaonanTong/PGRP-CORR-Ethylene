import pandas as pd
import numpy as np;
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# def f_interp( data ):
# 	t = [0, 0.25, 0.5, 1, 4, 12, 24]
# 	f = interp1d(t,data)
# 	x_new = np.arange(0, 24.25, 0.25)
# 	y_interp = f(x_new)
# 	return pd.Series(y_interp)

# db = pd.read_csv('DB6.csv',index_col = 'id')

# # ---------------------------------------------------
# # Generate DB contains interpolation information
# # ---------------------------------------------------
# data = db.loc[:,'t0':'t6']
# data = data.T
# # - Interpolation
# newDF = pd.DataFrame()
# for i in range(0,len(data.columns)):
# 	sample = data.iloc[:,i]
# 	sample_interp = f_interp(sample)
# 	newDF[sample.name] = sample_interp

# newDF = newDF.T
# newDF.rename(columns=lambda x: 'it'+str(x), inplace=True)
# result = pd.concat([db, newDF], axis=1, join='inner')

# # print result.head()
# result.to_csv('DB7.csv')
# # raw_input()
# data = newDF


# # ---------------------------------------------------
# # Plot all interpolated trajectories at a time
# # ---------------------------------------------------
# db = pd.read_csv('DB7.csv',index_col = 'id')
# db = db[db['ethylene']==1]
# # Plot All trajectories in ATP 
# for atp in range(2,6):
# 	data = db.loc[(db['ATP']==atp) | ((db['ATP']<atp) & db['TF']==1),'t0':'t6']
# 	# print len(data)
# 	# data = ( data - data.mean() ) / data.std()
# 	# data = data.cumsum()
# 	# plt.figure()
# 	data = data.T
# 	# - Interpolation
# 	newDF = pd.DataFrame()
# 	for i in range(0,len(data.columns)):
# 		sample = data.iloc[:,i]
# 		sample_interp = f_interp(sample)
# 		newDF[sample.name] = sample_interp

# 	newDF = ( newDF - newDF.mean() ) / newDF.std()
# 	data = newDF

# 	gh = data.plot(x = np.arange(0, 24.25, 0.25), title='Trajectories for Genes And Corresponding Potential Regulators, \n ATP = ' 
# 		+ str(atp), lw=2)
# 	gh.set_xlabel("Time (hour)")
# 	gh.set_ylabel("Normalized Expression Level")	
# 	gh.xaxis.set_major_locator(plt.FixedLocator([0,0.5,1,4,12,24]))
# 	gh.xaxis.set_minor_locator(plt.FixedLocator([0.25]))
# 	gh.legend(loc = 'best')
# 	gh.grid('on', which='minor', axis='x' )
# 	gh.grid('on', which='major', axis='x' )

# 	fig = plt.gcf()
# 	fig.set_size_inches(18.5,10.5)
# 	fig.savefig('Img/Trajectories_RealScale_atp_'+str(atp)+'.pdf',bbox_inches='tight')


# # ------------------------------------------
# # Plot All trajectories in ATP But Highleights Targets 
# # ------------------------------------------
# db = pd.read_csv('DB7.csv',index_col = 'id')
# for atp in range(2,6):
# 	data_TF = db.loc[( ( db['ethylene']==1 ) & ( db['ATP']==atp ) & ( db['TF']==1 ) )
# 	 | ( ( db['ethylene']==1 ) & (db['ATP']<atp) & (db['TF']==1) ),'it0':'it96']
# 	# print len(data_TF)
# 	data_nTF = db.loc[( ( db['ethylene']==1 ) & (db['ATP']==atp) & (db['TF']==0) ),'it0':'it96']

# 	data_TF = data_TF.T
# 	data_nTF = data_nTF.T

# 	data_TF = ( data_TF - data_TF.mean() ) / data_TF.std()
# 	data_nTF = (data_nTF - data_nTF.mean() )/data_nTF.std()

# 	# Plot Data
# 	gh = data_TF.plot(x = np.arange(0, 24.25, 0.25),
# 		lw=2,
# 		style = 'lightgrey')
# 	data_nTF.plot(x = np.arange(0, 24.25, 0.25),
# 		title='Trajectories for Targets And Corresponding Potential Regulators, \n ATP = ' 
# 		+ str(atp), 
# 		lw=2, style = 'blue', ax=gh)

# 	gh.set_xlabel("Time (hour)")
# 	gh.set_ylabel("Normalized Expression Level")	
# 	gh.xaxis.set_major_locator(plt.FixedLocator([0,0.5,1,4,12,24]))
# 	gh.xaxis.set_minor_locator(plt.FixedLocator([0.25]))
# 	gh.legend(loc = 'best')
# 	gh.grid('on', which='minor', axis='x' )
# 	gh.grid('on', which='major', axis='x' )

# 	# plt.legend((gh,gh2),('Potential Regulators', 'Targets Of Interest'), loc = 'best')

# 	# Output Figures
# 	fig = plt.gcf()
# 	fig.set_size_inches(18.5,10.5)
# 	fig.savefig('Img/Trajectories_Highlighted_RealScale__atp_'+str(atp)+'.pdf',bbox_inches='tight')
	
#########################################################
#
#
#
#########################################################
# -------------------------------------------------------
# Focus on AT1G04310
# Img/Trajectories_RealScale__AT1G04310.pdf
# -------------------------------------------------------
db = pd.read_csv('DB7.csv',index_col = 'id')
atp = 3
data_TF = db.loc[( ( db['ethylene']==1 ) & ( db['ATP']==atp ) & ( db['TF']==1 ) ) 
	| ( ( db['ethylene']==1 ) & (db['ATP']<atp) & (db['TF']==1) ),'it0':'it96']
# print len(data_TF)
data_nTF = db.loc[( ( db['ethylene']==1 ) & (db['ATP']==atp) & (db['TF']==0) ),'it0':'it96']

data_TF = data_TF.T
data_nTF = data_nTF.T

data_TF = ( data_TF - data_TF.mean() ) / data_TF.std()
data_nTF = (data_nTF - data_nTF.mean() )/data_nTF.std()

# print data_nTF

xt = np.arange(0, 24.25, 0.25)
# # Plot
gh = data_TF.plot(x = np.arange(0, 24.25, 0.25),
	lw=2 , style = 'lightgrey')


data_nTF.AT1G04310.to_frame().plot(x = np.arange(0, 24.25, 0.25),
	title='Trajectories for Targets And Corresponding Potential Regulators, \n ATP = ' 
	+ str(atp), 
	lw=2, style = 'blue', ax=gh)


# gh = data_nTF..plot(x = np.arange(0, 24.25, 0.25),	
# 	title='Trajectories for AT1G04310 And Corresponding Potential Regulators',
# 	lw=4)

gh.set_xlabel("Time (hour)")
gh.set_ylabel("Normalized Expression Level")	
gh.xaxis.set_major_locator(plt.FixedLocator([0,0.5,1,4,12,24]))
gh.xaxis.set_minor_locator(plt.FixedLocator([0.25]))
# gh.xaxis.set_major_locator(plt.FixedLocator([0,2,4,16,48,96]))
# gh.xaxis.set_minor_locator(plt.FixedLocator([1]))
# gh.legend(loc = 'best')
gh.grid('on', which='minor', axis='x' )
gh.grid('on', which='major', axis='x' )

# Output Figures
fig = plt.gcf()
fig.set_size_inches(18.5,10.5)
fig.savefig('RealScale/Trajectories_RealScale__AT1G04310.pdf',bbox_inches='tight')


# raw_input('Press to continue..')

#########################################################
# Interpolation and Pearson Correlation
#########################################################
# -------------------------------------------------------
# Focus on AT1G04310 and AT5G13910, Interp and Pearson Correlation
# -------------------------------------------------------

gh = data_TF.plot(x = np.arange(0, 24.25, 0.25), 
	title='Trajectories for AT1G04310 And Corresponding Potential Regulators'
	, lw=2 , legend=False , style = 'lightgrey')
data_nTF.AT1G04310.to_frame().plot(x = np.arange(0, 24.25, 0.25), lw=4,ax=gh, legend = True)
data_TF.AT5G13910.to_frame().plot(x = np.arange(0, 24.25, 0.25), lw=4,ax=gh, legend = True)

gh.set_xlabel("Time (hour)")
gh.set_ylabel("Normalized Expression Level")	
gh.xaxis.set_major_locator(plt.FixedLocator([0,0.5,1,4,12,24]))
gh.xaxis.set_minor_locator(plt.FixedLocator([0.25]))
gh.legend(loc = 'best')
gh.grid('on', which='minor', axis='x' )
gh.grid('on', which='major', axis='x' )

# gh.legend(loc = 'best')

# Output Figures
fig = plt.gcf()
fig.set_size_inches(18.5,10.5)
fig.savefig('RealScale/Trajectories_RealScale__AT1G04310_Interp_PC.pdf',bbox_inches='tight')


# -------------------------------------------------------
# Focus on AT1G04310 and AT5G13910, Interp and 4 highest Pearson Correlation
gh = data_TF.plot( x = np.arange(0, 24.25, 0.25), title='Trajectories for AT1G04310 And Corresponding Potential Regulators'
	, lw=2 , legend=False , style = 'lightgrey')
data_nTF.AT1G04310.to_frame().plot(x = np.arange(0, 24.25, 0.25), lw=4,ax=gh, legend = True)
data_TF.AT5G13910.to_frame().plot( x = np.arange(0, 24.25, 0.25), lw=4,ax=gh, legend = True)
data_TF.AT4G16750.to_frame().plot( x = np.arange(0, 24.25, 0.25), lw=4,ax=gh, legend = True)
data_TF.AT1G25560.to_frame().plot( x = np.arange(0, 24.25, 0.25), lw=4,ax=gh, legend = True)
data_TF.AT5G47230.to_frame().plot( x = np.arange(0, 24.25, 0.25), lw=4,ax=gh, legend = True)

gh.set_xlabel("Time (hour)")
gh.set_ylabel("Normalized Expression Level")	
gh.xaxis.set_major_locator(plt.FixedLocator([0,0.5,1,4,12,24]))
gh.xaxis.set_minor_locator(plt.FixedLocator([0.25]))
gh.legend(loc = 'best')
gh.grid('on', which='minor', axis='x' )
gh.grid('on', which='major', axis='x' )

# gh.legend(loc = 'best')

# Output Figures
fig = plt.gcf()
fig.set_size_inches(18.5,10.5)
fig.savefig('RealScale/Trajectories_RealScale_AT1G04310_Interp_PC_Highest4.pdf',bbox_inches='tight')



#########################################################
# Non-Interpolation and Pearson Correlation
#########################################################
# -------------------------------------------------------
# Focus on AT1G04310 and AT3G23240, and 7 time point Pearson Correlation
gh = data_TF.plot(x = np.arange(0, 24.25, 0.25), title='Trajectories for AT1G04310 And Corresponding Potential Regulators'
	, lw=2 , legend=False , style = 'lightgrey')
data_nTF.AT1G04310.to_frame().plot(x = np.arange(0, 24.25, 0.25), lw=4,ax=gh, legend = True)
data_TF.AT3G23240.to_frame().plot(x = np.arange(0, 24.25, 0.25), lw=4,ax=gh, legend = True)

gh.set_xlabel("Time (hour)")
gh.set_ylabel("Normalized Expression Level")	
gh.xaxis.set_major_locator(plt.FixedLocator([0,0.5,1,4,12,24]))
gh.xaxis.set_minor_locator(plt.FixedLocator([0.25]))
gh.legend(loc = 'best')
gh.grid('on', which='minor', axis='x' )
gh.grid('on', which='major', axis='x' )

# gh.legend(loc = 'best')

# Output Figures
fig = plt.gcf()
fig.set_size_inches(18.5,10.5)
fig.savefig('RealScale/Trajectories_RealScale_AT1G04310_PC.pdf',bbox_inches='tight')

# -------------------------------------------------------
# Focus on AT1G04310, and 7 time point 4 highest Pearson Correlation
gh = data_TF.plot(x = np.arange(0, 24.25, 0.25), title='Trajectories for AT1G04310 And Corresponding Potential Regulators'
	, lw=2 , legend=False , style = 'lightgrey')
data_nTF.AT1G04310.to_frame().plot(x = np.arange(0, 24.25, 0.25),lw=4,ax=gh, legend = True)
data_TF.AT3G23240.to_frame().plot(x = np.arange(0, 24.25, 0.25), lw=4,ax=gh, legend = True)
data_TF.AT1G28370.to_frame().plot(x = np.arange(0, 24.25, 0.25), lw=4,ax=gh, legend = True)
data_TF.AT5G13910.to_frame().plot(x = np.arange(0, 24.25, 0.25), lw=4,ax=gh, legend = True)
data_TF.AT5G47220.to_frame().plot(x = np.arange(0, 24.25, 0.25), lw=4,ax=gh, legend = True)

gh.set_xlabel("Time (hour)")
gh.set_ylabel("Normalized Expression Level")	
gh.xaxis.set_major_locator(plt.FixedLocator([0,0.5,1,4,12,24]))
gh.xaxis.set_minor_locator(plt.FixedLocator([0.25]))
gh.legend(loc = 'best')
gh.grid('on', which='minor', axis='x' )
gh.grid('on', which='major', axis='x' )

# Output Figures
fig = plt.gcf()
fig.set_size_inches(18.5,10.5)
fig.savefig('RealScale/Trajectories_RealScale_AT1G04310_PC_Highest4.pdf',bbox_inches='tight')



#########################################################
# DBN
#########################################################
# -------------------------------------------------------
# Focus on AT1G04310 DBN, 1 Parent with Highest Score
# -------------------------------------------------------

gh = data_TF.plot(x = np.arange(0, 24.25, 0.25), title='Trajectories for AT1G04310 And Corresponding Potential Regulators'
	, lw=2 , legend=False , style = 'lightgrey')
data_nTF.AT1G04310.to_frame().plot(x = np.arange(0, 24.25, 0.25), lw=4,ax=gh, legend = True)
data_TF.AT1G74840.to_frame().plot(x = np.arange(0, 24.25, 0.25), lw=4,ax=gh, legend = True)

gh.set_xlabel("Time (hour)")
gh.set_ylabel("Normalized Expression Level")	
gh.xaxis.set_major_locator(plt.FixedLocator([0,0.5,1,4,12,24]))
gh.xaxis.set_minor_locator(plt.FixedLocator([0.25]))
gh.legend(loc = 'best')
gh.grid('on', which='minor', axis='x' )
gh.grid('on', which='major', axis='x' )

# Output Figures
fig = plt.gcf()
fig.set_size_inches(18.5,10.5)
fig.savefig('RealScale/Trajectories_RealScale_AT1G04310_DBN_Highest.pdf',bbox_inches='tight')


# -------------------------------------------------------
# Focus on AT1G04310 DBN, 1 Parent with 4 Highest Score
# -------------------------------------------------------

gh = data_TF.plot(x = np.arange(0, 24.25, 0.25), title='Trajectories for AT1G04310 And Corresponding Potential Regulators'
	, lw=2 , legend=False , style = 'lightgrey')
data_nTF.AT1G04310.to_frame().plot(x = np.arange(0, 24.25, 0.25), lw=4,ax=gh, legend = True)
data_TF.AT1G74840.to_frame().plot(x = np.arange(0, 24.25, 0.25), lw=4,ax=gh, legend = True)
data_TF.AT3G61630.to_frame().plot(x = np.arange(0, 24.25, 0.25), lw=4,ax=gh, legend = True)
data_TF.AT1G64380.to_frame().plot(x = np.arange(0, 24.25, 0.25), lw=4,ax=gh, legend = True)
data_TF.AT4G11140.to_frame().plot(x = np.arange(0, 24.25, 0.25), lw=4,ax=gh, legend = True)

gh.set_xlabel("Time (hour)")
gh.set_ylabel("Normalized Expression Level")	
gh.xaxis.set_major_locator(plt.FixedLocator([0,0.5,1,4,12,24]))
gh.xaxis.set_minor_locator(plt.FixedLocator([0.25]))
gh.legend(loc = 'best')
gh.grid('on', which='minor', axis='x' )
gh.grid('on', which='major', axis='x' )

# Output Figures
fig = plt.gcf()
fig.set_size_inches(18.5,10.5)
fig.savefig('RealScale/Trajectories_RealScale_AT1G04310_DBN_Highest4.pdf',bbox_inches='tight')


