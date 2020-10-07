import os
import numpy as np
import scipy.stats
import pandas as pd
import matplotlib.pyplot as plt

#############################################################
code_directory = '/home/joshs/GitHub/neuropixels_platform_paper'
###################################################

dataDir = os.path.join(code_directory, 'data')

regions = ['LGd','V1','LM','RL','LP','AL','PM','AM']

hierScore = [-0.5150279628298357,
             -0.35733209934482374,
             -0.09388855125761343,
             -0.05987132463908328,
             0.10524780962600731,
             0.15221797920142832,
             0.32766807486511995,
             0.440986074378801]

hierColors = np.array([[217,141,194], # LGd
                       [129,116,177], # V1
                       [78,115,174], # LM
                       [101,178,201], # RL
                       [88,167,106], #LP
                       [202,183,120], # AL
                       [219,132,87], # PM
                       [194,79,84]] # AM
                     ).astype(float)
hierColors /= 255


changeModData = pd.read_csv(os.path.join(dataDir,'change_modulation_data.csv'))

decodingData = pd.read_csv(os.path.join(dataDir,'decoding_data.csv'))

decoderAccuracyVsNCells = pd.read_csv(os.path.join(dataDir,'decoder_accuracy_vs_number_of_cells.csv'))


# Figure4e
fig = plt.figure(facecolor='w',figsize=(6,5))
ax = fig.add_subplot(1,1,1)
d = [changeModData['Time To First Spike'][changeModData['Region']==r] for r in regions]
mn = [np.nanmean(regionData) for regionData in d]
ci = [np.percentile([np.nanmean(np.random.choice(regionData,len(regionData),replace=True)) for _ in range(5000)],(2.5,97.5)) for regionData in d]
slope,yint,rval,pval,stderr = scipy.stats.linregress(hierScore,mn)
x = np.array([min(hierScore),max(hierScore)])
ax.plot(x,slope*x+yint,'--',color='k')
r,p = scipy.stats.pearsonr(hierScore,mn)
title = 'Pearson: r = '+str(round(r,2))+', p = '+str(round(p,3))
r,p = scipy.stats.spearmanr(hierScore,mn)
title += '\nSpearman: r = '+str(round(r,2))+', p = '+str(round(p,3))
for h,m,c,clr in zip(hierScore,mn,ci,hierColors):
    ax.plot(h,m,'o',mec=clr,mfc=clr,ms=6)
    ax.plot([h,h],c,color=clr)
for side in ('right','top'):
    ax.spines[side].set_visible(False)
ax.tick_params(direction='out',top=False,right=False,labelsize=10)
ax.set_xticks([-0.4,-0.2,0,0.2,0.4])
ax.set_xlabel('Anatomical hierarchy score',fontsize=12)
ax.set_ylabel('Time to first spike (ms)',fontsize=12)
ax.set_title(title,fontsize=8)
plt.tight_layout()


# Figure 4f
fig = plt.figure(facecolor='w',figsize=(6,6))
ax = fig.add_subplot(1,1,1)
for state,fill,fitClr in zip(('Active','Passive'),(True,False),('k','0.5')):
    d = [changeModData['Change Modulation '+state][changeModData['Region']==r] for r in regions]
    mn = [np.nanmean(regionData) for regionData in d]
    ci = [np.percentile([np.nanmean(np.random.choice(regionData,len(regionData),replace=True)) for _ in range(5000)],(2.5,97.5)) for regionData in d]
    slope,yint,rval,pval,stderr = scipy.stats.linregress(hierScore,mn)
    x = np.array([min(hierScore),max(hierScore)])
    ax.plot(x,slope*x+yint,'--',color=fitClr)
    r,p = scipy.stats.pearsonr(hierScore,mn)
    if state=='Active':
        title = ''
    else:
        title +='\n'
    title += 'Pearson ('+state+'): r = '+str(round(r,2))+', p = '+str(round(p,3))
    r,p = scipy.stats.spearmanr(hierScore,mn)
    title += '\nSpearman ('+state+'): r = '+str(round(r,2))+', p = '+str(round(p,3))
    for i,(h,m,c,clr) in enumerate(zip(hierScore,mn,ci,hierColors)):
        mfc = clr if fill else 'none'
        lbl = state if i==0 else None
        ax.plot(h,m,'o',mec=clr,mfc=mfc,ms=6,label=lbl)
        ax.plot([h,h],c,color=clr)
for side in ('right','top'):
    ax.spines[side].set_visible(False)
ax.tick_params(direction='out',top=False,right=False,labelsize=10)
ax.set_ylim([0,0.36])
ax.set_xticks([-0.4,-0.2,0,0.2,0.4])
ax.set_xlabel('Anatomical hierarchy score',fontsize=12)
ax.set_ylabel('Change Modulation Index',fontsize=12)
ax.set_title(title,fontsize=8)
ax.legend(loc='upper left')
plt.tight_layout()


# Figure 4h
fig = plt.figure(facecolor='w',figsize=(5,5))
ax = fig.add_subplot(1,1,1)
lbl = 'Correleation of decoder prediction and mouse behavior'
mn = []
for i,(r,h,clr) in enumerate(zip(regions,hierScore,hierColors)):
    d = decodingData[lbl][decodingData['Region']==r]
    m = np.mean(d)
    s = np.std(d)/(len(d)**0.5)
    ax.plot(h,m,'o',mec=clr,mfc=clr)
    ax.plot([h,h],[m-s,m+s],color=clr)
    mn.append(m)
slope,yint,rval,pval,stderr = scipy.stats.linregress(hierScore,mn)
x = np.array([min(hierScore),max(hierScore)])
ax.plot(x,slope*x+yint,'--',color='0.5')
r,p = scipy.stats.pearsonr(hierScore,mn)
title = 'Pearson: r = '+str(round(r,2))+', p = '+str(round(p,3))
r,p = scipy.stats.spearmanr(hierScore,mn)
title += '\nSpearman: r = '+str(round(r,2))+', p = '+str(round(p,3))
for side in ('right','top'):
    ax.spines[side].set_visible(False)
ax.tick_params(direction='out',top=False,right=False,labelsize=10)
ax.set_xlabel('Anatomical hierarchy score',fontsize=12)
ax.set_ylabel(lbl,fontsize=12)
ax.set_title(title,fontsize=8)
plt.tight_layout()


# Extended Data Figure 9k
fig = plt.figure(facecolor='w',figsize=(6,6))
ax = fig.add_subplot(1,1,1)
for trials,fill,fitClr in zip(('Hit','Miss'),(True,False),('k','0.5')):
    d = [changeModData['Change Modulation '+trials][changeModData['Region']==r] for r in regions]
    mn = [np.nanmean(regionData) for regionData in d]
    ci = [np.percentile([np.nanmean(np.random.choice(regionData,len(regionData),replace=True)) for _ in range(5000)],(2.5,97.5)) for regionData in d]
    slope,yint,rval,pval,stderr = scipy.stats.linregress(hierScore,mn)
    x = np.array([min(hierScore),max(hierScore)])
    ax.plot(x,slope*x+yint,'--',color=fitClr)
    r,p = scipy.stats.pearsonr(hierScore,mn)
    if trials=='Hit':
        title = ''
    else:
        title +='\n'
    title += 'Pearson ('+trials+'): r = '+str(round(r,2))+', p = '+str(round(p,3))
    r,p = scipy.stats.spearmanr(hierScore,mn)
    title += '\nSpearman ('+trials+'): r = '+str(round(r,2))+', p = '+str(round(p,3))
    for i,(h,m,c,clr) in enumerate(zip(hierScore,mn,ci,hierColors)):
        mfc = clr if fill else 'none'
        lbl = trials if i==0 else None
        ax.plot(h,m,'o',mec=clr,mfc=mfc,ms=6,label=lbl)
        ax.plot([h,h],c,color=clr)
for side in ('right','top'):
    ax.spines[side].set_visible(False)
ax.tick_params(direction='out',top=False,right=False,labelsize=10)
ax.set_ylim([0,0.36])
ax.set_xticks([-0.4,-0.2,0,0.2,0.4])
ax.set_xlabel('Anatomical hierarchy score',fontsize=12)
ax.set_ylabel('Change Modulation Index',fontsize=12)
ax.set_title(title,fontsize=8)
ax.legend(loc='upper left')
plt.tight_layout()


# Extended Data Figure 9l-n
for rate in ('Pre-change Response','Change Response','Baseline Rate'):
    fig = plt.figure(facecolor='w',figsize=(6,5))
    ax = fig.add_subplot(1,1,1)
    ymax = 0
    for state,fill in zip(('Active','Passive'),(True,False)):
        d = [changeModData[rate+' '+state][changeModData['Region']==r] for r in regions]
        mn = [np.nanmean(regionData) for regionData in d]
        ci = [np.percentile([np.nanmean(np.random.choice(regionData,len(regionData),replace=True)) for _ in range(5000)],(2.5,97.5)) for regionData in d]
        r,p = scipy.stats.pearsonr(hierScore,mn)
        if state=='Active':
            title = ''
        else:
            title +='\n'
        title += 'Pearson ('+state+'): r = '+str(round(r,2))+', p = '+str(round(p,3))
        r,p = scipy.stats.spearmanr(hierScore,mn)
        title += '\nSpearman ('+state+'): r = '+str(round(r,2))+', p = '+str(round(p,3))
        for i,(h,m,c,clr) in enumerate(zip(hierScore,mn,ci,hierColors)):
            mfc = clr if fill else 'none'
            lbl = state if i==0 else None
            ax.plot(h,m,'o',mec=clr,mfc=mfc,ms=6,label=lbl)
            ax.plot([h,h],c,color=clr)
            ymax = max(ymax,c[1])
    for side in ('right','top'):
        ax.spines[side].set_visible(False)
    ax.tick_params(direction='out',top=False,right=False,labelsize=10)
    ax.set_ylim([0,1.05*ymax])
    ax.set_xticks([-0.4,-0.2,0,0.2,0.4])
    ax.set_xlabel('Anatomical hierarchy score',fontsize=12)
    ax.set_ylabel(rate+' (spikes/s)',fontsize=12)
    ax.set_title(title,fontsize=8)
    ax.legend(loc='lower right')
    plt.tight_layout()


# Extended Data Figure 9o
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(1,1,1)
d = np.array(decoderAccuracyVsNCells.iloc[:,2:])
mn = np.nanmean(d,axis=0)
sem = np.nanstd(d,axis=0)/(np.sum(~np.isnan(d),axis=0)**0.5)
ncells = np.array(decoderAccuracyVsNCells.columns[2:]).astype(int)
ax.plot(ncells,mn,'ko-')
for x,m,s in zip(ncells,mn,sem):
    ax.plot([x,x],[m-s,m+s],'k')
for side in ('right','top'):
    ax.spines[side].set_visible(False)
ax.tick_params(direction='out',top=False,right=False,labelsize=12)
ax.set_xticks(np.arange(0,100,10))
ax.set_xlim([0,max(ncells)+5])
ax.set_ylim([0.5,1])
ax.set_xlabel('Number of Cells',fontsize=12)
ax.set_ylabel('Decoder Accuracy',fontsize=12)
plt.tight_layout()
    

# Extended Data Figure 9p
fig = plt.figure(facecolor='w',figsize=(5,5))
ax = fig.add_subplot(1,1,1)
lbl = 'Decoder accuracy'
mn = []
for i,(r,h,clr) in enumerate(zip(regions,hierScore,hierColors)):
    d = decodingData[lbl][decodingData['Region']==r]
    m = np.mean(d)
    s = np.std(d)/(len(d)**0.5)
    ax.plot(h,m,'o',mec=clr,mfc=clr)
    ax.plot([h,h],[m-s,m+s],color=clr)
    mn.append(m)
r,p = scipy.stats.pearsonr(hierScore,mn)
title = 'Pearson: r = '+str(round(r,2))+', p = '+str(round(p,3))
r,p = scipy.stats.spearmanr(hierScore,mn)
title += '\nSpearman: r = '+str(round(r,2))+', p = '+str(round(p,3))
for side in ('right','top'):
    ax.spines[side].set_visible(False)
ax.tick_params(direction='out',top=False,right=False,labelsize=10)
ax.set_ylim([0.5,1])
ax.set_xlabel('Anatomical hierarchy score',fontsize=12)
ax.set_ylabel(lbl,fontsize=12)
ax.set_title(title,fontsize=8)
plt.tight_layout()

