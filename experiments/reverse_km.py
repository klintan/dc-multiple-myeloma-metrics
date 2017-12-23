from lifelines.datasets import load_gbsg2
from lifelines import KaplanMeierFitter

df = load_gbsg2()
kmf = KaplanMeierFitter()


#reverse event
#df['cens'] = df['cens'].astype(bool)
#df['cens'] = ~df['cens']
#df['cens'] = df['cens'].astype(int)
# https://cemsiis.meduniwien.ac.at/fileadmin/msi_akim/CeMSIIS/KB/volltexte/Schemper_Smith_1996_CCT.pdf
# this one actually recommends just swapping 
# https://www.graphpad.com/support/faq/determining-the-median-followup-time-in-survival-anlaysis/

# what is the median follow up time
# https://www2.le.ac.uk/departments/health-sciences/research/biostats/youngsurv/pdf/MShanyinde.pdf
# The median follow-up is an indicator of how ‘mature’ your survival data is (e.g. how many months on ‘average’ the patients were followed since randomisation into the study).

# http://myweb.uiowa.edu/pbreheny/7210/f15/notes/11-5.pdf
kmf.fit(df['time'], event_observed=df['cens'], reverse=True)

#print(kmf.survival_function_)
print(kmf.median_)
