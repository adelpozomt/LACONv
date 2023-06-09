#!/usr/bin/python

import sys, re, shlex , os, string, urllib, time, math, random, subprocess, shutil

import ConfigParser

import optparse

from pybedtools import BedTool

import pandas

import scipy.stats as stats

import numpy as np

from itertools import product

######################################################################

class OptionParser(optparse.OptionParser):

    def check_required (self, opt):

        option = self.get_option(opt)

        atrib = getattr(self.values, option.dest)
        
        if atrib is None:
#            self.error("%s option not supplied" % option)
            return False
        else:
            return True

#######################################################################

def predict_CNV(df_all,sample):
    
    df_sample = df_all.loc[:,sample]
    
    df_remaining = df_all.loc[:,df_all.columns != sample]
    
    for interval in list(df_sample.index):
        
        dosis_sample = df_sample.loc[interval]
        
        rates_remaining = list(df_remaining.loc[interval,:].values)
        
        zscore = stats.zmap([dosis_sample],rates_remaining)
        
        if dosis_sample <= 0.7 and zscore <= -2:
            return "CNV_loss"
        elif dosis_sample >= 1.2 and zscore >= 2:
            return "CNV_gain"
        else:
            return "None" 
    

def calculate_dosis_distribution(df_select):
    
    l_dosis = []
    
    for interval in list(df_select.index):
        
        dosis = df_select.loc[interval]
        
        if type(dosis) == np.float64:
            l_dosis.append(dosis)
        else:
            dosis = list(df_select.loc[interval].values)
            l_dosis.extend(dosis)
            
        
    return l_dosis


def calculate_zscore_distribution(df_select,sample):
    
    df_select_sample = df_select.loc[:,sample]
    
    df_select_remaining = df_select.loc[:,df_select.columns != sample]
    
    l_zscore   = []
    
    for interval in list(df_select_sample.index):
        
        rates_sample = [df_select_sample.loc[interval]]
        rates_remaining = list(df_select_remaining.loc[interval,:].values)
        
        zscore = stats.zmap(rates_sample,rates_remaining)
        
        l_zscore.append(zscore)
    
    return l_zscore

#######################################################################

def build_dataframe_allratios(l_samples,hash_query_beds):
    
    inicio = True
    
    df_all = None
    
    hash_index = {}
    
    for sample in l_samples:
        
        bed_sample = hash_query_beds[sample]
        
        if inicio:
            
            l_index = map(lambda intv: "%s:%s-%s" % (intv[0],intv[1],intv[2]), BedTool(bed_sample))
            
            hash_index = dict(map(lambda i: (i,True), l_index))
            
            df_all = pandas.DataFrame(columns=l_samples,index=l_index)
            
            inicio = False
            
        df_all[sample] = map(lambda intv: float(intv[3]), BedTool(bed_sample))
        
    return df_all,hash_index

#######################################################################

def build_dataframe_zscore(df_dosis_all):
    
    l_samples = df_dosis_all.columns
    
    l_intervals = df_dosis_all.index
    
    df_zscore_all = pandas.DataFrame(columns=l_samples,index=l_intervals)
    
    for sample in l_samples: 
        
        l_zscore = []
        
        for interval in l_intervals:
                
            rates_sample = [df_dosis_all.loc[interval,sample]]
            rates_remaining = list(df_dosis_all.loc[interval,df_dosis_all.columns != sample].values)
                   
            zscore = float(stats.zmap(rates_sample,rates_remaining))
            
            l_zscore.append(zscore)
            
        df_zscore_all[sample] = l_zscore 
        
    return df_zscore_all

#######################################################################

def run(argv=None):
    
    if argv is None: argv = sys.argv    
   
    parser = OptionParser(add_help_option=True,description="The script performs CNV estimation within the regions of interest")
    
    parser.add_option("--query_list",default=None,help="File with a list of beds to be tested. Each bed contains the intervals determined by the CNV tool. The file is tab-separated. First column must contain sample name and second includes the path of corresponding bed file",dest="f_query")
    parser.add_option("--target_list",default=False,help="File with a list of beds that contains the intervals detected by mlpa. The file is tab-separated. First column must contain sample name and second includes the path of corresponding bed file",dest="f_target")
    parser.add_option("--s",default=False,help="File with the list of samples",dest="f_samples")
    parser.add_option("--o",default=False,help="Path for output files",dest="output_path")
    parser.add_option("--t",default=False,help="Type: l(loss) g(gain)",dest="cnv_type")
    
    # Se leen las opciones aportadas por el usuario
    (options, args) = parser.parse_args(argv[1:])

    if len(argv) == 1:
        sys.exit(0)
    
    if not parser.check_required("--query_list"):
        raise IOError('adjust_threshold: The file with query intervals has not been provided')
    
    if not parser.check_required("--target_list"):
        raise IOError('adjust_threshold: The file with target intervals has not been provided')
    
    if not parser.check_required("--s"):
        raise IOError('adjust_threshold: The file with the list of input samples has not been provided')
        
    if not parser.check_required("--o"):
        raise IOError('adjust_threshold: The output path has not been provided')
               
    
    output_path = options.output_path
    
    if not os.path.exists(output_path):
        raise IOError("adjust_threshold: The output path does not exist. It must be created before launching the script")
    
    f_query = options.f_query
    
    if not os.path.exists(f_query):
        raise IOError("adjust_threshold: The file with the list of query bed files does not exist. %s" % (f_query))
    
    fi = open(f_query,'r')
    l_query = map(lambda x: x.strip().split('\t'), fi.readlines())
    fi.close()
    
    hash_query_beds = {}
    
    for (sample,f) in l_query:
        
        if not os.path.exists(f):
            raise IOError("adjust_threshold: The bed file does not exist. %s" % (f))
        
        if os.path.getsize(f) > 0:
            hash_query_beds[sample] = f
    
    f_target = options.f_target
    
    if not os.path.exists(f_target):
        raise IOError("adjust_threshold: The file with the list of target bed files does not exist. %s" % (f_target))
    
    fi = open(f_target,'r')
    l_target = map(lambda x: x.strip().split('\t'), fi.readlines())
    fi.close()
    
    hash_target_beds = {}
    
    for (sample,f) in l_target:
        
        if not os.path.exists(f):
            raise IOError("CNV_validacion: The bed file does not exist. %s" % (f))
        
        hash_target_beds[sample] = f
    
    f_samples = options.f_samples
    
    if not os.path.exists(f_samples):
        raise IOError("adjust_threshold: The file with the list of samples does not exist. %s" % (f_samples))
    
    fi = open(f_samples,'r')
    l_samples = map(lambda x: x.strip(), fi.readlines())
    fi.close()
    
    #df_all,hash_intervals = build_dataframe_allratios(l_samples,hash_query_beds)
    
    #df_all.to_csv(os.path.join(output_path,"dosis_all_ratios.tsv"),sep="\t")
    
    #df_zscore_all = build_dataframe_zscore(df_all)        
    
    #df_zscore_all.to_csv(os.path.join(output_path,"dosis_all_zscore.tsv"),sep="\t")
    
    if options.cnv_type == None:
        raise IOError("adjust_threshold: The cnv type has not been provided")
    
    cnv_type = None
    
    l_thresholds_zscore_loss = range(-11,1)
    l_thresholds_zscore_gain = range(1,12)  
    l_threshols_dosis_loss = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]
    l_threshols_dosis_gain = [1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9]
    
    
    if options.cnv_type == "l":
        cnv_type = "CNV_loss"
        l_threshold_dosis  = l_threshols_dosis_loss
        l_threshold_zscore = l_thresholds_zscore_loss
    else:
        cnv_type = "CNV_gain"
        l_threshold_dosis = l_threshols_dosis_gain
        l_threshold_zscore = l_thresholds_zscore_gain
        
    
    df_dosis_all = pandas.read_table(os.path.join(output_path,"dosis_all_ratios.tsv"),sep="\t",header=0,index_col=0)
    df_zscore_all = pandas.read_table(os.path.join(output_path,"dosis_all_zscore.tsv"),sep="\t",header=0,index_col=0)         
    
    l_samples = list(df_dosis_all.columns)
    
    l_thresholds = map(lambda (x,y): (x,y), product(l_threshold_dosis,l_threshold_zscore))
    
    hash_predictions = {}
    
    for sample in l_samples:
        
        mlpa_sample = hash_target_beds.get(sample,None)
        
        if mlpa_sample <> None:
            bed_mlpa = BedTool(mlpa_sample)
            
            if cnv_type == "CNV_loss":
                l_ = filter(lambda intv: intv[6].lower()=="deletion", bed_mlpa)
                if l_ <> []:
                    bed_mlpa_cnv = BedTool(l_)
                else:
                    bed_mlpa_cnv = None
                    
            else:
                l_ = filter(lambda intv: intv[6].lower()=="duplication", bed_mlpa)
                if l_ <> []:
                    bed_mlpa_cnv = BedTool(l_)
                else:
                    bed_mlpa_cnv = None
        else:
            bed_mlpa = None
            bed_mlpa_cnv = None
        
        if bed_mlpa_cnv <> None:
            P = bed_mlpa_cnv.count()
        else:
            P = 0
            
        num_pred = 0
        
        dosis_sample = df_dosis_all.loc[:,sample]         
        
        for (threshold_dosis,threshold_zscore) in l_thresholds:
                
            hash_predictions.setdefault((threshold_dosis,threshold_zscore),{})
        
            hash_predictions[(threshold_dosis,threshold_zscore)][sample] = None
            
            l_intervals_2 = []
        
            if cnv_type == "CNV_loss":
                
                l_intervals = list(dosis_sample.loc[dosis_sample <= threshold_dosis].index)
                
                if l_intervals <> []: 
                    zscores_sample = df_zscore_all.loc[l_intervals,sample]
                    l_intervals_2 = list(zscores_sample.loc[zscores_sample <= threshold_zscore].index)
                    
            else:
                l_intervals = list(dosis_sample.loc[dosis_sample >= threshold_dosis].index)
                if l_intervals <> []: 
                    zscores_sample = df_zscore_all.loc[l_intervals,sample]
                    l_intervals_2 = list(zscores_sample.loc[zscores_sample >= threshold_zscore].index)
                
            num_pred = 0
            TP = 0
            FN = 0
            FP = 0
                
            if l_intervals_2 == [] and bed_mlpa_cnv == None: ## No hay predicciones ni intervalos determinados por mlpa
                continue
            elif l_intervals_2 == [] and bed_mlpa_cnv <> None: ## No hay predicciones y si intervalos determinado por mlps
                FN = P
            elif l_intervals <> []: 
            
                num_pred = len(l_intervals_2)
                    
                if bed_mlpa_cnv == None:
                    FP = len(l_intervals_2) # no hay intervalos determinados por mlpa
                else:
                    if l_intervals_2 <> []: ## hay predicciones
                        #l_intervals_laconv = map(lambda x: re.findall("(chr)*\s*([0-9]{1,2}|X|Y|MT)\s*(-|:)?\s*(\d+)\s*(MB|M|K)?\s*(-|:)?\s*(\d+)\s*(MB|M|K)?",x), l_intervals_2)
                        l_intervals_laconv = map(lambda y: ("chr%s" % y[0][0],y[0][2],y[0][5]), map(lambda x: re.findall("\\S([0-9]{1,2}|X|Y|MT)\\s*(-|:)?\\s*(\\d+)\\s*(MB|M|K)?\\s*(-|:)?\\s*(\\d+)\\s*(MB|M|K)?",x), l_intervals_2))
                        bed_laconv = BedTool(l_intervals_laconv)
                    
                        TP = bed_mlpa_cnv.intersect(bed_laconv).count()
                        FN = bed_mlpa_cnv.intersect(bed_laconv,v=True).count()
                        FP = bed_laconv.intersect(bed_mlpa_cnv,v=True).count()
                    else: ## no hay prediciones
                        FN = P
                
            hash_predictions[(threshold_dosis,threshold_zscore)][sample] = (P,num_pred,TP,FN,FP)
        
    print "Dosis\tzscore\tP\tpredict\tTP\tFN\tFP" 

    for (t_dosis,t_zscore) in sorted(hash_predictions.keys()):
        
        total_P  = 0
        total_pred = 0
        total_TP = 0
        total_FN = 0
        total_FP = 0
        
        for sample in l_samples:
            
            if not hash_predictions[(t_dosis,t_zscore)].has_key(sample):
                continue
            
            prediction = hash_predictions[(t_dosis,t_zscore)][sample]
            
            if prediction <> None:
                
                total_P  += int(prediction[0])
                total_pred += int(prediction[1])
                total_TP += int(prediction[2])
                total_FN += int(prediction[3])
                total_FP += int(prediction[4])
                
                if (prediction[0]+prediction[1]) <> (prediction[2]+prediction[3]+prediction[4]):
                    a = 0
                
        
        print "%1.2f\t%d\t%d\t%d\t%d\t%d\t%d" % (t_dosis,t_zscore,total_P,total_pred,total_TP,total_FN,total_FP)
        
        
    """
    l_zscore_TP = []
    l_dosis_TP  = []
    l_zscore_TN = []
    l_dosis_TN  = []
    
    
    for sample in l_samples:
        
        print sample
        
        if hash_target_beds.has_key(sample): # the sample has True intervals with dosis changes
            
            bed_P      = BedTool(hash_target_beds[sample])
            bed_laconv = BedTool(hash_query_beds[sample])
            
            bed_laconv_TP = bed_laconv.intersect(bed_P)
            bed_laconv_TN = bed_laconv.intersect(bed_P,v=True)
            
            l_index_TP = []
            
            for intv in BedTool(bed_laconv_TP):
                
                key_intv = "%s:%s-%s" % (intv[0],intv[1],intv[2])
                
                if hash_intervals.has_key(key_intv):
                    l_index_TP.append(key_intv)
                else:
                    key_intv = "%s:%d-%s" % (intv[0],int(intv[1])-1,intv[2])
                    if hash_intervals.has_key(key_intv):
                        l_index_TP.append(key_intv)
                    else:
                        print "Problems with key: %s" % (key_intv)
                        
            l_index_TN = []
            
            for intv in BedTool(bed_laconv_TN):
                
                key_intv = "%s:%s-%s" % (intv[0],intv[1],intv[2])
                
                if hash_intervals.has_key(key_intv):
                    l_index_TN.append(key_intv)
                else:
                    key_intv = "%s:%d-%s" % (intv[0],int(intv[1])-1,intv[2])
                    if hash_intervals.has_key(key_intv):
                        l_index_TN.append(key_intv)
                    else:
                        print "Problems with key: %s" % (key_intv)
                        
            df_all_P = df_all.loc[l_index_TP,:]
            df_all_N = df_all.loc[l_index_TN,:]
            
            l_zscore_P_s = calculate_zscore_distribution(df_all_P,sample)
            
            l_dosis_P_s = calculate_dosis_distribution(df_all_P.loc[:,sample])
            
            l_dosis_TP.extend(l_dosis_P_s)
            
            l_dosis_N_s = calculate_dosis_distribution(df_all_P.loc[:,df_all_P.columns != sample])
            
            l_dosis_TN.extend(l_dosis_N_s)
            
            l_zscore_TP.extend(l_zscore_P_s)
            
            l_zscore_N_s = calculate_zscore_distribution(df_all_N,sample)
            
            l_zscore_TN.extend(l_zscore_N_s)
            
            l_dosis_N_s = calculate_dosis_distribution(df_all_N.loc[:,sample])
            
            l_dosis_TN.extend(l_dosis_N_s)
            
        else: # all the intervals determined are FP
            
            l_dosis_N_s = calculate_dosis_distribution(df_all.loc[:,sample])
            
            l_dosis_TN.extend(l_dosis_N_s)
            
            l_zscore_N_s = calculate_zscore_distribution(df_all,sample)
            
            l_zscore_TN.extend(l_zscore_N_s)
                    
    bins_zscore = range(-11,11)
    
    histP, bin_edgesP = np.histogram(l_zscore_TP,bins_zscore)
    histN, bin_edgesN = np.histogram(l_zscore_TN,bins_zscore)
        
    print "\nDistribution of zscore in positive and negative intervals"
    print "pval\tfreq.P\tfreq.N"
    print "\n".join(map(lambda (b,f1,f2): "%1.3f\t%d\t%d" % (b,f1,f2) , zip(bin_edgesP[1:],histP,histN)))
    
    bins_dosis = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,2]
    
    histP, bin_edgesP = np.histogram(l_dosis_TP,bins_dosis)
    histN, bin_edgesN = np.histogram(l_dosis_TN,bins_dosis)
    
    print "\nDistribution of interval dosis in positive and negative intervals"
    print "pval\tfreq.P\tfreq.N"
    print "\n".join(map(lambda (b,f1,f2): "%1.3f\t%d\t%d" % (b,f1,f2) , zip(bin_edgesP[1:],histP,histN)))
    """
    
############################################################################333

if __name__=='__main__':
    
    run()
