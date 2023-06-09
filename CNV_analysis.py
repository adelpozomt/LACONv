'''
Created on March 2019

@author: adelpozo
'''

#!/usr/bin/python

import sys, re, shlex , os, string, urllib, time, math, random, subprocess, shutil

from operator import itemgetter

from itertools import groupby

import ConfigParser

from os import path as osp

import optparse

from os import path as osp

from subprocess import Popen , PIPE

import logging

import numpy

from scipy.stats import zmap

from pybedtools import BedTool

import pandas

import gzip

import glob

######################################################################

class OptionParser(optparse.OptionParser):

	def check_required (self, opt):

		option = self.get_option(opt)

		atrib = getattr(self.values, option.dest)
		
		if atrib is None:
#			self.error("%s option not supplied" % option)
			return False
		else:
			return True


#######################################################################

def read_cfg_file(cfg_filename):

	l_mandatory_sections = ['REFERENCE','OUTPUT','BDs','SOFTWARE','MARKERS']
	
	fi = open(cfg_filename,'r')

	config = ConfigParser.ConfigParser()
	config.readfp(fi)
	
	l_sections = config.sections()
	
	for s in l_mandatory_sections:
		if s not in l_sections:
			raise AttributeError('CNV_analysis.read_cfg_file: section not found in cfg file: %s' % (s))

	hash_cfg = {}
	
	for s in l_sections:
		for field in config.options(s):
			hash_cfg[field] = config.get(s,field)

	fi.close()

	return hash_cfg

#######################################################################

def get_chromosome_sex(sample,bam_file,bed_file_chrom,samtools_path,output_path,ref_fasta=None):

	fi = open(bed_file_chrom,'r')
	l_regions = map(lambda r: r.strip().split('\t'), fi.readlines())
	fi.close()

	chrY_regions = []
	chrX_regions = []
	
	l_pos_chr = groupby(map(lambda (i,r): (i,r[0]), enumerate(l_regions)),lambda x: x[1])
	
	for chrom,l_pos in l_pos_chr:
		if chrom == "chrX" or chrom == "X":
			chrX_regions = map(lambda (pos,c): l_regions[pos], l_pos)
		elif chrom == "chrY" or chrom == "Y":
			chrY_regions = map(lambda (pos,c): l_regions[pos], l_pos)
	
	mean_depth_chrX = 0
	
	sex_output_file = os.path.join(output_path,sample+'.gender')
	
	fo_g = open(sex_output_file,'w')
	fo_g.write("Interval\tTotal reads\tMean depth\n")
	
	for region in chrX_regions:
		
		chrm  = region[0]
		start = int(region[1])
		end   = int(region[2])

		if ref_fasta <> None:
			args_samtools = [samtools_path,'depth',bam_file,'--reference',ref_fasta,'-r',"%s:%d-%d" % (chrm,start,end)]
		else:
			args_samtools = [samtools_path,'depth',bam_file,'-r',"%s:%d-%d" % (chrm,start,end)]
		
		samtools_sal = Popen(args_samtools,stdin=PIPE, stdout=PIPE, stderr=PIPE,close_fds=True,bufsize=1)
		(str_stdout,logdata) = samtools_sal.communicate()
		samtools_sal.wait()
		
		if logdata <> "":
			if logdata.lower().find('fail')<>-1 or logdata.lower().find('error')<>-1:
				raise IOError('CNV_analysis.get_chromosome_sex: error in samtools depth: %s' % (logdata))
		
		total_depth = mean_depth_intv = 0
		
		if str_stdout <> "":
			
			total_depth = sum(map(lambda (chrm,pos,depth): int(depth), map(lambda l: l.split('\t'),str_stdout.strip().split('\n'))))
			
			mean_depth_intv = float(total_depth)/(end-start+1)
		
			mean_depth_chrX += mean_depth_intv
			
		fo_g.write("%s:%d-%d\t%d\t%1.2f\n" % (chrm,start,end,total_depth,mean_depth_intv))
	
	mean_depth_chrX = mean_depth_chrX/len(chrX_regions)
	
	mean_depth_chrY = 0
	
	for region in chrY_regions:
		
		chrm  = region[0]
		start = int(region[1])
		end   = int(region[2])
		
		if ref_fasta <> None:
			args_samtools = [samtools_path,'depth',bam_file,'--reference',ref_fasta,'-r',"%s:%d-%d" % (chrm,start,end)]
		else:
			args_samtools = [samtools_path,'depth',bam_file,'-r',"%s:%d-%d" % (chrm,start,end)]
	
		samtools_sal = Popen(args_samtools,stdin=PIPE, stdout=PIPE, stderr=PIPE,close_fds=True,bufsize=1)
		(str_stdout,logdata) = samtools_sal.communicate()
		samtools_sal.wait()
		
		total_depth = mean_depth_intv = 0
		
		if str_stdout <> "":
			
			total_depth = sum(map(lambda (chrm,pos,depth): int(depth), map(lambda l: l.split('\t'),str_stdout.strip().split('\n'))))
			
			mean_depth_intv = float(total_depth)/(end-start+1)
		
			mean_depth_chrY += mean_depth_intv
			
		fo_g.write("%s:%d-%d\t%d\t%1.2f\n" % (chrm,start,end,total_depth,mean_depth_intv))
			
	fo_g.close()
		
	mean_depth_chrY = mean_depth_chrY/len(chrY_regions)
	
	#print mean_depth_chrX,mean_depth_chrY
	
	if mean_depth_chrY >= 6 and mean_depth_chrX > 0:
		return 'XY'
	elif mean_depth_chrX > 2:
		return 'XX'
	else:
		return "UN"
	
#######################################################################

def __check_interval_integrity(cov_list,bed_analisis):
	
	### The intervals of bed analysis are the reference of the analysis
	### In case of being testing a subset of the intervals within the intervals of coverage files, the script must continue without any error
	### The function checks whether the intervals of bed analysis are included in the coverage file
	fi = open(bed_analisis,'r')
	l_ = map(lambda x: x.strip().split('\t'), fi.readlines())
	fi.close()
	
	if len(l_[0]) >= 6:
		hash_intervals = dict(map(lambda x: ((x[0],x[1],x[2]),(x[3],x[4],x[5])), l_))
	elif len(l_[0]) >= 5:
		hash_intervals = dict(map(lambda x: ((x[0],x[1],x[2]),(x[3],x[4],'-')), l_))
	elif len(l_[0]) >= 4:
		hash_intervals = dict(map(lambda x: ((x[0],x[1],x[2]),(x[3],'-','-')), l_))
	else:
		hash_intervals = dict(map(lambda x: ((x[0],x[1],x[2]),('-','-','-')), l_))
	
	cov_file = cov_list[0]
		
	fi2 = open(cov_file,'r')
	l_cov = map(lambda x: x.strip().split('\t'), fi2.readlines())
	fi2.close()
	
	l_intervals = []
	
	hash_intv_check = {}
	
	for i,line in enumerate(l_cov): ### Coverage files are always sorted!!!
			
		chrm  = line[0]
		start = line[1]
		end   = line[2]
			
		if hash_intervals.has_key((chrm,start,end)):
			l_intervals.append((i,chrm,start,end)+hash_intervals[chrm,start,end])
			hash_intv_check.setdefault((chrm,start,end),True)
			
	for (chrm,start,end) in sorted(hash_intervals.keys()):
		if not hash_intv_check.has_key((chrm,start,end)):
			raise IOError('CNV_analysis.__check_interval_integrity: the interval (%s:%s-%s) in coverage file does not appear in coverage files' % (chrm,start,end))
		
	return cov_list,l_intervals
	

def parse_cov_files(cov_list,bed_analisis):
	
	cov_list,l_intervals = __check_interval_integrity(cov_list,bed_analisis)
	
	num_intervals = len(l_intervals)
	num_samples   = len(cov_list)
	
	mat_coverage = numpy.zeros((num_intervals,num_samples),float)
	
	for index_s,cov_file in enumerate(cov_list):
		
		fi2 = open(cov_file,'r')
		l_cov = map(lambda x: x.strip().split('\t'), fi2.readlines())
		fi2.close()
				
		#l_values_cov = map(lambda line: float(line[-1]), l_cov)
		l_values_cov = map(lambda intv: float(l_cov[intv[0]][-1]), l_intervals)
			
		mat_coverage[:,index_s] = l_values_cov
	
	return mat_coverage,map(lambda intv: intv[1:], l_intervals)

def perform_coverage_control_mosdepth(bam_list,output_path,mosdepth_path,fasta_file,bed_analysis,logger,ref_fasta,**kwargs):
		
	l_cov = []
	
	correct_overlap = kwargs.get('overlap',False)
	Q_threshold     = kwargs.get('base_qual',15)
		
	for bam_file in bam_list:
		
		fileName,extName = os.path.splitext(os.path.basename(bam_file))
		fileName = os.path.join(output_path,fileName+"_Mosdepthcoverage")
		
		if extName <> ".cram":
			args = [mosdepth_path,fileName,'--no-per-base','-t','12','-b',bed_analysis,'-F','4','-T','0,20,50','-Q',"%d" % (Q_threshold),bam_file]
		else:
			args = [mosdepth_path,fileName,'--no-per-base','-t','12','-b',bed_analysis,'-f',ref_fasta,'-F','4','-T','0,20,50','-Q',"%d" % (Q_threshold),bam_file]
			
		if correct_overlap:
			args = args + ['--fast-mode']
		
		try:
			mosdepth_sal = Popen(args,stdin=PIPE, stdout=PIPE, stderr=PIPE,close_fds=True,bufsize=1)
			(output,logdata) = mosdepth_sal.communicate()
			mosdepth_sal.wait()
		except:
			raise RuntimeError("CNV_analysis.perform_coverage_control_mosdepth: Error in mosdepth with bam file: %s" % (bam_file))
	
		if logdata.lower().find('error') <> -1:
			raise RuntimeError("CNV_analysis.perform_coverage_control_mosdepth: Error in mosdepth with bam file: %s\n%s" % (bam_file,logdata))
			
		covFileName_i = fileName+'.regions.bed.gz'
		covFileName_o = fileName+".regions.bed"
		
		f = gzip.open(covFileName_i, 'rb')
		file_content = f.read()
		f.close()
			
		fo = open(covFileName_o,'w')
		fo.write(file_content)
		fo.close()
		
		for f in glob.glob(fileName+'.*gz.csi'):
			os.remove(f) 
		for f in glob.glob(fileName+'*thresholds*'):
			os.remove(f)
		for f in glob.glob(fileName+'.*mosdepth*txt'):
			os.remove(f)
		for f in glob.glob(fileName+'.*regions.bed.gz'):
			os.remove(f)
				
		l_cov.append(covFileName_o)
		
	return l_cov

def __divide(mean_cov,median_reads,gc_median):
	
	if gc_median > 0:
		return mean_cov*(median_reads/gc_median)
	else:
		return float(0)
		

def __normalize_coverage_v2(coverage_mat,l_intervals,sample,index_sample,gc_content_bed,logger=None):
	
	### Creating a bed of gc values
	bed_gc = BedTool(gc_content_bed)
	
	### Creating a bed from coverage values from the matrix
	cov_values = zip(l_intervals,list(coverage_mat[:,index_sample]))
	bed_cov = BedTool("\n".join(map(lambda (intv,c): "%s\t%s\t%s\t%1.2f" % (intv[0],intv[1],intv[2],c), cov_values)), from_string=True)
		
	### Step 0: It is created a dict whose key is a shared gc value. So it is extracted a gc value and the list of intervals with this gc content
	hash_gc_2_intv = {}
	
	for intv_gc in bed_gc:
		
		gc_value = "%1.2f" % (float(intv_gc.name))
		
		hash_gc_2_intv.setdefault(gc_value,[]).append((intv_gc.chrom,intv_gc.start,intv_gc.stop))
		
	### Step 1: Then, it is calculated the median value of the coverage of intervals with the same gc value
	median_cov_by_gc = {}
	
	for gc in sorted(hash_gc_2_intv.keys()):
		
		l_intv = sorted(hash_gc_2_intv[gc],key=lambda x: int(x[0]) if x[0].isdigit() else x[0]) # get all the intervals with the same gc
		
		gc_intv = BedTool("\n".join(map(lambda intv: "%s\t%d\t%d" % (intv[0],intv[1],intv[2]), l_intv)),from_string=True) # It is created a bed file with all intervals that share same gc value
		
		bed_gc_in_cov = bed_cov.intersect(gc_intv,wa=True) # It is extracted the coverage of these intervals
		
		l_cov = map(lambda b: float(b.name), bed_gc_in_cov) # Coverage values 
		
		median_cov_by_gc[gc] = numpy.median(l_cov)
	
	### Step 2: It is calculated the median value of: 1) Average values of coverage and 2) total number of reads in each interval
	median_cov   = numpy.median(map(lambda intv: float(intv.name), bed_cov))
	median_reads = numpy.median(map(lambda intv: (intv.end-intv.start)*float(intv.name), bed_cov))
	
	if median_cov == float(0):
		raise RuntimeError("CNV_analysis.__normalize_coverage_v2: The median coverage of sample %s is zero" % (sample))
	
	if  median_reads == float(0):
		raise RuntimeError("CNV_analysis.__normalize_coverage_v2: The median number of reads in sample %s is zero" % (sample))
	
	### Step 3: finally, normalization is performed
	l_all_cov_gc = map(lambda intv: (intv[0],int(intv[1]),int(intv[2]),float(intv[3]),"%1.2f" % (float(intv[7]))), bed_cov.intersect(bed_gc,wb=True))
	
	l_gc_median = map(lambda (chrom,s,e,mean_cov,gc): median_cov_by_gc[gc], l_all_cov_gc)
	l_ratio_depth_cov = map(lambda (chrom,s,e,mean_cov,gc): mean_cov/median_cov, l_all_cov_gc)
	l_reads_cov_norm = map(lambda (chrom,s,e,mean_cov,gc): mean_cov*(e-s), l_all_cov_gc)
	
	l_normalized = map(lambda ((chrom,s,e,mean_cov,gc),gc_median,reads_cov,ratio_depth_cov): (chrom,s,e,reads_cov,mean_cov,__divide(mean_cov,median_reads,gc_median),ratio_depth_cov), zip(l_all_cov_gc,l_gc_median,l_reads_cov_norm,l_ratio_depth_cov))
	
	median_reads_cov_norm = numpy.nanmedian(map(itemgetter(5),l_normalized))
	
	### It is returned a list whose items are formed by the tuple (interval,num_reads,avg_reads,ratio reads_norm with repect mean of reads_norm,ratio_depth_cov)
	if median_reads_cov_norm == float(0):
		raise RuntimeError("CNV_analysis.__normalize_coverage_v2: The median value of ratios in sample %s is zero" % (index_sample))
	
	return map(lambda intv: (intv[0],intv[1],intv[2],intv[3],intv[4],float(intv[5])/median_reads_cov_norm,intv[6]), l_normalized)
	

def __normalize_coverage(cov_file,gc_content_bed):
		
	"""
	# fields of cov file: Target	total_coverage	average_coverage	17NR2154_total_cvg	17NR2154_mean_cvg	17NR2154_granular_Q1	17NR2154_granular_median	17NR2154_granular_Q3	17NR2154_%_above_15
	# The fields of interest are: Target, total_coverage and average_coverage
	fi = open(cov_file,'r')
	l_ = map(lambda x: x.strip().split('\t'), fi.readlines())[1:]
	fi.close()
	"""
	
	# Bed of mosdepth coverage file
	# Fields: chr start end annot cov
	bed_cov = BedTool(cov_file)
	
	#l_intv_all = map(lambda intv: (intv.chrom,intv.start,intv.stop), bed_cov)
	
	# Bed of gc values
	bed_gc = BedTool(gc_content_bed)
	
	# It is created a dict whose key is a common gc value of a list of intervals with this same gc content
	hash_gc_2_intv = {}
	
	for intv_gc in bed_gc:
		
		gc_value = "%1.2f" % (float(intv_gc.name))
		
		hash_gc_2_intv.setdefault(gc_value,[]).append((intv_gc.chrom,intv_gc.start,intv_gc.stop))
			
	# It is calculated the median of coverage of the intervals with the same gc value
	median_cov_by_gc = {}
	gc_table = {} 
	
	for gc in sorted(hash_gc_2_intv.keys()):
		
		l_intv = sorted(hash_gc_2_intv[gc])
		
		l_cov = []
		
		for intv in l_intv:
			
			gc_intv  = BedTool("%s %d %d" % (intv[0],intv[1],intv[2]),from_string=True)
			
			bed_gc_in_cov = bed_cov.intersect(gc_intv,wa=True)
			
			map(lambda b: gc_table.setdefault((b.chrom,b.start,b.stop),gc), bed_gc_in_cov)
			
			cov_intv = map(lambda b: float(b.score), bed_gc_in_cov)
			l_cov.extend(cov_intv)
		
		median_cov_by_gc[gc] = numpy.median(l_cov)
		
	### It is calculated the median value of the average values in each interval
	median_cov   = numpy.median(map(lambda intv: float(intv.score), bed_cov))
	median_reads = numpy.median(map(lambda intv: (intv.end-intv.start)*float(intv.score), bed_cov))
	
	### The normalization is performed
	cov_normalized = {}
	
	l_intv_all_passed = []
	
	for interval in bed_cov:
		
		key_interval = (interval.chrom,interval.start,interval.stop)
		
		if not gc_table.has_key(key_interval):
			continue
		
		if key_interval not in l_intv_all_passed: 
			l_intv_all_passed.append(key_interval)
		
		gc_content = gc_table[key_interval]
		
		mean_cov   = float(interval.score)
		reads_cov  = mean_cov*(interval.stop-interval.start)  
		gc_median  = median_cov_by_gc[gc_content]
		
		try:
			ratio_depth_cov = mean_cov/median_cov
		except Exception.DivisionByZero:
			ratio_depth_cov = float(0)
		
		try:
			reads_cov_norm = mean_cov*(median_reads/gc_median) #num.reads.norm.avg
		except ZeroDivisionError:
			raise AttributeError('CNV_analysis.__normalize_coverage: Error in division. mean_cov=%s median_reads=%s gc_median=%s' % (str(mean_cov),str(median_reads),str(gc_median)))
		
		cov_normalized[key_interval] = (reads_cov,mean_cov,reads_cov_norm,ratio_depth_cov)
		
	median_reads_cov_norm = numpy.nanmedian(map(itemgetter(2),cov_normalized.values()))
	
	# it returns a list whose items are formed by the tuple (interval,num_reads,avg_reads,ratio reads_norm with repect mean of reads_norm,ratio_depth_cov)
	if median_reads_cov_norm > 0: ### Need to be checked!!!!!
		return map(lambda intv: (intv[0],intv[1],intv[2],cov_normalized[intv][0],cov_normalized[intv][1],float(cov_normalized[intv][2])/median_reads_cov_norm,cov_normalized[intv][3]), l_intv_all_passed)
	else:
		return map(lambda intv: (intv[0],intv[1],intv[2],cov_normalized[intv][0],cov_normalized[intv][1],0,cov_normalized[intv][3]), l_intv_all_passed)

def perform_coverage_gc_normalization(coverage_mat,l_intervals,l_samples,output_path,gatk_path,java_path,fasta_file,bed_analysis_filtered,logger):
	
	#### Run GCContentByInterval from GATK tool
	gc_content_template = os.path.join(output_path,"gc_content_template.txt")
	
	if os.path.exists(gc_content_template):
		os.remove(gc_content_template)
		
	args = [java_path,'-jar',gatk_path,'-T','GCContentByInterval','-R',fasta_file,'-L',bed_analysis_filtered,'-o',gc_content_template]
	
	gatk_sal = Popen(args,stdin=PIPE, stdout=PIPE, stderr=PIPE,close_fds=True,bufsize=1)
	(_,logdata) = gatk_sal.communicate()
	gatk_sal.wait()
	
	if logdata.lower().find('error')<>-1:
		raise RuntimeError('CNV_analysis.perform_coverage_normalization_gatk_v2: Error in GATK:\n%s' % (logdata))
	
	#### Parsing GATK output to convert intervals into bed file
	fi = open(gc_content_template,'r')
	l_ = map(lambda x: x.strip().split('\t'), fi.readlines())
	fi.close() 
	
	hash_gc = dict(map(lambda (i,x): ((i,x[0]),x[1]),enumerate(l_)))
	
	gc_content_template_bed = os.path.splitext(gc_content_template)[0]+'.bed'
			
	fo = open(gc_content_template_bed,'w')
	
	for i,intv in sorted(hash_gc.keys()):
		
		lintv = [intv.split(':')[0]]+intv.split(':')[1].split('-')
		
		gc_value = hash_gc[(i,intv)]
		
		fo.write('%s\t%s\t%s\t%s\n' % (lintv[0],lintv[1],lintv[2],gc_value))
		
	fo.close()
		
	#### Normalization sample by sample
	#### The normalized values are included in a matrix
	mat_coverage_norm = numpy.zeros((len(l_intervals),len(l_samples)),float)
	
	for i,sample in enumerate(l_samples):
		
		# Cov_normalized elements: interval,num_reads,avg_reads,ratio reads_norm with repect mean of reads_norm,ratio_depth_cov
		cov_normalized = __normalize_coverage_v2(coverage_mat,l_intervals,sample,i,gc_content_template_bed)
		
		#pos 3: total_num_reads, 4: total_avg_coverage, 5: reads_norm_by_gc_and_cov, 6: cov_norm
		
		mat_coverage_norm[:,i] = map(itemgetter(6),cov_normalized)
		
		logger.info("Normalized sample %s %d/%d" % (sample,i,len(l_samples)))
	
	return mat_coverage_norm
		
#######################################################################

def remove_samples_low_coverage(coverage_mat,l_intervals,l_samples,threshold_sample):
	
	l_samples_removed = []
	
	hash_coverage = {}
	
	for index_s,sample in enumerate(l_samples):
		
		mean_cov = numpy.mean(coverage_mat[:,index_s])
		median_cov = numpy.median(coverage_mat[:,index_s])
		
		hash_coverage[(index_s,sample)] = (mean_cov,median_cov)
		
		if mean_cov < threshold_sample or median_cov < 1:
			l_samples_removed.append((index_s,os.path.basename(sample),mean_cov,median_cov))
	
	coverage_mat_del = numpy.delete(coverage_mat, map(itemgetter(0),l_samples_removed), axis=1) 

	return coverage_mat_del,l_samples_removed,hash_coverage

#######################################################################

def filter_intervals(coverage_mat,l_intervals,l_index_male,threshold_interval,coverage_path,logger):
	
	# The goal is to remove intervals poorly capture (i.e. below the threshold_interval)
	l_intervals_discarded = []
	l_intervals_discarded_index = []
	l_intervals_filtr = []
	
	for i,intv in enumerate(l_intervals):
		
		if intv[0].upper().find('Y') <> -1:
			median_cov = numpy.median(coverage_mat[i,l_index_male])
		else:
			median_cov = numpy.median(coverage_mat[i,:])
		
		lenght_inv = int(intv[2])-int(intv[1])
		
		if lenght_inv < 3:
			l_intervals_discarded.append((intv,lenght_inv,median_cov))
			l_intervals_discarded_index.append(i)
		else:
			if median_cov < threshold_interval:
				l_intervals_discarded.append((intv,lenght_inv,median_cov))
				l_intervals_discarded_index.append(i)
			else:
				l_intervals_filtr.append(intv)
		
	coverage_mat_filtr = numpy.delete(coverage_mat, l_intervals_discarded_index, axis=0)
	
	return coverage_mat_filtr,l_intervals_filtr,l_intervals_discarded

#######################################################################

def tag_bed(l_intervals_filtered,l_intervals_all,output_path):
	
	#hash_intervals_del = dict(map(lambda (x,y,z): (x,(y,z)), l_intervals_filtered))
	hash_intervals_del = dict(map(lambda (x,y,z): ((x[0],x[1],x[2]),(y,z)), l_intervals_filtered))
	
	l_all_intervals_filtered = []
	
	for intv in l_intervals_all:
		
		key_intv = (intv[0],intv[1],intv[2])
		
		if hash_intervals_del.has_key(key_intv):
			l_all_intervals_filtered.append(intv + ('not considered',str(hash_intervals_del[key_intv])))
		else:
			l_all_intervals_filtered.append(intv + ('.','.'))
			
	
	filename_filtered = output_path+'.tagged_removed_intervals.bed'
	
	fo = open(filename_filtered,'w')
	fo.write("#Chrom\tStart\tEnd\tgene\texon\tNM\tNot considered\tInterval length/Mean cov.\n")
	fo.write("\n".join(map(lambda intv: "\t".join(intv), l_all_intervals_filtered)))
	fo.close()
	
	return filename_filtered

#######################################################################

def __remove_chrom_sexual(analysis_bed,cnv_output_path):
	
	bed_root_name = os.path.splitext(os.path.basename(analysis_bed))[0]
	bed_analysis_filtered = os.path.join(cnv_output_path,bed_root_name+'.without_sex_chrom.bed')
	
	bed_intv = BedTool(analysis_bed)
	
	l_intervals = map(lambda intv: (intv[0],intv[1],intv[2],intv[3],intv[4],intv[5]), bed_intv)
	
	l_intervals_chrX = filter(lambda interval: interval[0].find("X")<>-1, l_intervals)
	l_intervals_chrY = filter(lambda interval: interval[0].find("Y")<>-1, l_intervals)
	
	bed_intv.subtract(BedTool(l_intervals_chrX+l_intervals_chrY)).saveas(bed_analysis_filtered)
	
	return bed_analysis_filtered
	

def __is_autosome(chrom):
	
	chrom_filtr = chrom.replace('chr','').split('_')[0]
	
	if chrom_filtr.isdigit():
		return True
	else:
		if chrom_filtr.find('X')<>-1 and chrom_filtr.find('Y')<>-1:
			return True
	
	return False 
		
		
def extract_chrom_autosomes(mat_all_ratios,l_intervals):
		
	l_intervals_autosomes = filter(lambda (i,interval): __is_autosome(interval[0]), enumerate(l_intervals))
		
	return map(itemgetter(0),l_intervals_autosomes)
	

def extract_chrom_sexual(mat_all_ratios,l_intervals):
	
	l_intervals_chrX = filter(lambda (i,interval): interval[0].find("X")<>-1, enumerate(l_intervals))
		
	l_intervals_chrY = filter(lambda (i,interval): interval[0].find("Y")<>-1, enumerate(l_intervals))
	
	return map(itemgetter(0),l_intervals_chrX),map(itemgetter(0),l_intervals_chrY)

#######################################################################
	
def calculate_dosis(coverage_mat_norm,l_intervals,l_samples,cnv_output_path,logger):
	
	num_intervals = len(l_intervals)
	num_samples   = len(l_samples)
	
	# normalization intra-sample
	norm_depth_mat = numpy.empty((num_intervals,num_samples),float)
	
	for i in range(num_samples): # index 'i' labels each sample
		
		l_rates_sample = coverage_mat_norm[:,i]
		
		mean_rate_sample = numpy.nanmean(l_rates_sample)
		
		norm_depth_mat[:,i] = map(lambda r: float(r)/mean_rate_sample, l_rates_sample)
		
	# normalization inter-sample (by interval)
	for j in range(num_intervals): # index 'i' labels each interval
		
		mean_rate_interval = numpy.nanmean(norm_depth_mat[j,:])
		
		try:
			norm_depth_mat[j,:] = norm_depth_mat[j,:]/mean_rate_interval
		except:
			norm_depth_mat[j,:] = 0
		
	
	#mat_filename = os.path.join(cnv_output_path,"dosis_ratios_all.mat")
	mat_filename = cnv_output_path
	
	fo = open(mat_filename,'w')
	fo.write("\t"+'\t'.join(l_samples)+'\n')
	for i,intv in enumerate(l_intervals):
		l_ = map(lambda x: "%1.3f" % (x), norm_depth_mat[i,:])
		fo.write("%s:%s-%s\t%s\n" % (intv[0],intv[1],intv[2],"\t".join(l_)))
	fo.close()
			
	return norm_depth_mat
	
#######################################################################

def __calculate_zscore(mat_all_ratios,index_sample):
	
	mat_sample = mat_all_ratios[:,[index_sample]]
	
	selector_samples = [x for x in range(mat_all_ratios.shape[1]) if x != index_sample]		 
	
	mat_remaining = mat_all_ratios[:, selector_samples]
	
	zscore_array = zmap(mat_sample,mat_remaining,axis=1)
	
	return zscore_array
	
def establish_deletions_and_duplications_zscore(mat_all_ratios,l_intervals,l_samples,threshold_del,threshold_dup,sig_level_del,sig_level_dup,logger):
	
	hash_index = dict(zip(l_intervals,range(len(l_intervals))))
	
	# First, the matrix of zscores is calculated
	mat_zscore = numpy.zeros((len(l_intervals),len(l_samples)),float)
	
	for i_s,sample in enumerate(l_samples):
		
		zscore_array = __calculate_zscore(mat_all_ratios,i_s)
		
		mat_zscore[:,i_s] = list(zscore_array) 
	
	# The, the values are inspected to get the intervals that supply thresholds contrains. The results are stored in a dict structure
	hash_intervals_cnv = {}
	
	for i_s,sample in enumerate(l_samples):
		
		hash_intervals_cnv[sample] = []
		
		mat_sample = mat_all_ratios[:,[i_s]]
			
		mat_remaining = mat_all_ratios[:, [x for x in range(mat_all_ratios.shape[1]) if x != i_s]]
		
		l_intervals_index_dosis_del = numpy.where((mat_sample<=threshold_del))[0]
		l_intervals_index_dosis_dup = numpy.where((mat_sample>=threshold_dup))[0]
		
		if l_intervals_index_dosis_del.size <> 0:
			
			mat_zscore_sample = mat_zscore[l_intervals_index_dosis_del,i_s]
			
			mat_remaining_del = mat_remaining[l_intervals_index_dosis_del,:]
			
			l_intervals_dosis = map(lambda i: l_intervals[i],list(l_intervals_index_dosis_del))
		
			l_intervals_index_zscore = numpy.where((mat_zscore_sample<=sig_level_del))[0]
			
			if l_intervals_index_zscore.size <> 0:
				
				for i_intv in numpy.nditer(l_intervals_index_zscore):
					
					intv = l_intervals_dosis[i_intv]
					
					i2_intv = hash_index[intv]
					
					(chrm,start,end,annot1,annot2,annot3) = l_intervals[i2_intv]
					
					dosis_sample = mat_sample[i2_intv]
					
					zscore_sample = mat_zscore_sample[i_intv]
															
					type_ = "CNV loss"
					
					other =  "|".join(map(lambda r: "%1.2f" % (r), mat_remaining_del[i_intv,:]))
					
					hash_intervals_cnv[sample].append((chrm,start,end,annot1,annot2,annot3,dosis_sample,zscore_sample,"*",type_,other))   
			
		if l_intervals_index_dosis_dup.size <> 0:
			
			mat_zscore_sample = mat_zscore[l_intervals_index_dosis_dup,i_s]
			
			mat_remaining_dup = mat_remaining[l_intervals_index_dosis_dup,:]
			
			l_intervals_dosis = map(lambda i: l_intervals[i],list(l_intervals_index_dosis_dup))
		
			l_intervals_index_zscore = numpy.where((mat_zscore_sample>=sig_level_dup))[0]
			
			if l_intervals_index_zscore.size <> 0:
				
				for i_intv in numpy.nditer(l_intervals_index_zscore):
					
					intv = l_intervals_dosis[i_intv]
					
					i2_intv = hash_index[intv]
					
					(chrm,start,end,annot1,annot2,annot3) = l_intervals[i2_intv]
					
					dosis_sample = mat_sample[i2_intv]
					
					zscore_sample = mat_zscore_sample[i_intv]
					
					type_ = "CNV gain"
					
					other =  "|".join(map(lambda r: "%1.2f" % (r), mat_remaining_dup[i_intv,:]))
					
					hash_intervals_cnv[sample].append((chrm,start,end,annot1,annot2,annot3,dosis_sample,zscore_sample,"*",type_,other)) 
			
	return hash_intervals_cnv,mat_zscore
	
#######################################################################

def calculate_overlapping(hash_intervals):
	
	hash_intervals_annot = {}
	
	hash_counts_loss = {}
	hash_counts_gain = {}
	
	for sample in sorted(hash_intervals.keys()):
			
			for (chrom,s,e,annot,annot2,annot3,rate,pval,sig,type_,other) in hash_intervals[sample]:
				
				if sig == '*' and type_ == "CNV loss":
					hash_counts_loss.setdefault((chrom,s,e,'loss'),[]).append(sample)
				elif sig == '*' and type_ == "CNV gain":
					hash_counts_gain.setdefault((chrom,s,e,'gain'),[]).append(sample)
					
	for sample in sorted(hash_intervals.keys()):
		
		hash_intervals_annot[sample] = []
		
		for (chrom,s,e,annot,annot2,annot3,rate,pval,sig,type_,other) in hash_intervals[sample]:
			
			count_overlap_less = "-"
			count_overlap_gain = "-"
			
			samples_overlap_less = "-"
			samples_overlap_gain = "-"
			
			if hash_counts_loss.has_key((chrom,s,e,'loss')) and sig == "*":
				
				l_samples_l = hash_counts_loss[(chrom,s,e,'loss')]
				
				count_overlap_less   = "%d" % (len(l_samples_l))
				samples_overlap_less = '|'.join(l_samples_l)
				
			if hash_counts_gain.has_key((chrom,s,e,'gain')) and sig == "*":
				
				l_samples_g = hash_counts_gain[(chrom,s,e,'gain')]
				
				count_overlap_gain   = "%d" % (len(l_samples_g))
				samples_overlap_gain = '|'.join(l_samples_g)
				
			hash_intervals_annot[sample].append((chrom,s,e,annot,annot2,annot3,sample,rate,pval,sig,type_,other,count_overlap_less,samples_overlap_less,count_overlap_gain,samples_overlap_gain))
	
	return hash_intervals_annot

#######################################################################

def write_reports_filtered(sample,cnv_output_path,tsv_analysis):
						
	output_xlsx_file = os.path.join(cnv_output_path,"%s_laconv_cnv.xlsx" % (sample))

	writer = pandas.ExcelWriter(output_xlsx_file, engine='xlsxwriter')
	#workbook = writer.book
	
	### writing configuration
	conf_filename = os.path.join(cnv_output_path,"LACONv.cfg")
	
	if not os.path.exists(conf_filename):
		raise IOError('CNV_analysis.write_reports: The conf file does not exist. %s' % (conf_filename))
	
	df_cfg = pandas.read_csv(conf_filename, delimiter='\t', na_values=['.'], header=0, index_col=None)
	df_cfg.to_excel(writer, index=False, sheet_name='LACONv config')
	
	### writing chromosomal gender
	tmp_path = os.path.join(cnv_output_path,'tmp')
	
	gender_filename = os.path.join(tmp_path,sample+'.gender')
	
	if not os.path.exists(gender_filename):
		raise IOError('CNV_analysis.write_reports: The gender file does not exist. %s' % (gender_filename))
	
	df_gender = pandas.read_csv(gender_filename, delimiter='\t', na_values=['.'], header=0, index_col=None,names=['Interval','Total reads','Mean depth'])
	df_gender.to_excel(writer, index=False, sheet_name='chromosomal sex markers')
	
	### writing intervals with cnvs
	df_sig = pandas.DataFrame(columns=['sample not considered for the analysis'])
	df_sig.to_excel(writer, index=False, sheet_name='CNV intervals')
	
	### writing all ratios
	df_all = pandas.DataFrame(columns=['sample not considered for the analysis'])
	df_all.to_excel(writer, index=False, sheet_name='Ratios all intervals')
	
	### writing removed intervals
	if not os.path.exists(tsv_analysis):
		raise IOError('CNV_analysis.write_reports: The tsv file does not exist. %s' % (tsv_analysis))
	
	df_intv = pandas.read_csv(tsv_analysis, delimiter='\t', na_values=['.'], header=0, index_col=None)
	df_intv.to_excel(writer, index=False, sheet_name='Intervals analyzed')
	
	writer.save()
	writer.close()
	
	return output_xlsx_file


def write_reports(sample,cnv_output_path,tsv_analysis):
	
	output_xlsx_file = os.path.join(cnv_output_path,"%s_CNV.xlsx" % (sample))

	writer = pandas.ExcelWriter(output_xlsx_file, engine='xlsxwriter')
	#workbook = writer.book
	
	### writing configuration
	conf_filename = os.path.join(cnv_output_path,"LACONv.cfg")
	
	if not os.path.exists(conf_filename):
		raise IOError('CNV_analysis.write_reports: The conf file does not exist. %s' % (conf_filename))
	
	df_cfg = pandas.read_csv(conf_filename, delimiter='\t', na_values=['.'], header=0, index_col=None)
	df_cfg.to_excel(writer, index=False, sheet_name='LACONv config')
	
	### writing chromosomal gender
	tmp_path = os.path.join(cnv_output_path,'tmp')
	
	gender_filename = os.path.join(tmp_path,sample+'.gender')
	
	if not os.path.exists(gender_filename):
		raise IOError('CNV_analysis.write_reports: The gender file does not exist. %s' % (gender_filename))
	
	df_gender = pandas.read_csv(gender_filename, delimiter='\t', na_values=['.'], header=0, index_col=None,names=['Interval','Total reads','Mean depth'])
	df_gender.to_excel(writer, index=False, sheet_name='chromosomal sex markers')
	
	### writing intervals with cnvs
	bed_sig_ratios = os.path.join(cnv_output_path,"%s.cnv.bed" % (sample))
	
	num_cnv = len(open(bed_sig_ratios).readlines())
		
	if not os.path.exists(bed_sig_ratios):
		df_sig = pandas.DataFrame(columns=['no CNVs determined in this sample'])
	elif num_cnv == 0:
		df_sig = pandas.DataFrame(columns=['no CNVs determined in this sample'])
	else:
		df_sig = pandas.read_csv(bed_sig_ratios, delimiter='\t', na_values=['.'], header=None, index_col=None,names=['Chrom','Start','End','Gene','Exon','Transcript','Sample','Dosis','zscore','is CNV','CNV type','Num. Total Samples','Num. Samples CNV loss','Samples CNV loss','Num. Samples CNV gain','Samples CNV gain','Dosis pool'])
	
	df_sig.to_excel(writer, index=False, sheet_name='CNV intervals')
	
	### writing all ratios
	bed_all_ratios = os.path.join(cnv_output_path,"%s_all_intervals.bed" % (sample))
		
	if not os.path.exists(bed_all_ratios):
		raise IOError('CNV_analysis.write_reports: The bed file does not exist. %s' % (bed_all_ratios))
	
	if os.stat(bed_all_ratios).st_size <> 0:
		df_all = pandas.read_csv(bed_all_ratios, delimiter='\t', na_values=['.'], header=None, index_col=None,names=['Chrom','Start','End','Annot','Dosis','Zscore'])
		df_all.to_excel(writer, index=False, sheet_name='Ratios all intervals')
	
	### writing removed intervals
	if not os.path.exists(tsv_analysis):
		raise IOError('CNV_analysis.write_reports: The tsv file does not exist. %s' % (tsv_analysis))
	
	df_intv = pandas.read_csv(tsv_analysis, delimiter='\t', na_values=['.'], header=0, index_col=None)
	df_intv.to_excel(writer, index=False, sheet_name='Intervals analyzed')
	
	writer.save()
	
	writer.close()
	
	return output_xlsx_file

#######################################################################

def __sum_overlap(c_less,c_gain):
	
	if c_less == '-':
		c_less = 0
	if c_gain == '-':
		c_gain = 0
		
	return int(c_less)+int(c_gain)

def __sort_key(x):
	
	chrom = x[0]

	if chrom[3:].isdigit():
		return int(chrom[3:]),int(x[1]),int(x[2])
	else:
		return x[0],int(x[1]),int(x[2])
	
#######################################################################

def run(argv=None):

	if argv is None: argv = sys.argv    

	parser = OptionParser(add_help_option=True,description="The script performs CNV estimation within the regions of interest")
	
	parser.add_option("--cfg",default=None,help="Config file with the complete information of the target regions and paths of the files needed for the calling",dest="f_cfg")
	parser.add_option("--i",default=None,help="File with a list of bam files",dest="f_listbam")
	parser.add_option("--f",default=None,action="append",help="Bam file",dest="f_bam")
	parser.add_option("--b",default=None,help="Bed file with the intervals to perform coverage",dest="f_bed")
	parser.add_option("--o",default=None,help="Output path. The results are written in subdirectory laconv/",dest="output_path")
	parser.add_option("--no-cov",action="store_true",default=False,help="Flag to disconnect coverage calculus because the files are already generated",dest="no_cov")
	parser.add_option("--no-sexual",action="store_true",default=False,help="Flag to disconnect analysis of CNVs in chromosome Y",dest="no_sexual")
	parser.add_option("--no-overlap",action="store_true",default=False,help="Flag to correct mate overlaps in coverage counting",dest="no_overlap")
	parser.add_option("--qbase",default=15,help="minimum base quality. Default=15",dest="qualbase_threshold")
	parser.add_option("--t-del",default=-2,help="zscore threshold for deletions. Default=-2",dest="zscore_threshold_del")
	parser.add_option("--t-dup",default=2,help="zscore threshold for duplications. Default=2",dest="zscore_threshold_dup")
	parser.add_option("--rdel",default=0.6,help="Deletion rate threshold. Default=0.8",dest="del_threshold")
	parser.add_option("--rdup",default=1.4,help="Duplication rate threshold. Default=1.4",dest="dup_threshold")
	parser.add_option("--min_depth_intv",default=15,help="Mean depth by interval. Default=10",dest="depth_threshold_interval")
	parser.add_option("--min_depth_sample",default=20,help="Mean depth by sample. Default=20",dest="depth_threshold_sample")
	parser.add_option("--clean",action="store_true",default=False,help="Flag to remove the temporal directory",dest="clean")


	# Se leen las opciones aportadas por el usuario
	(options, args) = parser.parse_args(argv[1:])

	if len(argv) == 1:
		sys.exit(0)

	if not parser.check_required("--cfg"):
		raise IOError('CNV_analysis: The cfg file has not been provided')
	if not parser.check_required("--b"):
		raise IOError('CNV_analysis: The bed file has not been provided')
			
	cfg_file = options.f_cfg

	if not os.path.exists(cfg_file):
		raise IOError('CNV_analysis: The cfg file does not exist\n%s' % (cfg_file))

	hash_cfg = read_cfg_file(cfg_file)
		
	ref_fasta = hash_cfg.get('ref_fasta','None')
	
	if not os.path.isfile(ref_fasta):
		raise IOError('CNV_analysis: The file does not exist. %s' % ref_fasta)
	
	analysis_bed = options.f_bed
	
	if not os.path.isfile(analysis_bed):
		raise IOError('CNV_analysis: The file does not exist. %s' % analysis_bed)
	
	java_path = hash_cfg.get("java_path",'None')
	
	if not os.path.exists(java_path):
		raise IOError('CNV_analysis: The path of java does not exist: %s' % (java_path))

	samtools_path = hash_cfg.get("samtools_path",'None')
	
	if not os.path.exists(samtools_path):
		raise IOError('CNV_analysis: The path of samtools tool does not exist: %s' % (samtools_path))
	
	gatk_path = hash_cfg.get('gatk_path','None')
	
	if not os.path.exists(gatk_path):
		raise IOError('CNV_analysis: The path of gatk tool does not exist: %s' % (gatk_path))
	
	mosdepth_path = hash_cfg.get("mosdepth_path",'None')

	if not os.path.exists(mosdepth_path):
		raise IOError('CNV_analysis: The path of mosdepth tool does not exist: %s' % (mosdepth_path))

	bed_file_chrom = hash_cfg.get('chrom_sexual','None')
	
	if not os.path.exists(bed_file_chrom):
		raise IOError("CNV_analysis: The bed of chromosomal sex markers has not be provided\n%s" % bed_file_chrom)
	
	cfg_output_path = hash_cfg.get('cnv_path',None)
	
	if cfg_output_path <> None and options.output_path == None:
		cnv_output_root_path = cfg_output_path
	elif options.output_path <> None:
		cnv_output_root_path = options.output_path 
	else:
		raise IOError('CNV_analysis: The path of CNVs has not been provided')
	
	if not os.path.exists(cnv_output_root_path):
		raise IOError('CNV_analysis: The path of CNVs does not exist. %s' % cnv_output_root_path)
		
	cnv_output_path = os.path.join(cnv_output_root_path,'laconv')

	if not os.path.exists(cnv_output_path):			
		os.mkdir(cnv_output_path)
		
	tmp_output_path = os.path.join(cnv_output_path,'tmp')
	
	if not os.path.exists(tmp_output_path):			
		os.mkdir(tmp_output_path)
	
	coverage_path = os.path.join(cnv_output_root_path,'coverage')

	if not os.path.exists(coverage_path):
		os.mkdir(coverage_path)
	
	if not os.path.isfile(samtools_path):
		raise IOError('CNV_analysis: The file does not exist. %s' % samtools_path)

	if not os.path.isfile(gatk_path):
		raise IOError('CNV_analysis: The file does not exist. %s' % gatk_path)
	
	#Configure logger
	formatter = logging.Formatter('%(asctime)s - %(module)s - %(levelname)s - %(message)s')
	console = logging.StreamHandler()
	console.setFormatter(formatter)
	console.setLevel(logging.INFO)

	log_file = os.path.join(cnv_output_path,"LACONv.log")

	if os.path.exists(log_file):
		os.remove(log_file)

	logging.basicConfig(filename=log_file,level=logging.INFO,format='%(asctime)s,%(msecs)d %(levelname)-8s %(message)s',datefmt='%Y-%m-%d:%H:%M:%S')

	logger = logging.getLogger("preprocess")
	
	# set log Handler just once
	if ( len( logger.handlers ) == 0 ):
		logger.addHandler(console)	

	logger.info("LACONv last modification: 04/05/2021\n")
	
	#### Global parameters
	DEPTH_THRESHOLD_BY_INTERVAL = int(options.depth_threshold_interval)
	DEPTH_THRESHOLD_BY_SAMPLE = int(options.depth_threshold_sample)
	SIG_LEVEL_DEL   = float(options.zscore_threshold_del)
	SIG_LEVEL_DUP   = float(options.zscore_threshold_dup)
	THRESHOLD_DEL   = float(options.del_threshold)
	THRESHOLD_DUP   = float(options.dup_threshold)
	QUAL_BASE       = int(options.qualbase_threshold)
	
	logger.info("LACONv parameter configuration:\n")
	if options.no_cov:
		logger.info("Minimum base quality in coverage calculus: not applied")
		logger.info("Considered overlap correction of read mate: not applied")
	else:
		logger.info("Minimum base quality in coverage calculus: %d" % (QUAL_BASE))
		logger.info("Considered overlap correction of read mate: %s" % (options.no_overlap))
	logger.info("zscore threshold to establish CNVs (duplications): %1.2f" % (SIG_LEVEL_DUP))
	logger.info("Minimal depth by sample to be considered: %dX" % (DEPTH_THRESHOLD_BY_SAMPLE))
	logger.info("Minimal depth in genomic intervals to be considered: %dX" % (DEPTH_THRESHOLD_BY_INTERVAL))
	logger.info("Dosis rate threshold for deletions: %1.2f" % (THRESHOLD_DEL))
	logger.info("Dosis rate threshold for duplications: %1.2f" % (THRESHOLD_DUP))
	logger.info("zscore threshold to establish CNVs (deletions): %1.2f" % (SIG_LEVEL_DEL))
	logger.info("zscore threshold to establish CNVs (duplications): %1.2f\n" % (SIG_LEVEL_DUP))
	
	config_filename = os.path.join(cnv_output_path,"LACONv.cfg")
	
	fcfg = open(config_filename,"w")
	fcfg.write("LACONv.1.0 configuration parameter\tValue\n")
	if options.no_cov:
		fcfg.write("Minimum base quality in coverage calculus: not applied\n")
		fcfg.write("Considered overlap correction of read mate: not applied\n")
	else:
		fcfg.write("Minimum base quality in coverage calculus: %d\n" % (QUAL_BASE))
		fcfg.write("Considered overlap correction of read mate: %s\n" % (options.no_overlap))
	fcfg.write("Minimal depth by sample to be considered:\t%dX\n" % (DEPTH_THRESHOLD_BY_SAMPLE))
	fcfg.write("Minimal depth in genomic intervals to be considered:\t%dX\n" % (DEPTH_THRESHOLD_BY_INTERVAL))
	fcfg.write("Dosis rate threshold for deletions:\t%1.2f\n" % (THRESHOLD_DEL))
	fcfg.write("Dosis rate threshold for duplications:\t%1.2f\n" % (THRESHOLD_DUP))
	fcfg.write("zscore threshold to establish CNVs (deletions):\t%1.3f\n" % (SIG_LEVEL_DEL))
	fcfg.write("zscore threshold to establish CNVs (duplications):\t%1.3f\n" % (SIG_LEVEL_DUP))
	fcfg.close()
	
	logger.info("Written parameters to configuration file: %s" % (config_filename))
	
	#### List of final files
	l_final_files = []
			
	#### First, it is checked if the bams exit
	l_bams	= []
	l_samples = []
	
	if options.f_listbam <> None:
		
		if not os.path.exists(options.f_listbam):
			raise IOError('CNV_analysis: The file with the list of bam files does not exist')
		
		fi = open(options.f_listbam,'r')
		l_bams_raw = map(lambda l: l.strip(), fi.readlines())
		fi.close()
	
	elif options.f_bam <> None:
		l_bams_raw = options.f_bam 
	else:
		raise IOError("CNV_analysis: No bam file has been provided")
			
	for bam_file in l_bams_raw:
			
		if bam_file[0] == '#':
			continue
			
		if not os.path.exists(bam_file) and not options.no_cov:
			raise IOError("CNV_analysis: The bam file does not exist:\n%s" % (bam_file))
			
		l_bams.append(bam_file)
			
		basename_bam = os.path.basename(bam_file)  
			
		pos = basename_bam.find("_align.realign.recal.bam")
			
		if pos <> -1:
			l_samples.append(basename_bam[:pos])
		else:
			l_samples.append(os.path.splitext(basename_bam)[0])
						
	logger.info("The paths of bams files are checked\n")
		
	#### Checking of chromosomal sex
	hash_gender = {}
	
	if options.no_sexual:
		
		analysis_bed_no_sex = __remove_chrom_sexual(analysis_bed,tmp_output_path)
		analysis_bed =  analysis_bed_no_sex
		
		logger.info("Chromosomal sex of bam files are not checked due to be disconnected")
		logger.info("Created a bam file with no sexual chromosomes: %s" % (analysis_bed_no_sex))
		
	else:
		logger.info("Checking the chromosomal sex of bam files")
	
		for i,(bam_file,s) in enumerate(zip(l_bams,l_samples)):
		
			if os.path.splitext(bam_file)[1] == '.cram':
				chrom_sex = get_chromosome_sex(s,bam_file,bed_file_chrom,samtools_path,tmp_output_path,ref_fasta)
			else:		
				chrom_sex = get_chromosome_sex(s,bam_file,bed_file_chrom,samtools_path,tmp_output_path)
		
			if chrom_sex == "XX":
				hash_gender.setdefault((i,s),"F")
			elif chrom_sex == "XY":
				hash_gender.setdefault((i,s),"M")
			else:
				hash_gender.setdefault((i,s),"X")
			
		logger.info("The chromosomal sex of samples are calculated:\n")
		logger.info("\n".join(map(lambda k: "%s-->%s" % (k[1],hash_gender[k]), sorted(hash_gender.keys()))))
										
	if options.no_cov:
		
		l_cov_output = []
		
		for bam_file in l_bams:
		
			fileName = os.path.splitext(os.path.basename(bam_file))[0]
			cov_fileName = os.path.join(coverage_path,fileName+"_Mosdepthcoverage.regions.bed")
			if os.path.exists(cov_fileName):
				l_cov_output.append(cov_fileName)
			else:
				raise IOError("CNV_analysis: The coverage file does not exist. %s" % (cov_fileName))
	
	else:
		# Before calling for CNV, it is necessary to create GATK coverage files , normalize them and do the estimation CNV 
		# GATK coverage average calling			
		logger.info("Starting coverage analysis...")
		
		l_cov_output = perform_coverage_control_mosdepth(l_bams,coverage_path,mosdepth_path,ref_fasta,analysis_bed,logger,ref_fasta,overlap=options.no_overlap)
	
		logger.info("Coverage finished!!")
		
	logger.info("Parsing coverage files...")
	
	coverage_mat,l_intervals = parse_cov_files(l_cov_output,analysis_bed)
	
	# Sample filtering: removing samples with a mean and median coverage below a threshold
	logger.info("Excluding the samples with mean coverage below the value %d..." % (DEPTH_THRESHOLD_BY_SAMPLE))
	
	coverage_mat_filtered,l_samples_removed,hash_coverage = remove_samples_low_coverage(coverage_mat,l_intervals,l_samples,DEPTH_THRESHOLD_BY_SAMPLE)
	
	filename_cov = os.path.join(coverage_path,'sample_coverage.tsv')
	fo = open(filename_cov,'w')
	fo.write('Sample\tmean coverage\tMedian coverage\n')
	fo.write('\n'.join(map(lambda (i,s): "%s\t%1.2f\t%1.2f" % (s,hash_coverage[(i,s)][0],hash_coverage[(i,s)][1]), sorted(hash_coverage.keys()))))
	fo.close()
	
	logger.info("Saving sample coverage in file: %s" % (filename_cov))
	
	logger.info("Done!. It has been excluded from the analysis %d samples due to low coverage.\nSamples:\n%s" % (len(l_samples_removed),"\n".join(map(lambda x: "%s (%1.2fX)" % (x[1],x[2]), l_samples_removed))))

	l_samples_filtered = map(lambda s: s, l_samples)
	
	map(lambda (i,s,c,m): l_samples_filtered.remove(s),  l_samples_removed)
	
	if len(l_samples_filtered) <= 1:
		raise RuntimeError("CNV_analysis: There are not enough samples to perform CNV analysis. Number=%s" % (len(l_samples_filtered)))
			
	# Removing those intervals with a median value less than threshold
	logger.info("Removing intervals with median values of depth across the samples below the threshold %dX or less than 2 bp" % (DEPTH_THRESHOLD_BY_INTERVAL))
	
	l_female_index = []
	l_male_index   = []
	
	l_samples_filtered_f = []
	l_samples_filtered_m = []
			
	if hash_gender.keys() <> []:
		hash_sample_2_index = dict(map(lambda (i,s): (s,i), hash_gender.keys()))
		l_samples_filtered_index = map(lambda s: (hash_sample_2_index[s],s), l_samples_filtered)
				
		for j,(i,s) in enumerate(l_samples_filtered_index):
			if hash_gender[(i,s)]=="M":
				l_male_index.append(j)
				l_samples_filtered_m.append(s)
			elif hash_gender[(i,s)]=="F":
				l_female_index.append(j)
				l_samples_filtered_f.append(s)
					
	coverage_mat_filtered_2,l_intervals_filtered,l_intervals_removed = filter_intervals(coverage_mat_filtered,l_intervals,l_male_index,DEPTH_THRESHOLD_BY_INTERVAL,coverage_path,logger)
			
	intervals_remove_filename = os.path.join(cnv_output_path,"intervals_removed_by_coverage.tsv")
	l_intervals_removed_str = map(lambda (intv,len_,cov): "%s\t%s\t%s\t%s\t%d\t%1.2fX" % (intv[0],intv[1],intv[2],intv[3],len_,cov), l_intervals_removed)
	fo = open(intervals_remove_filename,"w")
	fo.write("Chrom\tStart\tEnd\tAnnot.\tMean cov.\n")
	fo.write("\n".join(l_intervals_removed_str))
	fo.close()
	
	l_final_files.append(intervals_remove_filename)
	
	logger.info("Done!. It has been deleted %d intervals. Intervals are written in file %s" % (len(l_intervals_removed),intervals_remove_filename))
	
	# GC Normalization
	logger.info("Starting GC normalization of depth of coverage...")
	
	bed_root_name = os.path.splitext(os.path.basename(analysis_bed))[0]
	
	tagged_bed = os.path.join(cnv_output_path,bed_root_name)
	
	analysis_bed_tagged = tag_bed(l_intervals_removed,l_intervals,tagged_bed)
	
	l_final_files.append(analysis_bed_tagged)
	
	### It is created a new bed file with the intervals that are remaining after filtering
	bed_analysis_filtered = os.path.join(cnv_output_path,bed_root_name+'.without_filtered_intervals.bed')
	bed_filtered = BedTool("\n".join(map(lambda intv: "\t".join(intv), l_intervals_filtered)), from_string=True)
	bed_filtered.saveas(bed_analysis_filtered)
	
	l_final_files.append(bed_analysis_filtered)

	coverage_mat_filtered_norm = perform_coverage_gc_normalization(coverage_mat_filtered_2,l_intervals_filtered,l_samples_filtered,coverage_path,gatk_path,java_path,ref_fasta,bed_analysis_filtered,logger)
	
	# Splitting chrom
	logger.info("Splitting the intervals into autosomes and sexual chromosomes\n")
	
	l_intervals_filtered_autosomes_index = extract_chrom_autosomes(coverage_mat_filtered_norm,l_intervals_filtered)
	
	l_intervals_filtered_chrX_index = l_intervals_filtered_chrY_index = []
	
	if not options.no_sexual:
		l_intervals_filtered_chrX_index,l_intervals_filtered_chrY_index = extract_chrom_sexual(coverage_mat_filtered_norm,l_intervals_filtered)
	
	# Dosis estimation, stratified by chromosome type
	logger.info("Starting the calculus of dosis from depth of coverage in autosomes\n")
					
	if l_intervals_filtered_autosomes_index <> []:
		
		l_intervals_filtered_autosomes = map(lambda i: l_intervals_filtered[i] , l_intervals_filtered_autosomes_index)
		mat_all_ratios_autosomes = calculate_dosis(coverage_mat_filtered_norm[l_intervals_filtered_autosomes_index,:],l_intervals_filtered_autosomes,l_samples_filtered,os.path.join(tmp_output_path,"dosis_ratios_autosomes.mat"),logger)

	mat_all_ratios_chrX_f = mat_all_ratios_chrX_m = mat_all_ratios_chrY_m = None
		
	l_intervals_filtered_chrX_intv = []
	l_intervals_filtered_chrY_intv = []
	
	if not options.no_sexual:
		
		logger.info("Starting the calculus of dosis from depth of coverage in chrX\n")
						
		if l_samples_filtered_f <> [] and l_intervals_filtered_chrX_index <> []:
			l_intervals_filtered_chrX_intv = map(lambda i: l_intervals_filtered[i] , l_intervals_filtered_chrX_index)
			mat_all_ratios_chrX_f = calculate_dosis(coverage_mat_filtered_norm[l_intervals_filtered_chrX_index,:][:,l_female_index],l_intervals_filtered_chrX_intv,l_samples_filtered_f,os.path.join(tmp_output_path,"dosis_ratios_chrX_females.mat"),logger)
		if l_samples_filtered_m <> [] and l_intervals_filtered_chrX_index <> []:
			l_intervals_filtered_chrX_intv = map(lambda i: l_intervals_filtered[i] , l_intervals_filtered_chrX_index)
			mat_all_ratios_chrX_m = calculate_dosis(coverage_mat_filtered_norm[l_intervals_filtered_chrX_index,:][:,l_male_index],l_intervals_filtered_chrX_intv,l_samples_filtered_m,os.path.join(tmp_output_path,"dosis_ratios_chrX_males.mat"),logger)
		
		logger.info("Starting the calculus of dosis from depth of coverage in chrY\n")
		
		if l_samples_filtered_m <> [] and l_intervals_filtered_chrY_index <> []:
			l_intervals_filtered_chrY_intv = map(lambda i: l_intervals_filtered[i], l_intervals_filtered_chrY_index)
			mat_all_ratios_chrY_m = calculate_dosis(coverage_mat_filtered_norm[l_intervals_filtered_chrY_index,:][:,l_male_index],l_intervals_filtered_chrY_intv,l_samples_filtered_m,os.path.join(tmp_output_path,"dosis_ratios_chrY_males.mat"),logger)
	
	# Calculus of significance level in autosomes chrm
	logger.info("Establishing the intervals with ratios that are significant in autosomes chromosomes\n")
			
	hash_intervals_autosomes,mat_zscores_autosomes = establish_deletions_and_duplications_zscore(mat_all_ratios_autosomes,l_intervals_filtered_autosomes,l_samples_filtered,THRESHOLD_DEL,THRESHOLD_DUP,SIG_LEVEL_DEL,SIG_LEVEL_DUP,logger)
			
	# Calculus of significance level in sexual chrm
	mat_zscores_chrX_f = mat_zscores_chrX_m = mat_zscores_chrY_m = None
	
	hash_intervals_chrX_females = {}
	hash_intervals_chrX_males   = {}
	hash_intervals_chrY_males   = {}
	
	if not options.no_sexual:
		logger.info("Establishing the intervals with ratios that are significant in sexual chromosomes\n")
				
		logger.info("Establishing in females the intervals with ratios that are significant in sexual chromosomes (chrX)\n")
		if l_samples_filtered_f <> [] and l_intervals_filtered_chrX_index <> [] and len(l_samples_filtered_f) > 1:		
			hash_intervals_chrX_females,mat_zscores_chrX_f = establish_deletions_and_duplications_zscore(mat_all_ratios_chrX_f,l_intervals_filtered_chrX_intv,l_samples_filtered_f,THRESHOLD_DEL,THRESHOLD_DUP,SIG_LEVEL_DEL,SIG_LEVEL_DUP,logger)
	
		logger.info("Establishing in males the intervals with ratios that are significant in sexual chromosomes (chrX and chrY)\n")
		if l_samples_filtered_m <> [] and l_intervals_filtered_chrX_index <> [] and len(l_samples_filtered_m) > 1:
			hash_intervals_chrX_males,mat_zscores_chrX_m = establish_deletions_and_duplications_zscore(mat_all_ratios_chrX_m,l_intervals_filtered_chrX_intv,l_samples_filtered_m,THRESHOLD_DEL,THRESHOLD_DUP,SIG_LEVEL_DEL,SIG_LEVEL_DUP,logger)
		if l_samples_filtered_m <> [] and l_intervals_filtered_chrY_index <> [] and len(l_samples_filtered_m) > 1:
			hash_intervals_chrY_males,mat_zscores_chrY_m = establish_deletions_and_duplications_zscore(mat_all_ratios_chrY_m,l_intervals_filtered_chrY_intv,l_samples_filtered_m,THRESHOLD_DEL,THRESHOLD_DUP,SIG_LEVEL_DEL,SIG_LEVEL_DUP,logger)
	
	# Saving the complete matrices
	l_matrices = []
	
	numpy.save(os.path.join(tmp_output_path,"dosis_ratios_autosomes.mat"),mat_all_ratios_autosomes)
	l_matrices.append(os.path.join(tmp_output_path,"dosis_ratios_autosomes.mat"))
	
	if type(mat_all_ratios_chrX_f) is numpy.ndarray:
		numpy.save(os.path.join(tmp_output_path,"dosis_ratios_chrX_females.mat"),mat_all_ratios_chrX_f)
		l_matrices.append(os.path.join(tmp_output_path,"dosis_ratios_chrX_females.mat"))
	if type(mat_all_ratios_chrX_m) is numpy.ndarray:
		numpy.save(os.path.join(tmp_output_path,"dosis_ratios_chrX_males.mat"),mat_all_ratios_chrX_m)
		l_matrices.append(os.path.join(tmp_output_path,"dosis_ratios_chrX_males.mat"))
	if type(mat_all_ratios_chrY_m) is numpy.ndarray:
		numpy.save(os.path.join(tmp_output_path,"dosis_ratios_chrY_males.mat"),mat_all_ratios_chrY_m)
		l_matrices.append(os.path.join(tmp_output_path,"dosis_ratios_chrY_males.mat"))
	
	numpy.save(os.path.join(tmp_output_path,"zscore_dosis_autosomes.mat"),mat_zscores_autosomes)
	l_matrices.append(os.path.join(tmp_output_path,"zscore_dosis_autosomes.mat"))
	
	if type(mat_zscores_chrX_f) is numpy.ndarray:
		numpy.save(os.path.join(tmp_output_path,"zscore_dosis_chrX_females.mat"),mat_zscores_chrX_f)
		l_matrices.append(os.path.join(tmp_output_path,"zscore_dosis_chrX_females.mat"))
	if type(mat_zscores_chrX_m) is numpy.ndarray:
		numpy.save(os.path.join(tmp_output_path,"zscore_dosis_chrX_males.mat"),mat_zscores_chrX_m)
		l_matrices.append(os.path.join(tmp_output_path,"zscore_dosis_chrX_males.mat"))
	if type(mat_zscores_chrY_m) is numpy.ndarray:
		numpy.save(os.path.join(tmp_output_path,"zscore_dosis_chrY_males.mat"),mat_zscores_chrY_m)
		l_matrices.append(os.path.join(tmp_output_path,"zscore_dosis_chrY_males.mat"))
	
	# Saving the ratios and zscores in all intervals
	for i_s,sample in enumerate(l_samples_filtered):
	
		bed_all_filename = os.path.join(cnv_output_path,"%s_all_intervals.bed" % (sample))
		
		fo = open(bed_all_filename,'w')
		
		if type(mat_all_ratios_autosomes) is numpy.ndarray:
	
			mat_sample = mat_all_ratios_autosomes[:,[i_s]]
			mat_sample_zscore = mat_zscores_autosomes[:,[i_s]]
		
			for j_int,intv in enumerate(l_intervals_filtered_autosomes):
			
				dosis  = mat_sample[j_int]
				zscore = mat_sample_zscore[j_int]
			
				fo.write("%s\t%s\t%s\t%s\t%1.3f\t%1.3f\n" % (intv[0],intv[1],intv[2],intv[3],dosis,zscore))
		
		fo.close()
		
		l_final_files.append(bed_all_filename)
	
	if type(mat_all_ratios_chrX_f) is numpy.ndarray:
				
		for if_s,sample in enumerate(l_samples_filtered_f):
				
			mat_sample_f = mat_all_ratios_chrX_f[:,[if_s]]
			if type(mat_zscores_chrX_f) is numpy.ndarray:
				mat_sample_zscore_f = mat_zscores_chrX_f[:,[if_s]]
			
			bed_all_filename = os.path.join(cnv_output_path,"%s_all_intervals.bed" % (sample))
			
			fo = open(bed_all_filename,'a')
				
			for jf_int,intv in enumerate(l_intervals_filtered_chrX_intv):
				
				dosis  = mat_sample_f[jf_int]
				
				if type(mat_zscores_chrX_f) is numpy.ndarray:
					zscore = "%1.3f" % (mat_sample_zscore_f[jf_int])
				else:
					zscore = '-'
				
				fo.write("%s\t%s\t%s\t%s\t%1.3f\t%s\n" % (intv[0],intv[1],intv[2],intv[3],dosis,zscore))
				
			fo.close()
		
	if type(mat_all_ratios_chrX_m) is numpy.ndarray:
		
		for im_s,sample in enumerate(l_samples_filtered_m):
			
			mat_sample_m_chrX = mat_all_ratios_chrX_m[:,[im_s]]
			if type(mat_zscores_chrX_m) is numpy.ndarray:
				mat_sample_zscore_m_chrX = mat_zscores_chrX_m[:,[im_s]]
			
			bed_all_filename = os.path.join(cnv_output_path,"%s_all_intervals.bed" % (sample))
			
			fo = open(bed_all_filename,'a')
			
			for jf_int,intv in enumerate(l_intervals_filtered_chrX_intv):
			
				dosis = "%1.3f" % (mat_sample_m_chrX[jf_int])
				
				if type(mat_zscores_chrX_m) is numpy.ndarray:
					zscore = "%1.3f" % (mat_sample_zscore_m_chrX[jf_int])
				else:
					zscore = "-"
			
				fo.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (intv[0],intv[1],intv[2],intv[3],dosis,zscore))
			
			fo.close()
	
	if type(mat_all_ratios_chrY_m) is numpy.ndarray:

		for im_s,sample in enumerate(l_samples_filtered_m):
			
			mat_sample_m_chrY = mat_all_ratios_chrY_m[:,[im_s]]
			if type(mat_zscores_chrY_m) is numpy.ndarray:
				mat_sample_zscore_m_chrY = mat_zscores_chrY_m[:,[im_s]]
			
			bed_all_filename = os.path.join(cnv_output_path,"%s_all_intervals.bed" % (sample))
			
			fo = open(bed_all_filename,'a')
			for jf_int,intv in enumerate(l_intervals_filtered_chrY_intv):
			
				dosis = "%1.3f" % (mat_sample_m_chrY[jf_int])
				if type(mat_zscores_chrY_m) is numpy.ndarray:
					zscore = "%1.3f" % (mat_sample_zscore_m_chrY[jf_int])
				else:
					zscore = "-"
			
				fo.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (intv[0],intv[1],intv[2],intv[3],dosis,zscore))
			
			fo.close()
			
	# merge of intervals
	logger.info("Merging intervals into samples\n")
	hash_intervals = {}
	
	for sample in l_samples_filtered:

		hash_intervals.setdefault(sample,[])
		
		hash_intervals[sample].extend(hash_intervals_autosomes.get(sample,[]))
		hash_intervals[sample].extend(hash_intervals_chrX_females.get(sample,[]))
		hash_intervals[sample].extend(hash_intervals_chrX_males.get(sample,[]))
		hash_intervals[sample].extend(hash_intervals_chrY_males.get(sample,[]))
	
	# overlapping counting
	logger.info("Enriching each intervals with overlapping counts")
	
	hash_intervals_annot = calculate_overlapping(hash_intervals)

	for sample in l_samples:

		if sample not in l_samples_filtered:
			
			fo = open(os.path.join(cnv_output_path,"%s_all_intervals.bed" % (sample)),'w')
			fo.close()
			
			fo2 = open(os.path.join(cnv_output_path,"%s.cnv.bed" % (sample)),'w')
			fo2.close()
			
			cnv_report_filename = write_reports_filtered(sample,cnv_output_path,analysis_bed_tagged)
			
			l_final_files.append(cnv_report_filename)
			
			continue
		
		logger.info("Writing the del/dup intervals of sample: %s" % (sample))
		
		#l_intv_sig = filter(lambda (chrom,s,e,annot,sample,rate,zscore,sig,type_,other,c_less,s_less,c_gain,s_gain): sig == '*', sorted(hash_intervals_annot[sample],key=lambda x: int(x[0][2:]) if x[0][2:].isdigit() else x[0]))
		
		l_intv = sorted(hash_intervals_annot[sample],key=lambda x: __sort_key(x))
		
		l_intv_sig = filter(lambda (chrom,s,e,annot,annot2,annot3,sample,rate,zscore,sig,type_,other,c_less,s_less,c_gain,s_gain): sig == '*', l_intv)
		
		l_intv_sig_str = map(lambda (chrom,s,e,annot,annot2,annot3,sample,rate,zscore,sig,type_,other,c_less,s_less,c_gain,s_gain): "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%1.2f\t%s\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s" % (chrom,s,e,annot,annot2,annot3,sample,rate,zscore,sig,type_,__sum_overlap(c_less,c_gain),c_less,s_less,c_gain,s_gain,other), l_intv_sig)
		
		bed_sig_filename = os.path.join(cnv_output_path,"%s.cnv.bed" % (sample))
		
		fo2 = open(bed_sig_filename,'w')
		fo2.write('\n'.join(l_intv_sig_str))
		fo2.close()
		
		l_final_files.append(bed_sig_filename)
					
		cnv_report_filename = write_reports(sample,cnv_output_path,analysis_bed_tagged)
		
		l_final_files.append(cnv_report_filename)
		
	exists_final_all = True
	
	l_final_files_error = []
	
	for filename_final in l_final_files:
		if not os.path.exists(filename_final):
			exists_final_all = False
			l_final_files_error.append(filename_final)
		
	if not exists_final_all:
		fo_e = open(os.path.join(cnv_output_path,"error.txt"),"w")
		fo_e.write("\n".join(l_final_files_error))
		fo_e.close()
			
		logger.info("ERROR: Not all final files exists")
	else:
		fo_ok = open(os.path.join(cnv_output_path,"final_files.txt"),"w")
		fo_ok.write("\n".join(sorted(l_final_files)))
		fo_ok.close()
		logger.info("Checked that all final files exists. Written \'final_files.txt\' file")
		
	if options.clean:
		shutil.rmtree(tmp_output_path)
		logger.info("Cleaned all temporal files")			
						
	
############################################################################333

if __name__=='__main__':
	
	try:
		run()
	except Exception as e: 
		print(e)
	
