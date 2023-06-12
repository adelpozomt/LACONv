'''
Created on January 2020

@author: adelpozo
'''

#!/usr/bin/python

import sys, re, shlex , os, string, urllib, time, math, random, subprocess, shutil

import optparse

import gzip

from pybedtools import BedTool

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

def get_exons_gtf(gtf_filename):
    
    l_genes = []
    
    bed_gtf = BedTool(gtf_filename)
    
    for intv in bed_gtf:
        
        feature = intv[2]
        
        if feature == "exon":
            
            chrom = intv[0]
            start = intv[3]
            end   = intv[4]
            strand = intv[6]
            
            attrib_field = intv[8][:-1].replace("\"","").strip()
        
            hash_attrib = dict(list(map(lambda x: x.strip().split(' '), attrib_field.split(';'))))
            
            #chr1    1167648 1168658 B3GALT6 exon_1  NM_080605.3
            gene = hash_attrib.get('gene_name','.')
            transcript_id = hash_attrib.get('transcript_id','.')
            exon_name = hash_attrib.get('exon_number','.')
            
            if exon_name <> '.':
                exon_name = "exon_%s" % (exon_name)
            
            l_genes.append((chrom,start,end,gene,exon_name,transcript_id))
            
    return BedTool(l_genes)

#######################################################################

def run(argv=None):
    
    if argv is None: argv = sys.argv    
   
    parser = OptionParser(add_help_option=True,description="The script annotate bed analysis file and create a well-formed one")
    
    parser.add_option("--b",default=None,help="Bed analysis file of the gene panel/exome",dest="input_bedfile")
    parser.add_option("--gtf",default=None,help="Gtf file with the annotations of the intervals",dest="input_gtffile")
    parser.add_option("--o",default=None,help="Output path where the formatted file must be written",dest="output_path")
        
    slope = 10
                        
    # Se leen las opciones aportadas por el usuario
    (options, args) = parser.parse_args(argv[1:])

    if len(argv) == 1:
        sys.exit(0)
    
    if not parser.check_required("--b"):
        raise IOError('preoprocess_analysis_bed: The bed input file has not been provided')

    if not parser.check_required("--gtf"):
        raise IOError('preoprocess_analysis_bed: The annotation input gtf file has not been provided')
    
    if not parser.check_required("--o"):
        raise IOError('preoprocess_analysis_bed: The output path has not been provided')
        
    input_bedfile = options.input_bedfile

    if not os.path.exists(input_bedfile):
        raise IOError('preoprocess_analysis_bed: The input bed file does not exist.\n%s' % (input_bedfile))
    
    input_gtffile = options.input_gtffile
    
    if not os.path.exists(input_gtffile):
        raise IOError('preoprocess_analysis_bed: The input gtf file does not exist.\n%s' % (input_gtffile))
    
    output_path = options.output_path
    
    if not os.path.exists(output_path):
        raise IOError('preoprocess_analysis_bed: The output path does not exist.\n%s' % (output_path))
    
    bed_gtf = get_exons_gtf(input_gtffile)
         
    l_bed = []
    
    for intv in bed_gtf.intersect(BedTool(input_bedfile)):
        
        chrom = intv[0]
        start = min(int(intv[1]),int(intv[2]))-slope
        end   = max(int(intv[1]),int(intv[2]))+slope
        annot = "%s|%s|%s" % (intv[3],intv[4],intv[6]) 
        
        l_bed.append((chrom,start,end,annot))
        
    output_filename = os.path.join(output_path,os.path.splitext(os.path.basename(input_bedfile))[0]+'.annot_s%d.bed' % (slope))
        
    bed_preproc = BedTool(l_bed).sort().merge(c=4,o='collapse',delim=';')
        
    bed_preproc.saveas(output_filename)
    
    print "Written the bed file: %s" % (output_filename)    
    
    
############################################################################333

if __name__=='__main__':
    
    run()


