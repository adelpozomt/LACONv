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

def run(argv=None):
    
    if argv is None: argv = sys.argv    
   
    parser = OptionParser(add_help_option=True,description="The script annotate bed analysis file and create a well-formed one")
    
    parser.add_option("--b",default=None,help="Bed analysis file of the gene panel/exome",dest="input_bedfile")
    parser.add_option("--env",default=None,help="Bed file with the annotations of the intervals",dest="input_envfile")
    parser.add_option("--o",default=None,help="Bed file with the annotations of the intervals",dest="output_path")
    
    
    slope = 10
                        
    # Se leen las opciones aportadas por el usuario
    (options, args) = parser.parse_args(argv[1:])

    if len(argv) == 1:
        sys.exit(0)
    
    if not parser.check_required("--b"):
        raise IOError('preoprocess_analysis_bed: The bed input file has not been provided')

    if not parser.check_required("--env"):
        raise IOError('preoprocess_analysis_bed: The annotation input bed file has not been provided')
    
    if not parser.check_required("--o"):
        raise IOError('preoprocess_analysis_bed: The output path has not been provided')
        
    input_bedfile = options.input_bedfile

    if not os.path.exists(input_bedfile):
        raise IOError('preoprocess_analysis_bed: The input bed file does not exist.\n%s' % (input_bedfile))
    
    input_envfile = options.input_envfile
    
    if not os.path.exists(input_envfile):
        raise IOError('preoprocess_analysis_bed: The input annotation bed file does not exist.\n%s' % (input_envfile))
    
    output_path = options.output_path
    
    if not os.path.exists(output_path):
        raise IOError('preoprocess_analysis_bed: The output path does not exist.\n%s' % (output_path))
    
    if input_envfile[-2:] == 'gz':
    
        with gzip.open(input_envfile, 'rb') as f:
            l_ = f.readlines()
            f.close()
    
        bed_env = BedTool(l_)
        
    else:
        
        bed_env = BedTool(input_envfile)
    
    l_bed = []
    
    for intv in bed_env.intersect(BedTool(input_bedfile)):
        
        if intv[4].find('intron') <> -1:
            continue
        
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


