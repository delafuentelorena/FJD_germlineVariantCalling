#!/usr/bin/env python


import sys
import os
from glob import glob
import subprocess
import argparse
from collections import Counter
import datetime


def joinFastq(basespaceF, fastqFolder, inputF, sampleFile, analysis):

	if basespaceF=="True":

		### Mount basespace with basemount
		sys.stdout.write("\nMounting Basespace...\n")
		basespaceDir='/mnt/genetica3/lorena/basespaceBioinfo' 
		subprocess.call(["basemount",basespaceDir])

		sys.stdout.write("\nMounting Basespace COMPLETE...\n")

		if not os.path.isdir(basespaceDir+'/Projects/'+inputF): 
			sys.stderr.write("ERROR: Project '%s' does not exist in basespace\n" %(inputF))
			## unmount basespace
			sys.exit()
		else:
			path = glob(basespaceDir+'/Projects/'+inputF+'/Samples/*/Files/*.fastq.gz')
			#print(path)
	
	else: 
		if not os.path.isdir(inputF): 
			sys.stderr.write("ERROR: Input folder '%s' does not exist\n" %(inputF))
			sys.exit()
		else:
			inputF=os.path.realpath(inputF)
			path = glob(inputF+'/*.fastq.gz')


	# checking and reading sample names in sample file (if provided)

	file_samples = list()
	if sampleFile != "all":
		if os.path.isfile(sampleFile):
			samplesFile = open(sampleFile, "r")
			[file_samples.append((line.split(",")[1]).strip()) for line in samplesFile]
		else:
			sys.stderr.write("ERROR: '%s' does not exist\n" %(sampleFile))
			sys.exit()


	# sample concatenation 

	print file_samples

	sample_names = list()
	for sample in path:
		
		sys.stdout.write("\nANALYSING SAMPLE %s\n" %(sample))

		file = os.path.basename(sample) 
		dnaid = file[0:file.find('_')]

		if dnaid not in sample_names and (analysis in ["cnv","all"] or sampleFile=="all" or dnaid in file_samples):
			
			sample_names.append(dnaid)	

			if basespaceF=="True":
				sys.stdout.write("\ni am here\n")
				foward = glob(basespaceDir+'/Projects/'+inputF+'/Samples/*/Files/' + dnaid + '*R1*.fastq.gz')
				reverse = glob(basespaceDir+'/Projects/'+inputF+'/Samples/*/Files/' + dnaid + '*R2*.fastq.gz')
				sys.stdout.write("\n%s\n" %(":".join(foward)))

			else:
				foward = glob(inputF+ '/' + dnaid + '*R1*.fastq.gz')
				reverse = glob(inputF+ '/' + dnaid + '*R2*.fastq.gz')
			
			if len(foward)==len(reverse) and len(foward)!= 0:
				sys.stdout.write("\n- SAMPLE: "+ dnaid+"\n")
				
				### FASTQ CONCATENATION	
				sys.stdout.write("\ni am here2222\n")

				foward_file = fastqFolder + dnaid + '_R1.fastq.gz'
				reverse_file = fastqFolder + dnaid + '_R2.fastq.gz'
				sys.stdout.write("\n%s\n" %(foward_file))

				sys.stdout.write(foward_file+"\n")
				sys.stdout.write(reverse_file+"\n")

				sys.stdout.write("Joined fastqs will be stored in '%s' folder\n" %(fastqFolder))

				sys.stdout.write('Joining forward fastq files of sample ' + dnaid +"\n")
				mycmd = 'cat %s > %s' %(' '.join(foward), foward_file)
				subprocess.call(mycmd, shell=True)

				sys.stdout.write('Joining reverse fastq files of sample ' + dnaid +"\n")
				mycmd2 = 'cat %s > %s' %(' '.join(reverse), reverse_file)
				subprocess.call(mycmd2, shell=True)

			else:
				sys.stderr.write("ERROR: Not fastq files found for sample '%s' or different number of reverse and foward fastq files.\n" %(sample_name))


	if basespaceF=="True":

		### Mount basespace with basemount
		sys.stdout.write("\nUnmounting Basespace...\n")
		subprocess.call(["basemount","--unmount", basespaceDir])



if __name__ == "__main__":
    
    import sys
    joinFastq(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4],sys.argv[5])


