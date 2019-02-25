#!/usr/bin/env python


if args.basespace:
	### Mount basespace with basemount
	sys.stdout.write("\nMounting Basespace...\n")
	basespaceDir='/mnt/genetica3/lorena/basespaceBioinfo' 
	mycmdbs = 'basemount %s' %(basespaceDir)
	subprocess.call(["basemount",basespaceDir])
	sys.stdout.write("\nestoy aqui...\n")

	if not os.path.isdir(basespaceDir+'/Projects/'+args.input): 
		sys.stderr.write("ERROR: Project '%s' does not exist in basespace\n" %(args.input))
		sys.exit()
	else:
		path = glob(basespaceDir+'/Projects/'+args.input+'/Samples/*/Files/*.fastq.gz')

else:
	if not os.path.isdir(args.input): 
		sys.stderr.write("ERROR: Input folder '%s' does not exist\n" %(args.input))
		sys.exit()
	else:
		args.input=os.path.realpath(args.input)
		path = glob(args.input+'/*.fastq.gz')
