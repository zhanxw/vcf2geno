#!/net/fantasia/home/zhanxw/python27/bin/python
import sys, os, re
try:
    from XiaoweiLib import myopen
except:
    sys.path = [re.sub(r'^/home/zhanxw/', '/net/fantasia/home/zhanxw/', x) for x in sys.path]
    sys.path.append('/net/nfsb/fantasia/home/zhanxw/mylib/Python/')
    from XiaoweiLib import myopen

def usage():
    print("%s -o outputPrefix seqFile mapFile" % sys.argv[0] )
    print(" print coverage per sample and coverage per marker")

# return (min, median, mean, max, sd)
import numpy
def getStat(arr):
    length = len(arr)
    if length == 0:
	return (0, 0, 0, 0, 0, 0, 0, 0)
    
    x = numpy.sort(numpy.array(arr))
    x_min = x[0]
    x_q1 = x[int(length/4)]
    x_median = x[int(length/2)]
    x_mean = numpy.mean(x)
    x_q3 = x[int(length*3/4)]
    x_max =  x[length-1]
    x_sd = numpy.std(x)
    return (length, x_min, x_q1, x_median, x_mean, x_q3, x_max, x_sd)
    
if __name__ == '__main__':
    try:
        import getopt
	optlist, args = getopt.getopt(sys.argv[1:], 'o:')
	optlist = dict(optlist)
        seqFile, mapFile = args
	outputPrefix = optlist['-o']
    except:
        usage()
	raise
        sys.exit(1)

    print >> sys.stderr, "Calculate per sample coverage ..."
    data = []  # coverage matrix, sample by marker
    samples = []
    for lineNo, ln in enumerate(myopen(seqFile)):
    	fd = ln.strip().split('\t')
    	coverage = [ int(i.split(' ')[0]) for i in fd[6:] ]
    	data.append(coverage)
    	samples.append(fd[1])

    sampleCov = open(outputPrefix + '.cov.sample', 'w')
    sampleCov.write('PersonId')
    sampleCov.write('\t')
    sampleCov.write('\t'.join(['length', 'min', 'q1', 'median', 'mean', 'q3', 'max', 'sd']))
    sampleCov.write('\n')
    for idx, s in enumerate(samples):
    	sampleCov.write('%s\t' % samples[idx])
    	sampleCov.write('%s\n' % '\t'.join(['%.3f' % i for i in getStat(data[idx])]))
    sampleCov.close()

    print >> sys.stderr, "Calculate per marker coverage ..."
    markers = [ln.strip() for ln in open(mapFile).xreadlines()]
    # data = [ [] for i in xrange(len(markers))]
    # for lineNo, ln in enumerate(myopen(seqFile)):
    # 	fd = ln.strip().split('\t')
    # 	coverage = [int(i.split(' ')[0]) for i in fd[6:]]
    # 	for idx, cov in enumerate(coverage):
    # 	    if cov > 0:
    # 		try:
    # 		    data[idx].append(cov)
    # 		except:
    # 		    import pdb
    # 		    pdb.set_trace()
    # 		    raise

    markerCov = open(outputPrefix + '.cov.marker', 'w')
    markerCov.write('\t'.join(['chrom', 'markerName', 'cm', 'pos', 'ref']))
    markerCov.write('\t')
    markerCov.write('\t'.join(['length', 'min', 'q1', 'median', 'mean', 'q3', 'max', 'sd']))
    markerCov.write('\n')
    for idx, m in enumerate(markers):
	markerCov.write('%s\t' % markers[idx])
	markerCov.write('%s\n' % '\t'.join(['%.3f' % i for i in getStat([s[idx] for s in data])]))
    markerCov.close()

