#!/net/fantasia/home/zhanxw/python27/bin/python
import sys, os, re
try:
    from XiaoweiLib import myopen
except:
    sys.path = [re.sub(r'^/home/zhanxw/', '/net/fantasia/home/zhanxw/', x) for x in sys.path]
    sys.path.append('/net/nfsb/fantasia/home/zhanxw/mylib/Python/')
    from XiaoweiLib import myopen

def usage():
    print("%s [-b bedFileToExcludeRegion] [-i idFile] -m mapaFile -o outputPrefix pileupFiles...." % sys.argv[0] )
    print("bedfile: format:chrom, beg, end, rsNumber")
    print("idFile: pileupNumber, GWAS_ID - optional")
    # print("mapaFile: map file with additional ref and alt column")
    # print("outputPrefix: output .coverage file and .refCount file")

# return counts (ref, alt)
def count(ref, reads):
    if len(ref) != 1:
        return (0, 0)
    # stripping ^ and $
    reads = re.sub(r'\^.', '', reads)
    reads = reads.replace('$', '')
    reads = reads.replace('$', '')
    
    # stripping insertion and deletion
    while reads.find('+') > 0:
        pos = reads.find('+')
        end = pos + 1
        while reads[end].isdigit():
            end += 1
        reads = reads[:pos] + reads[end:]
    while reads.find('-') > 0:
        pos = reads.find('-')
        end = pos + 1
        while reads[end].isdigit():
            end += 1
        reads = reads[:pos] + reads[end:]
    reads = reads.replace('*','')

    # count ref and alt 
    r, a = 0, 0
    bases = set(list('acgtnACGTN'))
    for i in reads:
        if i == '.' or i == ',':
            r += 1
            continue
        if i in bases:
            a += 1
            continue
        print >> sys.stderr, "Unrecognized base", reads
        break
    return (r, a)

def printHeader(colDict, mapRef, fout):
    row1= ['.' for i in xrange(6)]
    row2= ['.' for i in xrange(6)]
    for pos, rsId in colDict.iteritems():
        ref = mapRef[rsId]
        row1.append(rsId)
        row2.append(ref)
    fout.write('\t'.join(row1))
    fout.write('\n')
    fout.write('\t'.join(row2))
    fout.write('\n')       
    
if __name__ == '__main__':
    try:
        import getopt
	optlist, args = getopt.getopt(sys.argv[1:], 'b:i:m:o:')
	optlist = dict(optlist)
	arg_bedFile = optlist.get('-b', None)
	arg_idFile = optlist.get('-i', None)
        arg_mapFile = optlist['-m']
        arg_outPrefix = optlist['-o']
        fns = args
    except:
        usage()
	raise
        sys.exit(1)
    
    from XiaoweiLib import BedFile
    bedFile = BedFile()
    if arg_bedFile != None:
	print >> sys.stderr, "Load bed lines and unique regions:", bedFile.open(arg_bedFile, trimChrPrefix = True)
    else:
	print >> sys.stderr, "Skip loading bed file, all loci will be processed."
#     print bedFile.contain("1", 100) #196341101       196341250
#     print bedFile.contain("1", 196341101) #196341101       196341250
#     print bedFile.contain("1", 196341249)
#     print bedFile.contain("1", 196341250)
#     print bedFile.contain("1", 196341251)
#     print '----------'
#     print bedFile.contain("22", 33412592)
#     print bedFile.contain("22", 33412740)
#     print bedFile.contain("22", 33412741)
#     print bedFile.contain("22", 33412742)
# # 33412592        33412741
#     sys.exit(0)

    mapContent = [i.strip().split() for i in open(arg_mapFile).xreadlines()]
    pos = ['%s:%s' % ( i[0].replace('chr',''), i[1]) for i in mapContent]
    rsId = [ i[2] for i in mapContent]
    # print >> sys.stderr, "Open bed file", bedFile
    # bedContent = [i.split() for i in open(bedFile).xreadlines()]
    #pos = [ '%s:%s' % (i[0].replace('chr',''), i[2]) for i in bedContent]
    #rsId = [ i[3] for i in bedContent]
    from collections import OrderedDict
    colDict = OrderedDict(zip(pos, rsId))  # pos -> rsId e.g. 1:123->rs99
    excludePos = [bedFile.contain( i[0].replace('chr',''), i[1]) for i in mapContent]
    excludePos = set([pos[idx] for idx, v in enumerate(excludePos) if v])
    print >> sys.stderr, "Excluding %d on-site markers" % len(excludePos)

    print >> sys.stderr, "%d markers loaded" % len(colDict)
    print >> sys.stderr, "%d pileup files" % len(fns)
    #mapRef = dict([i.strip().split() for i in open(mapRefFile).xreadlines()]) # rsId -> ref
    #mapRef = dict([ ])
    #print >> sys.stderr, "%d refbases loaded" % len(mapRef)

    if arg_idFile != None:
        idContent = [i.strip().split() for i in open(arg_idFile).xreadlines()] # pileupId -> GWAS_ID
	pileupId = {i[0]:i[1:] for i in idContent}
        print >> sys.stderr, "%d sample id loaded" % len(pileupId)

    # covFile = open(outPrefix + '.coverage', 'w')
    # refCountFile = open(outPrefix + '.refCount', 'w')
    # printHeader(colDict, mapRef, covFile)
    # printHeader(colDict, mapRef, refCountFile)
    logFile = open(arg_outPrefix + '.log', 'w')

    seqFile  = open(arg_outPrefix + '.seq', 'w')
    for fn in fns:
        res = dict()
        for ln in myopen(fn):  # loop each pileup file
            try:
                chrom, pos, ref, depth, reads, quals = ln.strip().split()
            except:
                print >> sys.stderr, "File [ %s ] is empty or not valid" % fn
                break
            #print 'reads = ', reads,
            try:
                refCount, altCount = count(ref, reads)
            except:
		print >> sys.stderr, "Cannot parse pileup data, entering debug mode ..."
                import pdb
                pdb.set_trace()
            #print 'refCnt = ', refCount, 'altCnt=', altCount
            key = '%s:%s' % (chrom.replace('chr', ''), pos )
            if key not in colDict:
                continue
            res[key] = (refCount, altCount)

        # outputs
        if arg_idFile != None:
	    key = os.path.basename(fn).split('.')[0]
	    if key in pileupId:
		gwasId = pileupId[key]
		if len(gwasId) == 2:
		    seqFile.write('\t'.join(gwasId))
		elif len(gwasId) == 1:
		    seqFile.write('\t'.join([gwasId, gwasId]))
		else:
		    print >> sys.stderr, "Error: id file does not have correct entry for ", key
		    seqFile.write('\t'.join([gwasId, gwasId]))
        else:
            gwasId = os.path.basename(fn).replace('.pileup', '')
	    seqFile.write('\t'.join([gwasId, gwasId]))

        for k, v in colDict.iteritems():
            if k in excludePos:
                seqFile.write('\t0 0')
                logFile.write('Exclude\t%s\n' % k)
                continue
            if k in res:
                ref, alt = res[k]
                seqFile.write('\t%d %d' % (ref+alt, ref))
            else:
                #"no pile up at position"
                seqFile.write('\t0 0')
        seqFile.write('\n')

    seqFile.close()

    # mapFile  = open(arg_outPrefix + '.map', 'w')
    # for fd in mapContent:
    #     chrom, pos, rs, ref, alt = fd
    #     mapFile.write("%s\t%s\t%s\t%s\t%s\n" % (chrom, rs, mapDist, pos, ref))

    # mapFile.close()
    logFile.close()
