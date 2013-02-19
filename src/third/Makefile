# download all third party library and compile them

LIB = tabix pcre eigen gsl bzip2
all: $(LIB)

tabix: tabix-0.2.6.tar.bz2
	-mv tabix-0.2.6.tar.bz2* tabix-0.2.6.tar.bz2
	tar jvxf tabix-0.2.6.tar.bz2
	ln -s -f tabix-0.2.6 tabix
	sed -i 's:-D_USE_KNETFILE::' tabix/Makefile
	sed -i 's:-D_USE_KNTFILE::' tabix/Makefile
	(cd tabix-0.2.6; make; make libtabix.so.1)

bzip2: bzip2-1.0.6.tar.gz
	tar zvxf $<
	ln -s -f bzip2-1.0.6 bzip2
	(cd bzip2-1.0.6; make)

######################################################################
# download sources
bzip2-1.0.6.tar.gz:
	echo "obtain Bzip2 library"
	wget http://www.bzip.org/1.0.6/bzip2-1.0.6.tar.gz
clean:
	-rm -rf $(LIB)
deepclean:
	-rm -rf $(LIB)
	-rm -rf tabix-0.2.6
	-rm -rf pcre-8.21
	-rm -rf gsl-1.15
	-rm -rf eigen-3.1.0
	-rm -rf bzip2-1.0.6

