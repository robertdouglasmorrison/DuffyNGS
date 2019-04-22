#
# DuffyNGS  ver 1.8.0

VER = 1.8.0

DuffyNGS.tar.gz : 
	rm  -f ./DuffyNGS_*.tar.gz

	# make sure we compile on the current system
	rm  -f src/*.o
	rm  -f src/samtools/*.o

	${R_PROGRAM} CMD SHLIB --output=src/DuffyNGS.so src/samtools/*.c  src/*.c 
	${R_PROGRAM} CMD build  --force .
	${R_PROGRAM} CMD INSTALL ./DuffyNGS_${VER}.tar.gz

